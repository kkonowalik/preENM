# Batch cropping of raster files ---------------------------------------------

#' Crop and mask many raster files to M (or any polygon), safely and on disk
#'
#' Cuts every raster file in a set of folders (e.g. all CHELSA or
#' CHELSA-TraCE21k time slices) to the accessible area or another polygon and
#' writes the results to a new folder. Designed for large collections: the
#' work is done by GDAL directly from file to file, so memory use stays low
#' however large the rasters are, finished files are skipped when the run is
#' repeated, and one failing file does not stop the others.
#'
#' **Values are copied exactly.** Many climate datasets (CHELSA among them)
#' store integers with a scale and offset in the file metadata (e.g.
#' temperature x 10 in Kelvin). Cropping by reading the values into R and
#' writing them back in the original integer type can silently round them
#' and drop the scale and offset (we found errors of up to 0.95 deg C with
#' `terra::crop(..., datatype = datatype(r))`). Here the raw values, the data
#' type, the scale, offset and NoData value are kept, and the output grid is
#' aligned exactly to the input grid, so no cell is resampled.
#'
#' A cell is kept when its centre lies inside the polygon (`all_touched =
#' FALSE`, the usual rule) or when the polygon touches it at all
#' (`all_touched = TRUE`). Cells outside are set to NoData; if a file has no
#' NoData value, the extreme value of its data type is used (e.g. 65535 for
#' unsigned 16-bit integers, NaN for floating point).
#'
#' @param files Raster files and/or folders to search.
#' @param m The area to cut to: an `m_area` from [delimit_m()], polygons
#'   (`sf`, [terra::SpatVector]) or a path (for a GeoPackage written by
#'   [write_m_area()] its `M` layer is used).
#' @param out_dir Output folder (created if needed). Output files keep their
#'   names; it must not be an input folder.
#' @param variables Optional variable names to select, matched as whole
#'   parts of the file name between `_`, `-` or `.` (e.g. `c("bio01",
#'   "bio12", "scd")` selects `CHELSA_TraCE21k_bio01_-100_V1.0.tif` but not
#'   `..._bio010_...`). Case-insensitive.
#' @param pattern Regular expression for the files searched in folders.
#' @param recursive Search sub-folders.
#' @param mask Set cells outside the polygon to NoData (default). With
#'   `FALSE` the rasters are only cropped to the polygon's extent.
#' @param all_touched Keep every cell touched by the polygon, not only cells
#'   whose centre is inside.
#' @param overwrite Re-create output files that already exist (default:
#'   skip them, so an interrupted run can simply be repeated).
#' @param compress GeoTIFF compression (`"LZW"`, `"DEFLATE"`, `"ZSTD"` or
#'   `"NONE"`).
#' @param verbose Print progress.
#' @return Invisibly, a data.frame with `input`, `output` and `status`
#'   (`written`, `skipped (exists)`, `no overlap` or `failed: ...`).
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' src <- file.path(tempdir(), "layers")
#' dir.create(src, showWarnings = FALSE)
#' terra::writeRaster(sim$env[[1]], file.path(src, "CHELSA_bio01_1981.tif"),
#'                    overwrite = TRUE)
#' terra::writeRaster(sim$env[[2]], file.path(src, "CHELSA_bio12_1981.tif"),
#'                    overwrite = TRUE)
#' res <- crop_rasters(src, m, file.path(tempdir(), "layers_M"),
#'                     variables = "bio01")
#' res
crop_rasters <- function(files, m, out_dir, variables = NULL,
                         pattern = "\\.tif{1,2}$", recursive = TRUE,
                         mask = TRUE, all_touched = FALSE, overwrite = FALSE,
                         compress = "LZW", verbose = TRUE) {
  say <- function(...) if (verbose) message(...)
  is_dir <- dir.exists(files)
  missing <- files[!is_dir & !file.exists(files)]
  if (length(missing)) {
    stop("Not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  found <- c(files[!is_dir],
             unlist(lapply(files[is_dir], list.files, pattern = pattern,
                           full.names = TRUE, recursive = recursive,
                           ignore.case = TRUE)))
  found <- unique(found)
  if (!is.null(variables)) {
    rx <- paste0("(^|[_.-])(", paste(variables, collapse = "|"),
                 ")([_.-]|$)")
    found <- found[grepl(rx, basename(found), ignore.case = TRUE)]
  }
  if (!length(found)) stop("No raster files selected.", call. = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_norm <- normalizePath(out_dir)
  if (any(normalizePath(dirname(found)) == out_norm)) {
    stop("`out_dir` must not be a folder containing the input files.",
         call. = FALSE)
  }
  dup <- duplicated(basename(found))
  if (any(dup)) {
    stop("Several input files share a name (e.g. ",
         basename(found[dup][1]), "); they would overwrite each other in ",
         "`out_dir`.", call. = FALSE)
  }

  if (is.character(m) && length(m) == 1 && grepl("\\.gpkg$", m) &&
      "M" %in% sf::st_layers(m)$name) {
    m <- sf::st_read(m, layer = "M", quiet = TRUE)
  }
  poly <- if (inherits(m, "m_area")) m$m else as_sfc(m, "m")
  poly <- sf::st_make_valid(sf::st_union(poly))

  tmp <- tempfile("preENM_crop")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  cutlines <- list()                       # one cutline file per CRS
  cutline_for <- function(r) {
    key <- terra::crs(r)
    if (is.null(cutlines[[key]])) {
      f <- file.path(tmp, paste0("cut", length(cutlines) + 1, ".gpkg"))
      g <- sf::st_transform(poly, sf::st_crs(key))
      sf::st_write(sf::st_sf(id = 1, geometry = g), f, quiet = TRUE)
      cutlines[[key]] <<- list(file = f, ext = terra::ext(terra::vect(g)))
    }
    cutlines[[key]]
  }

  n <- length(found)
  status <- character(n)
  out_files <- file.path(out_dir, basename(found))
  say("Cropping ", n, " file(s) to ", out_dir)
  for (i in seq_len(n)) {
    f <- found[i]
    out <- out_files[i]
    if (file.exists(out) && !overwrite) {
      status[i] <- "skipped (exists)"
      say(sprintf("[%d/%d] skipped (exists): %s", i, n, basename(f)))
      next
    }
    status[i] <- tryCatch({
      r <- suppressWarnings(terra::rast(f))   # errors are still caught
      cl <- cutline_for(r)
      e <- terra::intersect(cl$ext, terra::ext(r))
      if (is.null(e)) {
        "no overlap"
      } else {
        e <- terra::align(e, r, snap = "out")
        opts <- c("-te", num_arg(c(terra::xmin(e), terra::ymin(e),
                                   terra::xmax(e), terra::ymax(e))),
                  "-tr", num_arg(terra::res(r)),
                  "-r", "near", "-overwrite",
                  "-co", paste0("COMPRESS=", compress), "-co", "TILED=YES",
                  "-co", "BIGTIFF=IF_SAFER")
        if (mask) {
          opts <- c(opts, "-cutline", cl$file)
          if (all_touched) opts <- c(opts, "-wo", "CUTLINE_ALL_TOUCHED=TRUE")
          if (is.na(source_nodata(f))) {
            opts <- c(opts, "-dstnodata", default_nodata(terra::datatype(r)[1]))
          }
        }
        if (file.exists(out)) unlink(out)
        sf::gdal_utils("warp", source = f, destination = out, options = opts,
                       quiet = TRUE)
        if (!file.exists(out)) stop("GDAL did not write the file")
        "written"
      }
    }, error = function(err) {
      if (file.exists(out)) unlink(out)
      paste("failed:", conditionMessage(err))
    })
    say(sprintf("[%d/%d] %s: %s", i, n, status[i], basename(f)))
  }
  res <- data.frame(input = found,
                    output = ifelse(status %in% c("written",
                                                  "skipped (exists)"),
                                    out_files, NA_character_),
                    status = status, stringsAsFactors = FALSE)
  tab <- table(sub(":.*", "", status))
  say("Done: ", paste(names(tab), tab, sep = " = ", collapse = ", "))
  if (any(startsWith(status, "failed"))) {
    warning(sum(startsWith(status, "failed")), " file(s) failed; see the ",
            "`status` column.", call. = FALSE)
  }
  invisible(res)
}

# NoData value of the first band of a raster file, or NA (internal)
source_nodata <- function(f) {
  info <- sf::gdal_utils("info", f, quiet = TRUE)
  m <- regmatches(info, regexpr("NoData Value=[^\n]*", info))
  if (!length(m)) return(NA_character_)
  sub("NoData Value=", "", m)
}

# NoData used when a file has none: the extreme of its data type (internal)
default_nodata <- function(datatype) {
  switch(datatype,
         INT1U = "255", INT1S = "-128", INT2U = "65535", INT2S = "-32768",
         INT4U = "4294967295", INT4S = "-2147483648",
         INT8U = "18446744073709551615", INT8S = "-9223372036854775808",
         "nan")
}

# Numbers as GDAL arguments: full precision and no padding (format() pads a
# vector to a common width, e.g. " 2", which newer GDAL versions reject)
num_arg <- function(x) sprintf("%.15g", x)
