# Taxonomic complexity ----------------------------------------------------------------

#' Taxonomic complexity of species names
#'
#' Records of taxonomically difficult species are more often misidentified.
#' Tessarolo et al. (2021) weighted records by the class of their
#' identifier; when that is unknown, as in most GBIF downloads, the
#' taxonomic complexity of the name is a possible proxy: species with many
#' synonyms or many subspecies and varieties have been split, lumped and
#' re-identified more often. `taxon_complexity()` counts, for each name,
#' its synonyms and its accepted infraspecific taxa in a taxonomic backbone
#' (World Flora Online for plants, or any Darwin Core table such as the
#' GBIF backbone), and turns them into a reliability weight for
#' [record_ignorance()].
#'
#' **Counts.** The name is matched to an accepted taxon of the backbone
#' (by exact name, or with `WorldFlora::WFO.match()` when the WorldFlora
#' package is installed and `fuzzy = TRUE`; a matched synonym is followed
#' to its accepted name). `n_synonyms` is the number of names whose
#' `acceptedNameUsageID` is that taxon and whose status is synonym;
#' `n_infraspecific` the number of accepted subspecies, varieties and forms
#' whose `parentNameUsageID` is that taxon.
#'
#' **Weight.** `reliability = 1 / (1 + log(1 + n))`, with `n` the sum of
#' the counts chosen in `count`: 1 for a name without synonyms or
#' infraspecific taxa, 0.59 with 1, 0.36 with 5, 0.29 with 10, 0.20 with 50.
#' This is a rule of thumb, not a published or validated weighting:
#' synonymy also reflects how much a group has been studied, not only how
#' hard it is to identify. Use it as one option and say so.
#'
#' @param names Species names (e.g. the unique values of the species column
#'   of the records).
#' @param backbone A taxonomic backbone in Darwin Core form, with columns
#'   `taxonID`, `scientificName`, `taxonRank`, `taxonomicStatus`,
#'   `acceptedNameUsageID` and `parentNameUsageID`: e.g. the World Flora
#'   Online backbone (`WorldFlora::WFO.download()`, then
#'   `WorldFlora::WFO.remember()`), or a path to its `classification.csv`.
#' @param count What counts towards complexity: `"synonyms"`,
#'   `"infraspecific"` or both (default).
#' @param fuzzy Match names with `WorldFlora::WFO.match()` (needs the
#'   WorldFlora package) instead of exactly.
#' @return A data.frame with `name`, `matched` (the accepted name, `NA` if
#'   not found), `n_synonyms`, `n_infraspecific` and `reliability` (`NA`
#'   when the name was not found).
#' @references Tessarolo, G., Ladle, R.J., Lobo, J.M., Rangel, T.F.,
#'   Hortal, J. (2021) Using maps of biogeographical ignorance to reveal the
#'   uncertainty in distributional data hidden in species distribution
#'   models. *Ecography* 44, 1743-1755. \doi{10.1111/ecog.05793}
#'
#'   Kindt, R. (2020) WorldFlora: An R package for exact and fuzzy
#'   matching of plant names against the World Flora Online taxonomic
#'   backbone data. *Applications in Plant Sciences* 8, e11388.
#'   \doi{10.1002/aps3.11388}
#' @seealso [record_ignorance()]
#' @export
#' @examples
#' bb <- data.frame(
#'   taxonID = c("t1", "t2", "t3", "t4", "t5", "t6"),
#'   scientificName = c("Aus bus", "Aus cus", "Aus bus var. dus",
#'                      "Xus bus", "Yus bus", "Aus ecus"),
#'   taxonRank = c("species", "species", "variety", "species", "species",
#'                 "species"),
#'   taxonomicStatus = c("Accepted", "Accepted", "Accepted", "Synonym",
#'                       "Synonym", "Accepted"),
#'   acceptedNameUsageID = c("", "", "", "t1", "t1", ""),
#'   parentNameUsageID = c("", "", "t1", "", "", ""))
#' taxon_complexity(c("Aus bus", "Aus cus", "Xus bus"), bb)
taxon_complexity <- function(names, backbone,
                             count = c("synonyms", "infraspecific"),
                             fuzzy = FALSE) {
  count <- match.arg(count, c("synonyms", "infraspecific"),
                     several.ok = TRUE)
  if (is.character(backbone) && length(backbone) == 1) {
    backbone <- utils::read.delim(backbone, quote = "", encoding = "UTF-8",
                                  colClasses = "character")
  }
  bb <- as.data.frame(backbone, stringsAsFactors = FALSE)
  need <- c("taxonID", "scientificName", "taxonRank", "taxonomicStatus",
            "acceptedNameUsageID", "parentNameUsageID")
  miss <- setdiff(need, names(bb))
  if (length(miss)) {
    stop("`backbone` lacks the Darwin Core column(s) ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  bb[need] <- lapply(bb[need], function(v) {
    v <- as.character(v)
    v[is.na(v)] <- ""
    v
  })
  status <- tolower(bb$taxonomicStatus)
  rank <- tolower(bb$taxonRank)
  accepted <- status == "accepted"
  nm <- unique(as.character(names))

  # match names to taxa of the backbone
  id <- if (fuzzy) {
    if (!requireNamespace("WorldFlora", quietly = TRUE)) {
      stop("Fuzzy matching needs WorldFlora: install.packages(",
           "\"WorldFlora\")", call. = FALSE)
    }
    m <- suppressMessages(WorldFlora::WFO.one(WorldFlora::WFO.match(
      spec.data = nm, WFO.data = bb, verbose = FALSE), verbose = FALSE))
    stats::setNames(as.character(m$taxonID), nm)
  } else {
    hit <- match(nm, bb$scientificName)
    # prefer an accepted name when a name occurs several times
    acc_hit <- match(nm, bb$scientificName[accepted])
    out <- ifelse(!is.na(acc_hit), bb$taxonID[accepted][acc_hit],
                  bb$taxonID[hit])
    stats::setNames(out, nm)
  }
  # follow synonyms to their accepted taxon
  row <- match(id, bb$taxonID)
  syn <- !is.na(row) & status[row] == "synonym" &
    nzchar(bb$acceptedNameUsageID[row])
  id[syn] <- bb$acceptedNameUsageID[row[syn]]
  row <- match(id, bb$taxonID)

  is_syn <- status == "synonym"
  is_infra <- accepted & rank %in% c("subspecies", "variety", "form",
                                     "subvariety", "subform", "infraspecies")
  n_syn <- vapply(id, function(t) {
    if (is.na(t)) return(NA_real_)
    sum(is_syn & bb$acceptedNameUsageID == t)
  }, numeric(1))
  n_inf <- vapply(id, function(t) {
    if (is.na(t)) return(NA_real_)
    sum(is_infra & bb$parentNameUsageID == t)
  }, numeric(1))
  n <- (if ("synonyms" %in% count) n_syn else 0) +
    (if ("infraspecific" %in% count) n_inf else 0)
  data.frame(name = nm,
             matched = ifelse(is.na(row), NA_character_,
                              bb$scientificName[row]),
             n_synonyms = unname(n_syn), n_infraspecific = unname(n_inf),
             reliability = unname(1 / (1 + log1p(n))),
             stringsAsFactors = FALSE, row.names = NULL)
}
