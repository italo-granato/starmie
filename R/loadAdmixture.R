# loadAdmixture.R
# Functions for parsing output from admixture

#' Read Admixture Output
#' @param  path directory containing admixture run information
#' @param  k_values integer vector of K values tried
#' @param  snp_map name of PLINK map file for genomic coordinates
#' @param  logfile include logfiles for each admixture run (required for diagnositc plots)
#' @export
readAdmixture <- function(path, k_values, snp_map, logfile = FALSE) {
  # i/o checks
  if (!is.character(path) ) stop("Path must be character variable.")
  if (!(is.character(snp_map)) ) stop("PLINK map file must be character variable.")
  if ( !dir.exists(path) ) stop("Path does not exist.")

  if(any(!is.finite(k_values) | !is.integer(k_values) | is.na(k_values))) {
    stop("k_values must all be finite integers")
  }
  # do the work

}
