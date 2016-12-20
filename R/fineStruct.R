# fineStruct.R
#' Constructor for fineStruct object
#'
#' @return an fineStruct object which is a list with 3 elements:
#'  nsamples: number of samples used
#'
#'  cfactor: estimate of the c parameter in the fineStructure model
#'
#'  chunkcounts_df: a data.frame of chunkcounts
#'
#'  dendro: a dendrogram object estimated by fineStructure
#'
#'  mcmc_df: a data.frame containing the mcmc output of fineStructure
#'
#' @export
fineStruct <- function() {
  structure(list(nsamples = NULL, chunkcounts_df = NULL, cfactor = NULL,
                 tree =  NULL, mcmc_df = NULL),
            class = "fineStruct")
}

#' @export
print.fineStruct <- function(x, ...) {
  cat("fineStruct object containing run information\n")
  cat("Model parameters:\n")
  cat(paste("  ", "No. samples:", x$nsamples, "\n"))
  cat(paste("  ", "No. markers:", x$nmarkers, "\n"))
}
