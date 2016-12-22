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
                 dendro =  NULL, mcmc_df = NULL, meancoincedence_matrix=NULL,
                 mappopchunk_matrix=NULL),
            class = "fineStruct")
}

#' @export
print.fineStruct <- function(x, ...) {
  cat("fineStruct object containing run information\n")
  cat(paste("  ", "No. samples:", x$nsamples, "\n"))
  cat("Chunk count dataframe:\n")
  print(head(x$chunkcounts_df, n=2))
  cat("Dendrogram:\n")
  print(x$dendro)
  cat("MCMC dataframe:\n")
  print(head(x$mcmc_df, n=2))
  if (!is.null(x$meancoincedence_matrix)) print("Includes mean coincedence matrix")
  if (!is.null(x$mappopchunk_matrix)) print("Includes population-by-population chunkcount matrix")
}
