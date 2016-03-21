#' Read Structure Output
#' @param  path directory containing output of STRUCTURE runs
#' @param  n_runs number of runs for each K value
#' @param logfile include STRUCTURE logs for each run (required if wanting to do MCMC plots)
#' @export
readStructure <- function(path, n_runs, logfile = FALSE) {
  # i/o checks
  if(!is.character(path)) {
    stop("Path must be character variable.")
  }
  if(!dir.exists(path)) {
    stop("Path does not exist.")
  }

  if(!is.finite(n_runs) && !is.integer(n_runs)) {
    stop("n_runs must be a finite integer")
  }

  # do the work

}
