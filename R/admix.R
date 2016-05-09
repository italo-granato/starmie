# admix.R
#' Constructor for admix object
#'
#' @return an admix object which is a list with 6 elements:
#'  K: number of clusters estimated by ADMIXTURE
#'  nsamples: number of samples used
#'  nmarkers: number of markers used
#'  Q_df: a data.frame of cluster membership probabilities
#'  P_df: a data.frame of estimated marker frequencies in each inferred population
#'  log_info: a data.frame containing the K, CVerror and logLik of the last model.
#' @export
admix <- function() {
  structure(list(K = NULL, nsamples = NULL, nmarkers = NULL,
                 Q_df =  NULL, P_df = NULL, log_info = NULL),
            class = "admix")
}

