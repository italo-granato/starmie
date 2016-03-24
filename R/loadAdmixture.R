# loadAdmixture.R
# Functions for parsing output from admixture


#' Read Admixture Output
#' @param  path directory containing admixture run information
#' @param  k_values integer vector of K values tried
#' @param  logfile include logfiles for each admixture run (required for diagnositc plots)
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#' @importFrom data.table fread
#' @export
loadAdmixture <- function(starmie_obj, path, k_values, logfile = FALSE) {
  # i/o checks
  if (!inherits(starmie_obj, "starmie")) stop("Not a valid starmie object.")
  if (!is.character(path) ) stop("Path must be character variable.")
  if ( !dir.exists(path) ) stop("Path does not exist.")
  if ( !is.logical(logfile) ) stop("logfile must be TRUE/FALSE")
  if(any(!is.finite(k_values) | !is.integer(k_values) | is.na(k_values))) {
    stop("k_values must all be finite integers")
  }
  # check whether sample id is available, and sample meta data is inputted
  # otherwise throw an errror.
  if (is.null(starmie_obj$sample_data)) {
    stop("Sample metadata is required.")
  }

  # do the work
  # at the moment just read in q_files
  q_files <- list.files(path, ".Q$", full.names = TRUE)
  read_qfiles <- function(qfile) {
    q_df <- fread(qfile, data.table=FALSE, header=FALSE)
    q_df$sample.id <- 1:nrow(q_df)
    q_df$K <- rep(ncol(q_df), nrow(df))
    tidyr::gather(q_df, cluster, probability, -sample.id, -K)
  }

  allQ <- dplyr::bind_rows(lapply(q_files, read_qfiles))

  starmie_obj$admixture_run <- allQ
  starmie_obj

}
