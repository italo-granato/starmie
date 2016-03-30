# loadAdmixture.R
# Functions for parsing output from admixture


#' Read Admixture Output
#' @param  starmie_obj an object of class \code{\link{starmie}} with sample_data filled
#' @param  file_list a vector of .Q file names for admixture
#' @param  logfile_list a vector of logfile names for each admixture run (required for diagnositc plots)
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#' @importFrom data.table fread
#' @export
#' @examples
#' my_starmie <- starmie()
#' # add in sample metadata
#' sample_file <- system.file("extdata/hapmap3_files", "hapmap3.fam", package="starmie")
#' my_starmie <- loadSampleData(my_starmie, sample_file)
#' # read in Q-files
#' hapmap_files <- system.file("extdata/hapmap3_files", package="starmie")
#' qfiles <- list.files(hapmap_files, pattern = ".Q$", full.names = TRUE)
#' my_starmie <- loadAdmixture(my_starmie, qfiles)
#'
loadAdmixture <- function(starmie_obj, file_list, logfile_list = NULL) {
  # i/o checks
  if (!inherits(starmie_obj, "starmie")) stop("Not a valid starmie object.")
  if (!all(is.character(file_list) & !is.na(file_list))) stop("file_list must be a character vector.")
  if (!is.null(logfile_list)) {
    if (!all(is.character(logfile_list) & !is.na(logfile_list))) stop("logfile_list must be a character vector.")
  }
  # check whether sample id is available, and sample meta data is inputted
  # otherwise throw an errror.
  if (is.null(starmie_obj$sample_data)) {
    stop("Sample metadata is required.")
    if(is.null(starmie_obj$sample_data$sample.id)) {
      stop("Sample identifier required")
    }
  }
  # do the work
  # at the moment just read in q_files, number of Q-files corresponds to K
  read_qfiles <- function(qfile) {
    q_df <- fread(qfile, data.table=FALSE, header=FALSE)
    stopifnot(nrow(q_df) == nrow(starmie_obj$sample_data))
    q_df$K <- rep(ncol(q_df), nrow(q_df))
    q_df$sample.id <- starmie_obj$sample_data$sample.id
    gather(q_df, cluster, probability, -sample.id, -K)
  }

  allQ <- bind_rows(lapply(file_list, read_qfiles))
  # change cluster character name into integer
  allQ$cluster <- as.integer(gsub("V", "", allQ$cluster))

  # read in logfiles
  if (!is.null(logfile_list)) {
    log_info <- bind_rows(lapply(logfile_list, read_logfiles))
  } else {
    log_info <- NULL
  }

  # output list containing q_info, and log_info
  starmie_obj$admixture_run <- list(q_info = allQ, log_info = log_info)
  starmie_obj

}

read_logfiles <- function(logfile) {
  # grab final logliklihood and CV error
  logfin <- readLines(logfile)
  loglik <- as.numeric(gsub("Loglikelihood: ", "",
                            logfin[grepl("^Loglikelihood:", logfin)]))
  cverror <- logfin[grepl("^CV error", logfin)]
  K <- as.integer(gsub(".*\\(K=(.*)\\).*", "\\1", cverror))
  cverror <- as.numeric(gsub(".*: ", "\\1", cverror))
  data.frame(K = K, logL = loglik, CVerror = cverror, stringsAsFactors = FALSE)
}
