# loadAdmixture.R
# Functions for parsing output from admixture

#' Read Admixture Output
#' @param  qfile a valid Q file from ADMIXTURE
#' @param pfile a corresponding P file from ADMXIXTURE
#' @param  logfile logfile from corresponding AMIXTURE run ()
#' @importFrom tidyr gather
#' @importFrom data.table fread
#' @return an \link{admix} object containing the output of of an admixture run
#' @export
#' @examples
#' qfin <- system.file("extdata/hapmap3_files", "hapmap3.2.Q", package = "starmie")
#' pfin <- system.file("extdata/hapmap3_files", "hapmap3.2.P", package = "starmie")
#' my_admix <- loadAdmixture(qfin, pfin)
#' # add log file
#' logfin <- system.file("extdata/hapmap3_files", "log2.out", package = "starmie")
#' my_admix <- loadAdmixture(qfin, pfin, logfin)
#'
loadAdmixture <- function(qfile, pfile, logfile = NULL) {
  # i/o checks
  if ( !(is.character(qfile) & length(qfile) == 1) ) stop("qfile must be a character vector of length 1.")
  if ( !(is.character(pfile) & length(pfile) == 1) ) stop("pfile must be a character vector of length 1.")

  if ( !is.null(logfile) ) {
    if ( !(is.character(logfile) & length(logfile) == 1) ) stop("logfile must be a character vector of length 1.")
  }
  # create new admixture object
  admix_obj <- admix()
  # qfile data reader
  q_df <- fread(qfile, data.table = FALSE, header = FALSE)
  K <- ncol(q_df)
  nsamples <- nrow(q_df)
  colnames(q_df) <- paste("Cluster", seq(1,ncol(q_df)))
  # pfile data reader
  p_df <- fread(pfile, data.table = FALSE, header = FALSE)
  nmarkers <- nrow(p_df)
  if(ncol(p_df) != K) {
    stop("Number of populations does not match between Q and P files")
  }

  # read in logfiles
  if ( !is.null(logfile) ) {
    log_info <- read_logfiles(logfile)
    if(log_info$K != K) {
      stop("logfile does not match input P and Q files")
    }
  } else {
    log_info <- NULL
  }

  # output list containing q_info, p_info and log_info

  admix_obj$K <- K
  admix_obj$nsamples <- nsamples
  admix_obj$nmarkers <- nmarkers
  admix_obj$Q_df <- q_df
  admix_obj$P_df <- p_df
  admix_obj$log_info <- log_info
  return(admix_obj)
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
