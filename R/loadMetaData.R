# loadMetaData.R
# read in sample metadata from a fam or txt file

loadSampleData <- function(starmie_obj, ...)
UseMethod("loadSampleData")
#' Load sample metadata
#'
#' @description Read in sample information from space delimted PLINK fam file
#' or plain text file.
#'
#' @param sample_file filename containing sample information
#' @param pop_identifer logical does file contain geographic location?
#' @details Assume that we have PLINK format FAM file with optional location identifier in the final column
#' @export
loadSampleData <- function(starmie_obj, sample_file, pop_identifier = FALSE) {
  # i/o checks
  if (!is.character(sample_file)) stop("sample_file must be character variable")
  if (!is.logical(pop_identifier)) stop("pop_identifier must be TRUE/FALSE")
  if (!inherits(starmie_obj, "starmie")) stop("Not a valid starmie object")
  # set up connection
  fin <- file(sample_file, "r")
  col_names <- c("family.id", "sample.id", "paternal.id",
                 "maternal.id", "sex", "phenotype")
  sample_data <- read.table(sample_file,
                            header = FALSE,
                            na.strings = c("NA", "-9"),
                            stringsAsFactors = FALSE)
  if (pop_identifier) {
    # check valid FAM file + population_id
    if (ncol(sample_data) != 7) stop("Sample file requires 7 columns to be valid")
    col_names[ncol(sample_data)] <- "population.id"

  } else {
    if(ncol(sample_data) != 6) stop("Sample file requires 6 columns to be valid")
  }
  colnames(sample_data) <- col_names
  starmie_obj$n_samples <- nrow(sample_data)
  starmie_obj$sample_data <- sample_data
  close(fin)
  starmie_obj
}

loadMarkerData <- function(starmie_obj, ...)
UseMethod("loadMarkerData")

#' Read in SNP map file
#'
#' @param starmie_obj object of class \code{\link{starmie}}
#' @param snp_file a PLINK map file to process
#' @param ploidy ploidy of population (default 2 )
#' @export
loadMarkerData <- function(starmie_obj, snp_file, ploidy = 2) {
  # i/o checks
  if (!inherits(starmie_obj, "starmie")) stop("Not a valid starmie object")
  if (!is.character(snp_file)) stop("snp_file must be character variable")
  if (!is.integer(ploidy) && !is.finite(ploidy)) stop("ploidy must be finite integer")
  fin <- file(snp_file, "r")
  col_names <- c("chromosome", "snp.id", "genetic_distance", "position")

  marker_data <- read.table(snp_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(marker_data) <- col_names
  close(fin)
  starmie_obj$n_markers <- nrow(marker_data)
  starmie_obj$ploidy <- ploidy
  starmie_obj$marker_data <- marker_data
  starmie_obj

}
