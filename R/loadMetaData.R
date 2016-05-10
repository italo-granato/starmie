# loadMetaData.R
# read in sample metadata from a fam or txt file

#' Load sample metadata
#'
#' @description Read in sample information from space delimted PLINK fam file
#' or plain text file.
#' @param sample_file filename containing sample information
#' @param pop_identifier logical does file contain population of origin?
#' @details Assume that we have PLINK format FAM file with optional location identifier in the final column
#' @export
loadSampleData <- function(sample_file, pop_identifier = FALSE) {
  # i/o checks
  if (!is.character(sample_file)) stop("sample_file must be character variable")
  if (!is.logical(pop_identifier)) stop("pop_identifier must be TRUE/FALSE")
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
  sample_data

}

#' Read in SNP map file
#'
#' @param snp_file a PLINK map file to process
#' @export
loadMarkerData <- function(snp_file) {
  # i/o checks
  if (!is.character(snp_file)) stop("snp_file must be character variable")
  fin <- file(snp_file, "r")
  col_names <- c("chromosome", "snp.id", "genetic_distance", "position")

  marker_data <- read.table(snp_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(marker_data) <- col_names
  close(fin)
  marker_data

}
