# starmie.R
# Place documentation for test datasets here.

#' Constructor for starmie object
#'
#'@description starmie object for storing STRUCTURE/ADMIXTURE run information
#'@details The \code{starmie} object is a list with 7 elements:
#'  n_samples: number of samples in marker set
#'  n_markers: number of SNPs/microsatellites in marker set
#'  ploidy: ploidy of samples
#'  sample_data: data.frame of sample meta information
#'  marker_data: data.frame of SNP/microsatellite meta information.
#'  structure_run: data.frame of STRUCTURE results
#'  adxmiture_run: data.frame of ADMIXTURE results
#'  clumpp: data.frame of CLUMPP results
#'@return a \code{starmie} object
#'@export
starmie <- function() {
 structure(list(n_samples = integer(), n_markers = integer(), ploidy = integer(),
                sample_data = NULL, marker_data = NULL,
                structure_run = NULL,  admixture_run = NULL, clumpp = NULL),
                class = "starmie")
}
