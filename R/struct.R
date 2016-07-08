# starmie.R
# Place documentation for test datasets here.

#' Constructor for struct object
#' @description struct object for storing structure run information
#' @details The \code{structure} object is a list with 11 elements:
#'  K: number of clusters estimated by Structure
#'  run_params: the run paramters given to the Structure program
#'  mem_df: assigned cluster mebership proportions
#'  alle_freqs: Estimated allele frequencies
#'  avg_dist_df: Average distance between individuals
#'  fit_stats_df: Model fit statistics
#'  fst_df: Fst values
#'  ancest_df: Inferred ancestory of individuals
#'  clust_allele_list: Cluster allele frequencies
#'  burn_df: Burn in MCMC iteration output
#'  nonburn_df: Main MCMC iteration output
#' @return a \code{structure} object
#' @export
struct <- function() {
  structure(list(K = NULL, run_params = NULL, mem_df = NULL,
                 alle_freqs = NULL, avg_dist_df = NULL,
                 fst_df=NULL, fit_stats_df = NULL,  ancest_df = NULL,
                 clust_allele_list = NULL,
                 burn_df=NULL, nonburn_df=NULL),
            class = "struct")
}

#' Accessor methods for structure object
#' @description Return K from a \code{\link{struct}}
#' @param structure_obj a code{\link{struct}} object
#' @export
getK <- function(structure_obj) {
  structure_obj$K
}

#' Accessor methods for structure object
#' @describeIn getK  Return the log posterior probability from a \code{\link{struct}}
#' @param structure_obj a code{\link{struct}} object
#' @export
getPosterior <- function(structure_obj){
  structure_obj$fit_stats_df[structure_obj$fit_stats_df$Statistic=="Estimated Ln Prob of Data",][2]
}

#' Accessor methods for structure object
#' @describeIn getQ Return the Q matrix from a \code{\link{struct}}
#' @param structure_obj a code{\link{struct}} object
#' @export
getQ <- function(structure_obj){
  Q <- data.matrix(structure_obj$ancest_df[,4:ncol(structure_obj$ancest_df)])
  rownames(Q) <- structure_obj$ancest_df$Label

  if (ncol(Q)==1){
    colnames(Q) <- "Cluster 1"
  }

  return(Q)
}

#' Acessor methods for structure object
#' @describeIn getMCMC Return non-burn in MCMC iterations.
#' @param structure_obj a code{\link{struct}} object
#' @export
getMCMC <- function(structure_obj) {
  mcmc_df <- data.frame(K = getK(structure_obj),
                        structure_obj$nonburn_df[, c("Rep#:", "Alpha", "Ln Like")])
  colnames(mcmc_df) <- c("K", "Iteration", "Alpha", "LogL")
  mcmc_df
}

#' Example structure objects
#' @param example_type a character string either "clumpp" or "mcmc_diagnostics" or "barplot"
#' @description load structure objects for different starmie functions
#' @export
exampleStructure <- function(example_type) {
  structure_files <- system.file("extdata/microsat_testfiles",
                                 package="starmie")

  if (example_type == "clumpp") {
    k3_out <- list.files(structure_files,
                         pattern = ".*K3.*out_f",
                         full.names = TRUE)

    #k3_log <- list.files(structure_files,
    #                     pattern = ".*K3.*log",
    #                     full.names = TRUE)
    lapply(k3_out, loadStructure)
  } else if (example_type == "mcmc_diagnostics") {
    k10_out <- list.files(structure_files,
                         pattern = ".*K10.*out_f",
                         full.names = TRUE)
    k10_log <- list.files(structure_files,
                         pattern = ".*K10.*log",
                         full.names = TRUE)
    mapply(loadStructure, k10_out, k10_log, SIMPLIFY = FALSE)
  } else if (example_type == "barplot") {
    k6_out <- list.files(structure_files,
                         pattern = "^locprior_K6.out_f",
                         full.names = TRUE)
    loadStructure(k6_out)

  }

}
