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
#' @description Return K from a \code{\link{struct_obj}}
#' @export
getK <- function(structure_obj) {
  structure_obj$K
}

#' Accessor methods for structure object
#' @description Return the log posterior probability from a \code{\link{struct_obj}}
#' @export
getPosterior <- function(structure_obj){
  structure_obj$fit_stats_df[structure_obj$fit_stats_df$Statistic=="Estimated Ln Prob of Data",][2]
}
