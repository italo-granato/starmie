# averagePairWiseSimilarityH.R
# Function for computing the average pairwise similarity between Q matrices (H) as defined in Jakobsson et al.

#' Average Q matrices.
#' @param Q_list A list of of Q matrices.
#' @export
#' @examples
#' # Read in Structure files
#' multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
#' Q_list <- lapply(multiple_runs_k10, getQ)
#' avgQ <- averagePairWiseSimilarityH(Q_list)
averagePairWiseSimilarityH <- function(Q_list){
  #i/o checks
  if(!all(unlist(lapply(Q_list, inherits, "matrix")))) stop("cluster runs must be a list of Q matrices")

  R <- length(Q_list)
  H <- 0
  for (i in 1:(R-1)){
    for (j in (i+1):R){
      H <- H + G(Q_list[[i]], Q_list[[j]])
    }
  }
  H <- 2/(R*(R-1)) * H
  return(H)
}


G <- function(Q_1, Q_2){
  W <- matrix(1, nrow(Q_1), ncol(Q_1))/ncol(Q_1)
  1-norm(Q_1-Q_2, type="F")/sqrt(norm(Q_1-W, type="F")*norm(Q_2-W, type="F"))
}
