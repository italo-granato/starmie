# clumpak.R
# Function for inferring modes within multiple structure runs.

#' Run the CLUMPP algorithms.
#' @param Q_list A list of of Q matrices.
#' @param method The method the algorithm uses to infer the correct permutations. One of 'greedy' or 'greedyLargeK' or 'stephens' or 'none'
#' @importFrom MCL mcl
#' @importFrom proxy simil
#' @export
#' @examples
#' # Read in Structure files
#' multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
#' Q_list <- lapply(multiple_runs_k10, getQ)
#' clumpak_results <- clumpak(Q_list)
clumpak <- function(Q_list, method="none"){

  # i/o checks
  if (!(method %in% c("greedy", "greedyLargeK", "stephens", "none"))) {
    stop("Not a valid CLUMPP method, please use on of: 'greedy', 'greedyLargeK' or 'stephens' or 'none'")
  }
  if(!all(unlist(lapply(Q_list, inherits, "matrix")))) stop("cluster runs must be a list of Q matrices")

  if (method!="none"){
    Q_list <- clumpp(Q_list, method)$Q_list
  }

  simMatrix <- as.matrix(simil(Q_list, method=G))
  diag(simMatrix) <- 1
  t <- calcThreshold(simMatrix)
  simMatrix[simMatrix<t] <- 0
  clusters <- mcl(simMatrix, addLoops = TRUE)$Cluster
  split(Q_list, clusters)
}

calcThreshold <- function(simMatrix){
  #We want to choose a threshold t such that the number of singletons
  # is less than 10% of nodes and the mean node degree is at least 50%
  # of the total number of nodes.
  thresholds <- as.vector(simMatrix)
  thresholds <- thresholds[order(thresholds)]

  is_valid <- unlist(lapply(thresholds, function(t){
    degrees <- apply(simMatrix, 1, function(r) sum(r[r>t]))
    ((sum(degrees==1)/length(degrees))<0.1) & (mean(degrees) >= 0.5*nrow(simMatrix))
  }))
  threshold <- thresholds[is_valid][sum(is_valid)]
  return(threshold)
}

G <- function(Q_1, Q_2){
  W <- matrix(1, nrow(Q_1), ncol(Q_1))/ncol(Q_1)
  1-norm(Q_1-Q_2, type="F")/sqrt(norm(Q_1-W, type="F")*norm(Q_2-W, type="F"))
}
