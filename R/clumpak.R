# clumpak.R
# Function for inferring modes withing multiple structure runs.

#' Run the CLUMPP algorithms.
#' @param Q_list A list of of Q matrices.
#' @param method The algorithm to use to infer the correct permutations. One of 'greedy' or 'greedyLargeK'
#' @import iterpc
#' @importFrom combinat permn
#' @importFrom purrr map_dbl
#' @importFrom MCL mcl
#' @importFrom proxy simil
#' @export
#' @examples
#' # Read in Structure files
#' structure_files <- system.file("extdata/microsat_testfiles", package="starmie")
#' structure_output_files <- list.files(structure_files, pattern = ".*K10.*out_f", full.names = TRUE)
#' structure_log_files <- list.files(structure_files, pattern = ".*K10.*log", full.names = TRUE)
#' structure_runs <- mapply(loadStructure, structure_output_files, structure_log_files, SIMPLIFY=FALSE)
#' Q_list <- lapply(structure_runs, getQ)
#' clumpak_results <- clumpak(Q_list)
#'

clumpak <- function(Q_list){
  simMatrix <- as.matrix(simil(Q_list, method=G))
  t <- calcThreshold(simMatrix)
  simMatrix[simMatrix<t] <- 0
  clusters <- mcl(simMatrix, addLoops = TRUE)$Cluster
  split(Q_list, clusters)
  }

calcThreshold <- function(simMatrix){
  min_edge <- min(simMatrix, na.rm = TRUE)
  max_edge <- max(simMatrix, na.rm = TRUE)

  thresholds <- seq(min_edge,max_edge,length.out=100)
  is_valid <- unlist(lapply(thresholds, function(t){
    degrees <- rowSums(simMatrix>t, na.rm=TRUE)
    (sum(degrees>0)/length(degrees)<0.1) & (mean(degrees) > 0.5*nrow(simMatrix))
  }))
  threshold <- thresholds[c(1:length(is_valid))[is_valid][-1]]
  return(threshold)
}

G <- function(Q_1, Q_2){
  W <- matrix(1, nrow(Q_1), ncol(Q_1))/ncol(Q_1)
  1-norm(Q_1-Q_2, type="F")/sqrt(norm(Q_1-W, type="F")*norm(Q_2-W, type="F"))
}
