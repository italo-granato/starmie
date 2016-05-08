# clumpp.R
# Functions for merging Q matrices from Structure runs.

#' Run the CLUMPP algorithms.
#' @param Q_list A list of of Q matrices.
#' @param method The algorithm to use to infer the correct permutations. One of 'greedy' or 'greedyLargeK'
#' @importFrom combinat permn
#' @importFrom purrr map_dbl
#' @export
#' @examples
#' # Read in Structure files
#' structure_files <- system.file("extdata/microsat_testfiles", package="starmie")
#' structure_output_files <- list.files(structure_files, pattern = ".*K10.*out_f", full.names = TRUE)
#' structure_log_files <- list.files(structure_files, pattern = ".*K10.*log", full.names = TRUE)
#' structure_runs <- mapply(loadStructure, structure_output_files, structure_log_files, SIMPLIFY=FALSE)
#' Q_list <- lapply(structure_runs, getQ)
#' clumpp(Q_list)

clumpp <- function(Q_list, method="greedy"){

  if (method=="greedy"){
    #Greedy clumpp algorithm
    for (i in 1:(length(Q_list)-1)){
      permuations <- permn(1:ncol(Q_list[[i+1]]))
      perm_scores <- map_dbl(permuations, J_perm, Q_list[[i+1]], Q_list[1:i])
      Q_list[[i+1]] <- Q_list[[i+1]][,permuations[[which.max(perm_scores)]]]
    }
  } else{
    #code for LargeKGreedy algorithm
  }

  return(Q_list)
}

G <- function(Q_1, Q_2){
  W <- matrix(1, nrow(Q_1), ncol(Q_1))/ncol(Q_1)
  1-norm(Q_1-Q_2, type="F")/sqrt(norm(Q_1-W, type="F")*norm(Q_2-W, type="F"))
}

J_perm <- function(perm, Q_x, Q_sub_list){
  J(Q_x[,perm], Q_sub_list)
}

J <- function(Q_x, Q_sub_list){
  sum(unlist(lapply(Q_sub_list, G, Q_x)))/length(Q_sub_list)
}

