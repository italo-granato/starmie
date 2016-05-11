# clumpp.R
# Functions for merging Q matrices from Structure runs.

#' Run the CLUMPP algorithms.
#' @param Q_list A list of of Q matrices.
#' @param method The algorithm to use to infer the correct permutations. One of 'greedy' or 'greedyLargeK'
#' @import iterpc
#' @importFrom combinat permn
#' @importFrom purrr map_dbl
#' @export
#' @examples
#' # Read in Structure files
#' structure_files <- system.file("extdata/microsat_testfiles", package="starmie")
#' structure_output_files <- list.files(structure_files, pattern = ".*K3.*out_f", full.names = TRUE)
#' structure_log_files <- list.files(structure_files, pattern = ".*K3.*log", full.names = TRUE)
#' structure_runs <- mapply(loadStructure, structure_output_files, structure_log_files, SIMPLIFY=FALSE)
#' Q_list <- lapply(structure_runs, getQ)
#' clumpp(Q_list)
#'
clumpp <- function(Q_list, method="greedy"){

  if (method=="greedy"){
    #Greedy clumpp algorithm
    K <- ncol(Q_list[[1]])
    if (K>8){
      Q_list <- iterativeGreedy(Q_list)
    }
    else{
      Q_list <- memoryGreedy(Q_list)
    }

  } else{
    #code for LargeKGreedy algorithm
    for (i in 1:(length(Q_list-1))){
      apply(combn(1:ncol(Q_list[[i]]), 2), 2, function(x) { c(G(Q_list[[i]][,x[1]], Q_list[[i+1]][,x[2]]),x[1],x[2])})

    }

  }
  return(Q_list)
}

memoryGreedy <- function(Q_list){
  for (i in 1:(length(Q_list)-1)){
    permuations <- permn(1:ncol(Q_list[[i+1]]))
    perm_scores <- map_dbl(permuations, J_perm, Q_list[[i+1]], Q_list[1:i])
    Q_list[[i+1]] <- Q_list[[i+1]][,permuations[[which.max(perm_scores)]]]
  }
  return(Q_list)
}

iterativeGreedy <- function(Q_list){
  for (i in 1:(length(Q_list)-1)){
    permuations <- iterpc(ncol(Q_list[[i+1]]), ordered = TRUE)

    max_perm <- 1:ncol(Q_list[[i+1]])
    max <- -Inf
    j <- 0
    while(j<getlength(permuations)) {
      value <- J_perm(getnext(permuations), Q_list[[i+1]], Q_list[1:i])
      if (value>max){
        max <- value
        max_perm <- getcurrent(permuations)
      }
      j <- j+1
    }
    Q_list[[i+1]] <- Q_list[[i+1]][,max_perm]
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
  sum(unlist(map_dbl(Q_sub_list, G, Q_x)))/length(Q_sub_list)
}

