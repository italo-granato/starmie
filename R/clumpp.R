# clumpp.R
# Functions for merging Q matrices from Structure runs.

#' Run the CLUMPP algorithms.
#' @param Q_list A list of of Q matrices.
#' @param method The algorithm to use to infer the correct permutations. One of 'greedy' or 'greedyLargeK'
#' @import iterpc
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
      permuations <- iter_wrapper(iterpc(ncol(Q_list[[i+1]]), ordered = TRUE))

      foreach(x=permuations, .combine=c) %do% J_perm(x, Q_list[[i+1]], Q_list[1:i])

      max_perm <- 1:ncol(Q_list[[i+1]])
      max <- c(-Inf)
      j <- 0
      repeat {
        value <- J_perm(getnext(permuations), Q_list[[i+1]], Q_list[1:i])
        if (value>max){
          max <- value
          max_perm <- getcurrent(permuations)
        }
        j <- j+1
      }
      Q_list[[i+1]] <- Q_list[[i+1]][,max_perm]
    }
  } else{
    #code for LargeKGreedy algorithm
    for (i in 1:(length(Q_list-1))){
      apply(combn(1:ncol(Q_list[[i]]), 2), 2, function(x) { c(G(Q_list[[i]][,x[1]], Q_list[[i+1]][,x[2]]),x[1],x[2])})

    }

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

