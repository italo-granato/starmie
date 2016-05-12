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
#' clump <- clumpp(Q_list)
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
    #Use LargeKGreedy algorithm
    Q_list <- largeKGreedy(Q_list)
  }
  return(Q_list)
}

memoryGreedy <- function(Q_list){
  #Faster but with a high memory footprint for large K
  for (i in 1:(length(Q_list)-1)){
    permuations <- permn(1:ncol(Q_list[[i+1]]))
    perm_scores <- map_dbl(permuations, J_perm, Q_list[[i+1]], Q_list[1:i])
    Q_list[[i+1]] <- Q_list[[i+1]][,permuations[[which.max(perm_scores)]]]
  }
  return(Q_list)
}

iterativeGreedy <- function(Q_list){
  #slower but memory efficient
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

largeKGreedy <- function(Q_list){
  #Initial iteration
  #calculate pairwise column comparisons
  column_pairs <- getall(iterpc(ncol(Q_list[[1]]), 2, ordered = TRUE, replace = TRUE))
  pair_comparisons <- apply(column_pairs, 1
                            , function(x) { list(G = G(Q_list[[1]][,x[1],drop=FALSE], Q_list[[2]][,x[2],drop=FALSE])
                                                 , Qy = x[1]
                                                 , Qz = x[2])})
  pair_comparisons_df <- data.frame(matrix(unlist(pair_comparisons), ncol=3, byrow=TRUE))
  #permute the second Q matrix
  Q_list[[2]] <- Q_list[[2]][,get_best_permutation(pair_comparisons_df)]

  if (length(Q_list)>2){
    for (i in 3:length(Q_list)){
      #remaining iterations
      pair_comparisons <- apply(column_pairs, 1
                                , function(x) { list(J = J_largeK(Q_list[1:(i-1)], Q_list[[i]], x[1], x[2])
                                                     , Qy = x[1]
                                                     , Qz = x[2])})
      pair_comparisons_df <- data.frame(matrix(unlist(pair_comparisons), ncol=3, byrow=TRUE))
      Q_list[[i]] <- Q_list[[i]][,get_best_permutation(pair_comparisons_df)]
    }
  }

  return(Q_list)
}

G <- function(Q_1, Q_2){
  W <- matrix(1, nrow(Q_1), ncol(Q_1))/ncol(Q_1)
  1-norm(Q_1-Q_2, type="F")/sqrt(norm(Q_1-W, type="F")*norm(Q_2-W, type="F"))
}

get_best_permutation <- function(pair_df){
  #returns the best permutation of columns based on a dataframe of pairwise comparisons
  K <- max(pair_df[,2])
  best_pairs_df <- data.frame(x=rep(0.0, K), y=rep(0,K), z=rep(0,K))
  for (i in 1:K){
    best_pairs_df[i,] <- pair_df[which.max(pair_df[,1]),]
    pair_df <- pair_df[pair_df[,2]!=pair_df[which.max(pair_df[,1]),2], ]
  }
  best_pairs_df <- best_pairs_df[order(best_pairs_df[,2]),]
  return(best_pairs_df[,3])
}

J_largeK <- function(Q_sub_list, Q_x, y, z){
  1/length(Q_sub_list) * sum(
    unlist(lapply(Q_sub_list, function(Q){
      G(Q[,y,drop=FALSE], Q_x[,z,drop=FALSE])
      }))
    )
}

J_perm <- function(perm, Q_x, Q_sub_list){
  J(Q_x[,perm], Q_sub_list)
}

J <- function(Q_x, Q_sub_list){
  sum(unlist(map_dbl(Q_sub_list, G, Q_x)))/length(Q_sub_list)
}

