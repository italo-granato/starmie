# bestK.R
# Functions for determining a suitable K value from multiple Structure runs.


#' Determine a suitable K value from multiple Structure runs
#' @param  structure_runs an object of class \code{\link{list}} with structure runs represented as \code{\link{struct}} objects.
#' @param  method the method used to calculate the best K
#' @param  plot whether of not to generate diagnostic plots
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' # read in Structure files
#' structure_files <- system.file("extdata/microsat_testfiles", package="starmie")
#' structure_output_files <- list.files(structure_files, pattern = "*.out_f", full.names = TRUE)
#' structure_log_files <- list.files(structure_files, pattern = ".*log", full.names = TRUE)
#' structure_runs <- mapply(loadStructure, structure_output_files, structure_log_files, SIMPLIFY=FALSE)
#'
bestK <- function(structure_runs, method="evano", plot=TRUE){
  #i/o checks
  if (!all(unlist(lapply(structure_runs, inherits, "struct")) & !is.na(structure_runs))) stop("structure_runs must be a vector of struct objects.")
  if (!(method %in% c("evano", "structure"))) stop("method must be one of 'evano' or 'structure'")
  if (!is.logical(plot)) stop("plot must be one of (TRUE/FALSE)")

  #order structure runs and collect by K
  posterior_probs = data.frame(
    K=unlist(lapply(structure_runs, getK)),
    pos_prob=unlist(lapply(structure_runs, getPosterior)))


  #If K are not sequential or do not have equal number of runs just perform the structure approach and produce a simpler plot
  Ks <- table(posterior_probs$K)

  if (!is.sequential(names(Ks)) |
      !all(Ks[[1]]==Ks) |
      method=="structure"){

    if (method!="structure") warning("WARNING! K values are not sequential or there are an uneven number of runs per K. Reverting to structure method.")

    #summarise data by K
    posterior_probs_summary <- group_by(posterior_probs, K) %>%
      summarize(
        runs = n(),
        meanK = mean(pos_prob),
        sdK = sd(pos_prob)
      )

    #generate plot
    gg <- ggplot(posterior_probs_summary, aes(x=K, y=meanK))
    gg <- gg + geom_point()
    gg <- gg + geom_errorbar(aes(ymax=meanK+sdK, ymin=meanK-sdK), width=0.1)
    gg <- gg + theme_bw()
    gg <- gg + scale_x_continuous(breaks = min(posterior_probs_summary$K):max(posterior_probs_summary$K))
    gg <- gg + ylab("Mean log posterior probability")

    if (plot) {suppressWarnings(print(gg))}

    return(posterior_probs_summary$K[which.max(posterior_probs_summary$meanK)])
  }

  #split by K and calculate the first and second derivatives
  split_posterior_probs <- split(posterior_probs, posterior_probs$K)
  split_posterior_probs <- split_posterior_probs[order(as.numeric(names(split_posterior_probs)))]

  prev <- 0
  posterior_probs_summary <- data.frame()
  for (k in 2:(length(split_posterior_probs)-1)){
    dK <- split_posterior_probs[[k]]$pos_prob - split_posterior_probs[[k-1]]$pos_prob
    ddK <- split_posterior_probs[[k+1]]$pos_prob-2*split_posterior_probs[[k]]$pos_prob+split_posterior_probs[[k-1]]$pos_prob
    LK <-  split_posterior_probs[[k]]$pos_prob

    if (sd(LK)==0) stop("No deviation between runs of the same K. Evano statistics cannot be computed.")

    K <- as.numeric(names(split_posterior_probs)[k])
    temp_df <- data.frame(K = K,
                          variable = "L(K)",
                          value = mean(LK),
                          sd = sd(LK))
    posterior_probs_summary <- rbind(posterior_probs_summary, temp_df)

    temp_df <- data.frame(K = K,
                          variable = "L'(K)",
                          value = mean(dK),
                          sd = sd(dK))
    posterior_probs_summary <- rbind(posterior_probs_summary, temp_df)

    temp_df <- data.frame(K = K,
                          variable = "L''(K)",
                          value = mean(abs(ddK)),
                          sd = sd(ddK))
    posterior_probs_summary <- rbind(posterior_probs_summary, temp_df)

    temp_df <- data.frame(K = K,
                          variable = "delta K",
                          value = mean(abs(ddK))/sd(LK),
                          sd = NA)
    posterior_probs_summary <- rbind(posterior_probs_summary, temp_df)

  }


  gg <- ggplot(posterior_probs_summary, aes(x=K, y=value))
  gg <- gg + geom_point()
  gg <- gg + geom_errorbar(aes(ymax=value+sd, ymin=value-sd))
  gg <- gg + facet_wrap(~variable, ncol=2, scales = "free_y")
  gg <- gg + theme_bw()
  gg <- gg + scale_x_continuous(breaks = min(posterior_probs_summary$K):max(posterior_probs_summary$K))

  if (plot) {suppressWarnings(print(gg))}

  best <- posterior_probs_summary[posterior_probs_summary$variable=="delta K",]
  best <- best$K[which.max(best$value)]

  return(best)

}

is.sequential <- function(x){
  x <- as.numeric(x)
  all(diff(x) == diff(x)[1])
}
