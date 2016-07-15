# bestK.R
# Functions for determining a suitable K value from multiple Structure runs.


#' Determine a suitable K value from multiple Structure runs
#' @param  runs a list of \code{\link{struct}} objects or \code{\link{admix}} objects.
#' @param  method the method used to calculate the best K
#' @param  make_plot whether of not to generate diagnostic plots
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' multi_K <- exampleStructure("multiple_runs")
#' # Run the evanno method and display diagnostic plots.
#' bestK(multi_K)
#' bestK(multi_K, "structure")
bestK <- function(runs, method="evanno", make_plot=TRUE){
  #i/o checks
  run_types <- unlist(lapply(runs, "class"))

  if ( !(all(run_types == "struct") | all(run_types == "admix")) )  stop("runs must be a vector of struct or admix objects.")
  if (!(method %in% c("evanno", "structure", "cverror"))) stop("method must be one of 'evanno' or 'structure' or 'cverror'")
  if (!is.logical(make_plot)) stop("make_plot must be one of TRUE or FALSE")

  #order structure runs and collect by K
  if ( all(run_types == "struct") ) {
    message("Creating diagnostic plots for structure runs.")
    posterior_probs = data.frame(
      K=unlist(lapply(runs, getK)),
      pos_prob=unlist(lapply(runs, getPosterior)))


    #If K are not sequential or do not have equal number of runs just
    # perform the structure approach and produce a simpler plot
    Ks <- table(posterior_probs$K)

    if ( !is.sequential(names(Ks)) | !all(Ks[[1]]==Ks) | method == "structure") {

      if (method!="structure") warning("WARNING! K values are not sequential or there are an uneven number of runs per K. Reverting to structure method instead.")
      .bestK_structure(posterior_probs, make_plot)

    } else {
      .bestK_evanno(posterior_probs, make_plot)
    }
  } else if ( all(run_types == "admix") ) {
    message("Creating diagnositc plots for admixture runs")

    log_df <- do.call("rbind", lapply(runs, function(y) y$log_info))
    if ( any(is.na(log_df)) ) {
      stop("Need log file information to produce diagnositc plots")
    }
    .bestK_admixture(log_df, make_plot)

  }

}

.bestK_evanno <- function(posterior_probs, make_plot) {
  # plot delta K
  #split by K and calculate the first and second derivatives

  split_posterior_probs <- split(posterior_probs, posterior_probs$K)
  split_posterior_probs <- split_posterior_probs[order(as.numeric(names(split_posterior_probs)))]
  # empty data frame for storing results
  output_variables <- c("L(K)", "L'(K)",  "L''(K)", "delta K")
  output_values <- list()
  for (k in 2:(length(split_posterior_probs)-1)) {
    dK <- split_posterior_probs[[k]]$pos_prob - split_posterior_probs[[k-1]]$pos_prob
    ddK <- split_posterior_probs[[k+1]]$pos_prob-2*split_posterior_probs[[k]]$pos_prob+split_posterior_probs[[k-1]]$pos_prob
    LK <-  split_posterior_probs[[k]]$pos_prob

    if ( sd(LK) == 0 ) stop("No deviation between runs of the same K. Evanno statistics cannot be computed.")

    output_values[[k]] <-  data.frame(K = as.numeric(names(split_posterior_probs)[k]),
                                                variable = output_variables,
                                                value = c(mean(LK), mean(dK), mean(abs(ddK)), mean(abs(ddK))/sd(LK)),
                                                sd = c(sd(LK), sd(dK), sd(ddK), NA))


  }
  posterior_probs_summary <- do.call("rbind", output_values)

  if (make_plot) {
    gg <- ggplot(posterior_probs_summary, aes(x=K, y=value)) +
      geom_point() +
      geom_errorbar(aes(ymax=value+sd, ymin=value-sd)) +
      facet_wrap(~variable, ncol=2, scales = "free_y") +
      theme_bw() +
      scale_x_continuous(breaks = min(posterior_probs_summary$K):max(posterior_probs_summary$K))
    suppressWarnings(print(gg))
  }

  best <- posterior_probs_summary[posterior_probs_summary$variable=="delta K",]
  best <- best$K[which.max(best$value)]
  message("Probs not the best K, though hey. lol.")
  return(best)

}
.bestK_structure <- function(posterior_probs, make_plot) {
  # look for change point in log likelihood plot
  #summarise data by K
  posterior_probs_summary <- group_by(posterior_probs, K) %>%
    summarize(
      runs = n(),
      meanK = mean(pos_prob),
      sdK = sd(pos_prob)
    )

  #generate plot
  if (make_plot) {
    gg <- ggplot(posterior_probs_summary, aes(x=K, y=meanK)) +
      geom_point() +
      geom_errorbar(aes(ymax=meanK+sdK, ymin=meanK-sdK), width=0.1) +
      theme_bw() +
      scale_x_continuous(breaks = min(posterior_probs_summary$K):max(posterior_probs_summary$K)) +
      ylab("Mean log posterior probability")
    suppressWarnings(print(gg))
  }

  message("Probs not the best K, though hey. lol.")
  return(posterior_probs_summary$K[which.max(posterior_probs_summary$meanK)])

}
.bestK_admixture <- function(log_df, make_plot) {
  message("Probs not the best K, though hey. lol.")

  if (make_plot) {
    gg <- ggplot(log_df, aes(x = K, y = CVerror)) + geom_point() + theme_bw()
    suppressWarnings(print(gg))
  }
  best <- log_df$K[which.min(log_df$CVerror)]
  return(best)
}

is.sequential <- function(x){
  x <- as.numeric(x)
  all(diff(x) == diff(x)[1])
}
