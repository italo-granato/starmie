# starmie-methods.R

# print method (could add in extra info)
print.starmie <- function(starmie_obj) {
  statement <- paste("starmie object with\n",
        starmie_obj$n_samples, "samples\n",
        starmie_obj$n_markers, "markers with ploidy equal to", starmie_obj$ploidy, "\n")
  if(!is.null(starmie_obj$admixture_run$q_info)) {
    statement <- paste(statement, "ADMIXTURE run with K =",
                       paste(unique(starmie_obj$admixture_run$q_info$K), collapse= " "))
  }
  cat(statement)
}

#' ggplot method for starmie object
#' @param starmie_obj an object of class \code{\link{starmie}}
#' @param method a character corresponding to either 'admixture' or 'structure'
#' @param type a character corresponding to type of plot to make (default 'all')
#' @import ggplot2
autoplot.starmie <- function(starmie_obj, method, type = "all", ...) {
  # i/o checks
  if (!method %in% c("admixture", "structure") && length(method) != 1) {
    stop("Invalid method command choose either structure or admixture")
  }

  # theme for starmie object
  starmie_theme <- theme_bw() + theme(panel.grid=element_blank(),
                                      panel.border=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      panel.grid.major=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank(),
                                      legend.position="none")

  # probably best to do case, switch here
  if (method == "admixture") {
    q_df <- get_admixture(starmie_obj)$q_info
    if (nrow(q_df) == 0) {
      stop("No ADMIXTURE clustering information in starmie object")
    }

    admixture_plot <- ggplot(q_df, aes(x = sample.id,
                                       y = probability,
                                       fill = factor(cluster))) +
      facet_grid(K ~ .) +
      geom_bar(stat = "identity") +
      scale_y_continuous(expand = c(0,0)) +
      xlab("Sample") +
      ylab("Proportion of cluster") +
      starmie_theme

    admixture_plot

  }

}

#' Return a tidy data frame from a starmie object for custom plotting
#'
#' @importFrom ggplot2 fortify
fortify.starmie <- function(starmie_obj, method) {

}

