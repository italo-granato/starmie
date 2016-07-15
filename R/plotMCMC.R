# plotMCMC.R
# Function for plotting MCMC chains for diagnostics

#' Plot STRUCTURE MCMC chains
#'
#' @param struct_list a list \code{\link{struct}} objects
#' @param make_plot logical print resulting plot default TRUE
#' @param facet logical facet by K default TRUE
#' @param use_logL logical plot log-likelihood (TRUE) or admixture coeffecient
#' @description Plot non-burn MCMC iterations of STRUCTURE for checking convergence.
#' If make_plot is set to FALSE a data.frame is returned containing the log likelihood
#' and alpha values over different K and runs and not plot is printed to the device.
#' @return If make_plot is TRUE a ggplot is printed to the screen and the
#' plot object and the data to generate it are returned. Otherwise,
#' a data.frame containing MCMC info it returned.
#' @import ggplot2
#' @importFrom dplyr mutate group_by bind_rows
#' @export
#' @examples
#' #Read in Structure files
#' multiple_runs_k10 <- exampleStructure("mcmc_diagnostics")
#' print(multiple_runs_k10)
#' results <- plotMCMC(multiple_runs_k10, make_plot = TRUE)
plotMCMC <- function(struct_list, make_plot = TRUE, use_logL = TRUE, facet = TRUE) {
  #i/o checks
  struct_class <- lapply(struct_list, class)
  if ( any(struct_class != "struct") ) {
    stop(" struct_list contains non struct objects ")
  }

  stopifnot(is.logical(make_plot))
  stopifnot(is.logical(use_logL))
  stopifnot(is.logical(facet))

  # generate data frame of mcmc diagnostics from non-burn iterations
  mcmc_df <- bind_rows(lapply(struct_list, getMCMC))
  # run number is kind of arbirtary here
  mcmc_df2 <- mutate(group_by(mcmc_df, K, Iteration),
                     run = row_number())

  if (make_plot) {
    if (use_logL) {
      gg <- ggplot(mcmc_df2,
                   aes(x = Iteration, y = LogL, colour = factor(run)))
        geom_line()
    } else {
      gg <- ggplot(mcmc_df2,
                   aes(x = Iteration, y = Alpha, colour = factor(run)))
    }

    gg <-  gg +
      geom_line() +
      scale_colour_hue("Run No.") +
      theme(legend.position = "bottom") + theme_bw()

    if (facet) {
      gg <- gg + facet_grid(K ~ . , scales = "free_y")
    }

    print(gg)
    return(list(mcmc_info = mcmc_df2, mcmc_plot = gg))
  }

  mcmc_df2

}
