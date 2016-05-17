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
#' @import ggplot2
#' @import dplyr
plotMCMC <- function(struct_list, make_plot = TRUE, use_logL = TRUE, facet = TRUE) {
  #i/o checks

  # generate data frame of mcmc diagnostics from non-burn iterations
  mcmc_df <- do.call("rbind", lapply(struct_list, getMCMC))
  # run number is kind of arbirtary here
  mcmc_df2 <- mcmc_df %>%
    group_by(K, Iteration) %>%
    mutate(run = row_number())

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
    return(mcmc_df2)
  }

  return(mcmc_df2)

}
