# plot-methods for visualising full Structure model
#

#' Plot prinicipal compenents between inferred clusters
#' @param structure_obj
#' @import ggplot2
#'
plotMDS <- function(structure_obj, method = "nnd") {
  # gather visual elements

  if (method == "nnd") {
    # use structure allele-freqs diveragnes as distance matrix
    dist_xy <- structure_obj$allele_freqs
    # set diags to 0
    diag(dist_xy) <- 0
    # grab cluster expected heterzygosities
    clust_eh <- structure_obj$avg_dist_df

    # compute mds coordinates
    mds_clust <- cmdscale(dist_xy)

    mds_df <- data.frame(clust_eh, PC1 = mds_clust[,1], PC2 = mds_clust[,2])

    ggplot(mds_df, aes(x = PC1, y = PC2, size = Avg.dist)) + geom_point() + theme_bw()


  } else {
    return(NULL)
  }

}
