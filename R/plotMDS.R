# plot-methods for visualising full Structure model
#

#' Plot prinicipal compenents between inferred clusters
#' @param structure_obj a \code{\link{struct}} object
#' @param method one of "nnd" or "jsd"
#' @details "nnd" uses the nucleotide distance matrix estimated by STRUCTURE
#' to construct the principal coordinates, sizing the points by the expected
#' heterozygosity within a cluster. "jsd" produces a principal coordinates
#' from the Jensen Shannon Divergence metric as used by the 'ldavis' package.
#' @import ggplot2
#' @importFrom  proxy dist
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

    ggplot(mds_df, aes(x = PC1, y = PC2, size = Avg.dist)) +
      geom_point() +
      ggrepel::geom_text_repel(
        aes(label = Cluster),
        size = 4,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      guides(size = guide_legend(title = "Expected Heterozygosity"))+
      theme_bw()


  } else if (method == "jsd") {
    Q <- getQ(structure_obj)
    dist_xy <- proxy::dist(t(Q), method = .JSD)
    mds_clust <- cmdscale(dist_xy)
    mds_df <- data.frame(Cluster = as.integer(sub("Cluster ", "", colnames(Q))),
                         Relative.Contribution = colSums(Q) / sum(Q),
                         PC1 = mds_clust[,1],
                         PC2 = mds_clust[,2])
    ggplot(mds_df, aes(x = PC1, y = PC2, size = Nk)) +
    geom_point() +
      ggrepel::geom_text_repel(
        aes(label = Cluster),
        size = 4,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      theme_bw()




  } else {
    stop("Not a valid method, must be either 'nnd' or 'ldavis'")
  }

}

.JSD<- function(x,y) sqrt(0.5 * .KLD(x, (x+y)/2) + 0.5 * .KLD(y, (x+y)/2))
.KLD <- function(x,y) sum(x * log(x/y))

