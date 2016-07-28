# plot-methods for visualising full Structure model
#' Plot allele frequencies within inferred cluster
#'
#' @param structure_obj a \code{\link{struct}} object
#' @param cluster_num an integer value for the cluster of interest
#' @param relevance relevance metric applied to ranking
#' @details This function plots the estimated allele frequencies within a genetic
#' cluster and contrasts them to the overall population allele frequencies.
#' @import ggplot2
#' @importFrom data.table setorder melt
#' @export
#' @examples
#' k6_data <- exampleStructure("barplot")
#' plotRanks(k6_data, cluster_num = 2)
plotRanks <- function(structure_obj, cluster_num, relevance = NULL) {
  # i/o checks
  if (!inherits(structure_obj, "struct"))
    stop("Not a structure object")
  if (cluster_num < 1 | cluster_num > getK(structure_obj))
    stop(paste("cluster_num must be between 1 and", getK(structure_obj)))

  # get overall term frequencies
  within_cluster_freqs <- getClusterAlleleFreqMat(structure_obj)[, c(1,2, cluster_num+2)]
  within_cluster_freqs <- data.frame(locus = within_cluster_freqs[,1],
                                     allele_num = within_cluster_freqs[,2],
                                     cluster_freq = within_cluster_freqs[,3])
  total_freqs  <- getCompleteAlleleFreqMat(structure_obj)
  total_freqs <- data.frame(locus = total_freqs[,1],
                            allele_num = total_freqs[,2],
                            population_freq = total_freqs[,3])

  combined_freqs <- merge(total_freqs, within_cluster_freqs)
  setorder(combined_freqs, -cluster_freq)
  cfreqs_tidy <- melt(combined_freqs, id.vars = c("locus", "allele_num"),
                      measure.vars = c("population_freq", "cluster_freq"),
                      variable.name = "category",
                      value.name = "freq")
  ggplot(cfreqs_tidy,
         aes(x = factor(locus),
             y = freq, colour = category)) +
    geom_pointrange(aes(ymin = 0, ymax = freq)) +
    geom_point() +
    coord_flip() +
    theme_bw()
}

#' Plot prinicipal coordinates between inferred clusters
#' @param structure_obj a \code{\link{struct}} object
#' @param method one of "nnd" or "jsd"
#' @details "nnd" uses the nucleotide distance matrix estimated by STRUCTURE
#' to construct the principal coordinates, sizing the points by the expected
#' heterozygosity within a cluster. "jsd" produces a principal coordinates
#' from the Jensen Shannon Divergence metric as used by the 'ldavis' package.
#' @import ggplot2
#' @importFrom  proxy dist
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples
#' k6_data <- exampleStructure("barplot")
#' plotMDS(k6_data)
#' plotMDS(k6_data, method = "jsd")
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
    ggplot(mds_df, aes(x = PC1, y = PC2, size = Relative.Contribution)) +
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

