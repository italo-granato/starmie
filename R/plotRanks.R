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
