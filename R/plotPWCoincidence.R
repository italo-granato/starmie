# plot the Pairwise coincidence, showing the MAP state

#' fineSTRUCTURE pairwise coincidence plot
#'
#' @param x a \code{\link{fineStruct}} object
#' @param cex (default=0.7) number indicating the amount by which plotting text and symbols should be scaled relative to the default.
#' @param ... Additional arguments to be passed to ape::plot.phylo
#' @export
#' @examples
#' # fineStruct example
#' chunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.chunkcounts.out", package = "starmie")
#' treefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.tree.xml", package = "starmie")
#' mcmcfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mcmc.xml", package = "starmie")
#' fineData <- loadFineStructure(chunkfile, treefile, mcmcfile)
#' plotPWCoincidence(fineData)
plotPWCoincidence <- function(x, ...) {
  UseMethod("plotPWCoincidence", x)
}

# #' @method plotFineStructDendro phylo
# #' @export
# plotFineStructDendro.phylo <- function(x, cex=0.7, ...) {
#   plot(x, direction = "right", label.offset = 0.1, cex = cex, ...)
#   nodelabels(x$node.label, frame = "none", adj = c(2,-0.5), cex=cex)
# }

#' @method plotPWCoincidence fineStruct
#' @export
plotPWCoincidence.fineStruct <- function(x, ...) {

  h <- as.hclust(x$dendro)
  plot(h)
  text(h, col = 2, cex = 0.7, font = NULL, float = 0.01,
       print.num = print.num)

  }




