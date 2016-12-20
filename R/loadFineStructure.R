# loadFineStructure.R
# Functions for parsing output from fineSTRUCTURE

#' Read fineSTRUCTURE Output
#' @param chunkfile a string containing a chunk counts file produce by fineStructure
#' @param treefile a string containing a tree file produce by fineStructure
#' @param mcmcfile a string containing a mcmc file produce by fineStructure
#' @importFrom data.table fread
#' @importFrom stringr str_split
#' @export
#' @examples
#' # read in example fineSTRUCTURE output
#' chunkfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.chunkcounts.out", package = "starmie")
#' treefile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.tree.xml", package = "starmie")
#' mcmcfile <- system.file("extdata/fine_structure_files", "EastAsiaSimple.EMlinked.mcmc.xml", package = "starmie")
#' fineData <- loadFineStructure(chunkfile, treefile, mcmcfile)
#' fineData
loadFineStructure <- function(chunkfile, treefile, mcmcfile){
  # i/o checks
  if (!is.na(chunkfile) & !is.character(chunkfile)) stop("chunkfile must be a string.")
  if (!is.na(treefile) & !is.character(treefile)) stop("treefile must be a string.")
  if (!is.na(mcmcfile) & !is.character(mcmcfile)) stop("mcmcfile must be a string.")
  if(!requireNamespace("XML", quietly = TRUE)) stop("XML package not installed, please install it")
  if(!requireNamespace("ape", quietly = TRUE)) stop("ape package not installed, please install it")

  #create new structure object
  fine_obj <- fineStruct()

  # do the work (much of this code is adapted from Daniel Lawson's R scripts https://people.maths.bris.ac.uk/~madjl/finestructure/index.html)
  #Load the chunk counts file
  fine_obj$cfactor <- as.numeric(str_split(readLines(chunkfile, n=1), " ", simplify = TRUE)[[2]])
  fine_obj$chunkcounts_df <- fread(chunkfile, skip=1, header=TRUE, data.table=FALSE)
  fine_obj$nsamples <- ncol(fine_obj$chunkcounts_df)-1

  #Load tree file
  xml <- xml2::read_xml(treefile)
  fine_obj$dendro <- ape::read.tree(
    text=xml2::xml_text(xml2::xml_find_all(xml,".//Tree"))
  )

  #Load mcmc file
  xml <- xml2::read_xml(mcmcfile)
  xml <- xml2::xml_find_all(xml,".//Iteration")
  fine_obj$mcmc_df <- data.frame(do.call(rbind
                        , lapply(xml, function(x) xml2::xml_text(xml2::xml_children(x)))
                        )
                        , stringsAsFactors = FALSE)
  fine_obj$mcmc_df[] <- lapply(fine_obj$mcmc_df
                               , function(x) type.convert(as.character(x), as.is=TRUE))
  colnames(fine_obj$mcmc_df) <- xml2::xml_name(xml2::xml_children(xml[[1]]))

  return(fine_obj)
}
