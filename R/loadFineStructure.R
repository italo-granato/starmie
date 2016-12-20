# loadFineStructure.R
# Functions for parsing output from fineSTRUCTURE

#' Read fineSTRUCTURE Output
#' @param chunkfile a string containing a chunk counts file produce by fineStructure
#' @param treefile a string containing a tree file produce by fineStructure
#' @param mcmcfile a string containing a mcmc file produce by fineStructure
#' @importFrom data.table fread
#' @import xml2
#' @importFrom stringr str_split
#' @importFrom ape read.tree
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

  #Load tree file
  xml <- read_xml(treefile)
  fine_obj$dendro <- read.tree(
    text=xml_text(xml_find_all(xml,".//Tree"))
  )

  #Load mcmc file
  xml <- read_xml(mcmcfile)
  xml <- xml_find_all(xml,".//Iteration")
  xml_attrs(xml[[1]])
  iterations <- xml$doc$children$outputFile[which(names(xml$doc$children$outputFile)=="Iteration")]
  names(iterations[[1]])
}


as.data.frame.myres<-function(txml){
  ## Converts our xml file format into a matrix, one row per iteration
  tmpits<-xml$doc$children$outputFile[which(names(xml$doc$children$outputFile)=="Iteration")]
  #	xmlChildren(tmpits[[1]])
  cnames<-names(tmpits[[1]])
  res<-as.data.frame(matrix(nrow=length(tmpits),ncol=length(cnames)))
  for(i in 1:length(cnames)) {
    res[,i]<-extractValue(xml,cnames[i])
  }
  colnames(res)<-cnames
  res[,which(!names(res) %in% c("Pop","P","Q"))]  <-apply(res[,which(!names(res) %in% c("Pop","P","Q"))],2,as.numeric)
  res
}

extractValue<-function(txml,v,getNames=FALSE){
  ## Important utility function for extracting element v from the xml
  if(getNames) num<-sapply(txml,getV,v="Number")
  res<-sapply(txml,getV,v=v)
  if(!getNames) names(res)<-NULL
  else names(res)<-num
  return(res)
}
getV<-function(it,v){
  ## extract element v from iteration it
  return(xmlValue(xmlChildren(it)[[v]]))
}




