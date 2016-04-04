# loadStructure.R
# Functions for parsing output from STRUCTURE

#' Read Structure Output
#' @param starmie object
#' @param  location of the STRUCTURE output file
#' @param logfile include STRUCTURE logs for each run (required if wanting to do MCMC plots)
#' @importFrom readr read_lines
#' @import stringr
#' @importFrom purrr map
#' @export
readStructure <- function(starmie_obj, filename, logfile = NULL) {
  # i/o checks
  if (!inherits(starmie_obj, "starmie")) stop("Not a valid starmie object.")
  if(!is.character(filename)) {
    stop("File must be character variable.")
  }
  if(!file.exists(filename)) {
    stop("File does not exist.")
  }

  # do the work
  s_f <- readr::read_lines(filename)

  #remove blank lines
  s_f <- s_f[s_f!=""]

  #Get run parameters
  run_lines <- s_f[(which(s_f=="Run parameters:")+1):(which(grepl("^RANDOMIZE.*", s_f))-1)]
  run_lines <- str_trim(run_lines)
  run_lines <- str_split_fixed(run_lines, " ", n=2)
  run_params <- data.frame(Parameter=run_lines[,2], Value=as.numeric(run_lines[,1]))
  pops <- run_params[run_params$Parameter=="populations assumed",][2]

  #Get membership proportion
  mem_lines <- s_f[(which(grepl("^Proportion of membership.*",s_f))+3):(which(grepl("^Allele-freq", s_f))-2)]
  mem_lines <- str_trim(mem_lines)
  mem_lines <- str_split_fixed(mem_lines, "\\s+", n=pops+2)
  mem_lines[,1] <- str_replace(mem_lines[,1], ":", "")
  mem_df <- data.matrix(data.frame(mem_lines[-1, ], stringsAsFactors = FALSE))
  colnames(mem_df) <- mem_lines[1,]

  #Get Allele-freq
  alle_lines <- s_f[(which(grepl("^Allele-freq", s_f))+4):which(grepl("^Average distances.*", s_f))-1]
  alle_lines <- str_trim(alle_lines)
  alle_lines <- str_split_fixed(alle_lines, "\\s+", n=pops+1)
  alle_freqs <- alle_lines[,2:ncol(alle_lines)]
  class(alle_freqs) <- "numeric"

  #Average distances
  avg_dist_lines <- s_f[(which(grepl("^Average distances.*", s_f))+1):(which(grepl("^Estimated Ln.*", s_f))-2)]
  avg_dist_lines <- str_trim(avg_dist_lines)
  avg_dist_lines <- str_split_fixed(avg_dist_lines, "  : ", n=2)
  avg_dist_lines[,1] <- str_replace(avg_dist_lines[,1], "cluster  ", "")
  class(avg_dist_lines) <- "numeric"
  avg_dist_df <- data.frame("Cluster"=avg_dist_lines[,1], "Avg.dist"=avg_dist_lines[,2])

  #Model fit stats
  fit_lines <- s_f[(which(grepl("^Estimated Ln.*", s_f))):(which(grepl("^Mean value of Fst_1.*", s_f))-1)]
  fit_lines <- str_trim(fit_lines)
  fit_lines <- str_split_fixed(fit_lines, "= ", n=2)
  fit_lines[,1] <- str_trim(fit_lines[,1])
  fit_stats_df <- data.frame(Statistic=fit_lines[,1], Value=as.numeric(fit_lines[,2]))

  #Fst values
  fst_lines <- s_f[which(grepl("^Mean value of Fst_1.*", s_f)):(which(grepl("^Inferred ancestry of.*", s_f))-1)]
  fst_lines <- str_trim(fst_lines)
  fst_lines <- str_split_fixed(fst_lines, "= ", n=2)
  fst_lines[,1] <- str_trim(fst_lines[,1])
  fst_lines[,1] <- str_replace(fst_lines[,1], "Mean value of Fst_", "")
  fst_df <- data.frame(Fst.Group=as.numeric(fst_lines[,1]), Value=as.numeric(fst_lines[,2]))

  #Inferred ancestory individuals
  ances_lines <- s_f[(which(grepl("^Inferred ancestry of.*", s_f))+1):(which(grepl("^Estimated Allele Frequencies .*", s_f))-1)]
  ances_lines <- str_trim(ances_lines)
  ances_lines <- str_split_fixed(ances_lines, "[\\s:]+", n=6)
  header <- ances_lines[1,][1:3]
  ances_lines <- ances_lines[2:nrow(ances_lines),2:ncol(ances_lines)]
  ances_lines[,2] <- str_replace_all(ances_lines[,2], "[()]", "")
  class(ances_lines) <- "numeric"
  ancest_df <- data.frame(ances_lines)
  colnames(ancest_df)[1:3] <- header

  #Cluster allele frequencies
  clust_allel_lines <- s_f[(which(grepl("^First column gives.*", s_f))+1):(which(grepl("^Values of parameters used.*", s_f))-1)]
  pos <- which(grepl("^Locus .*", clust_allel_lines))
  clust_allel_lines <- unname(split(clust_allel_lines, cumsum(seq_along(clust_allel_lines) %in% pos)))
  clust_allele_list <- purrr::map(clust_allel_lines, function(x){
    list(Locus=str_split(x[[1]], "\\s+")[[1]][2]
         , AlleleNumber=str_split(x[[2]], "\\s+")[[1]][1]
         , MissingDataPercentage=str_split(x[[3]], "\\s+")[[1]][1]
         , FreqMatrix=str_split_fixed(str_trim(x[4:length(x)]), "\\s+", n=4)
    )
  })

  starmie_obj$structure_run <- list(pops=pops
                                    , run_params=run_params
                                    , mem_df=mem_df
                                    , alle_freqs=alle_freqs
                                    , avg_dist_df=avg_dist_df
                                    , fit_stats_df=fit_stats_df
                                    , fst_df=fst_df
                                    , ancest_df=ancest_df
                                    , clust_allele_list=clust_allele_list)
  if(!is.null(logfile)){
    # i/o checks
    if(!is.character(logfile)) {
      stop("File must be character variable.")
    }
    if(!file.exists(logfile)) {
      stop("File does not exist.")
    }

    #do more work
    l_f <- readr::read_lines(logfile)
    #remove blank lines
    l_f <- l_f[l_f!=""]

    #Get burn and non-burn in iterations as seperate data frames.
    ##NOTE: This relies heavily on the current format.
    intervals <- unlist(lapply(which(grepl("^ Rep#:   Lambda   Alpha.*",l_f)), function(x){ x[1]:(x[1]+11) }))
    burn_lines <- l_f[intervals]
    burn_lines <- burn_lines[grepl("[0-9]|R.*", burn_lines)]
    burn_lines <- str_trim(burn_lines)
    mid <- which(grepl(".*Ln Like  Est Ln P.*",burn_lines))[1]
    nonburn_lines <- burn_lines[mid:length(burn_lines)]
    burn_lines <- burn_lines[1:mid-1]

    burn_header <- unlist(str_split(burn_lines[1], "\\s{2,}"))
    burn_lines <- burn_lines[!grepl("Rep#:   Lambda   Alpha.*", burn_lines)]
    burn_lines <- str_split_fixed(burn_lines, "\\s+", n=pops*2+4)
    burn_lines[,1] <- gsub(":","",burn_lines[,1])
    burn_df <- data.matrix(data.frame(burn_lines, stringsAsFactors = FALSE))
    colnames(burn_df) <- burn_header

    nonburn_header <- unlist(str_split(nonburn_lines[1], "\\s{2,}"))
    nonburn_lines <- nonburn_lines[!grepl("Rep#:.*", nonburn_lines)]
    nonburn_lines <- str_split_fixed(nonburn_lines, "\\s+", n=pops*2+5)
    nonburn_lines[,1] <- gsub(":","",nonburn_lines[,1])
    nonburn_df <- data.matrix(data.frame(nonburn_lines, stringsAsFactors = FALSE))
    colnames(nonburn_df) <- nonburn_header

    starmie_obj$structure_run["burn_df"] = burn_df
    starmie_obj$structure_run["nonburn_df"] = nonburn_df
  }

  starmie_obj
}







