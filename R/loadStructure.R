# loadStructure.R
# Functions for parsing output from STRUCTURE

#' Read Structure Output
#' @param  path directory containing output of STRUCTURE runs
#' @param  n_runs number of runs for each K value
#' @param logfile include STRUCTURE logs for each run (required if wanting to do MCMC plots)
#' @export
readStructure <- function(path, n_runs, logfile = FALSE) {
  # i/o checks
  if(!is.character(path)) {
    stop("Path must be character variable.")
  }
  if(!dir.exists(path)) {
    stop("Path does not exist.")
  }

  if(!is.finite(n_runs) && !is.integer(n_runs)) {
    stop("n_runs must be a finite integer")
  }

  # do the work

}


library(readr)
library(stringr)
library(purrr)


s_f <- read_lines("locprior.out_f")

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
header <- ances_lines[1,][1:5]
ances_lines <- ances_lines[2:nrow(ances_lines),2:ncol(ances_lines)]
ances_lines[,2] <- str_replace_all(ances_lines[,2], "[()]", "")
class(ances_lines) <- "numeric"
ancest_df <- data.frame(ances_lines)
colnames(ancest_df) <- header

#Cluster allele frequencies
clust_allel_lines <- s_f[(which(grepl("^First column gives.*", s_f))+1):(which(grepl("^Values of parameters used.*", s_f))-1)]
pos <- which(grepl("^Locus .*", clust_allel_lines))
clust_allel_lines <- unname(split(clust_allel_lines, cumsum(seq_along(clust_allel_lines) %in% pos)))
clust_allele_list <- purrr::map(clust_allel_lines, function(x){
  list(Locus=str_split(x[[1]][1], "\\s+")[[1]][2]
       , AlleleNumber=str_split(x[[1]][2], "\\s+")[[1]][1]
       , MissingDataPercentage=str_split(x[1][3], "\\s+")[[1]][1]
       , FreqMatrix=str_split_fixed(str_trim(x[1][4:length(x)]), "\\s+", n=4)
  )
})





