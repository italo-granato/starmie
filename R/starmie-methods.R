# starmie-methods.R

# print method (could add in extra info)
print.starmie <- function(starmie_obj) {
  statement <- paste("starmie object with\n",
        starmie_obj$n_samples, "samples\n",
        starmie_obj$n_markers, "markers with ploidy equal", starmie_obj$ploidy, "\n")
  cat(statement)
}

# ggplot method not sure if have to use ggproto here, look into this
