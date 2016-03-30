# test-loadAdmixture.R
context("I/O checks for loadAdmixture")

new_starmie <- starmie()
sample_file <- system.file("extdata/hapmap3_files", "hapmap3.fam", package="starmie")
hapmap_files <- system.file("extdata/hapmap3_files", package="starmie")
qfiles <- list.files(hapmap_files, pattern = ".Q$", full.names = TRUE)
logfiles <- list.files(hapmap_files, pattern = ".out$", full.names = TRUE)

test_that("loadAdmixture I/O checks", {
  expect_error(loadAdmixture(list(), qfiles), "Not a valid starmie object.")
  expect_error(loadAdmixture(new_starmie, c(NA, qfiles)), "file_list must be a character vector.")
  expect_error(loadAdmixture(new_starmie, c(1,2,3)), "file_list must be a character vector.")
  expect_error(loadAdmixture(new_starmie, qfiles, c(NA, logfiles)), "logfile_list must be a character vector.")
  expect_error(loadAdmixture(new_starmie, qfiles, logfiles), "Sample metadata is required.")
})


