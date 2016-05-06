# test-bestK.R
context("checks for bestK")

structure_files <- system.file("extdata/microsat_testfiles", package="starmie")
structure_output_files <- list.files(structure_files, pattern = ".*out_f", full.names = TRUE)
structure_log_files <- list.files(structure_files, pattern = ".*log", full.names = TRUE)
structure_runs <- mapply(loadStructure, structure_output_files, structure_log_files, SIMPLIFY=FALSE)

test_that("loadStructure I/O checks", {
  expect_error(bestK(c(1,2,3)), "structure_runs must be a vector of struct objects.")
  expect_error(bestK(structure_runs, method=NA), "method must be one of 'evano' or 'structure'")
  expect_error(bestK(structure_runs, plot="K"), "plot must be one of TRUE or FALSE")
})

