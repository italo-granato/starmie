# test-loadStructure.R
context("checks for loadStructure")

new_starmie <- starmie()
sample_file <- system.file("extdata/hapmap3_files", "hapmap3.fam", package="starmie")
structure_files <- system.file("extdata/microsat_testfiles", package="starmie")
structure_output_files <- list.files(structure_files, pattern = ".*out_f", full.names = TRUE)
structure_log_files <- list.files(structure_files, pattern = ".*log", full.names = TRUE)

test_that("loadStructure I/O checks", {
  expect_error(loadStructure(list(), structure_output_files), "Not a valid starmie object.")
  expect_error(loadStructure(new_starmie, c(NA, structure_output_files)), "file_list must be a character vector.")
  expect_error(loadStructure(new_starmie, c(1,2,3)), "file_list must be a character vector.")
  expect_error(loadStructure(new_starmie, structure_output_files, c(NA, structure_log_files)), "logfile_list must be a character vector.")
  expect_error(loadStructure(new_starmie, structure_output_files, structure_log_files), "Sample metadata is required.")
})

test_that("loadStructure population mismatch check", {
  expect_error(loadStructure(my_starmie, structure_output_files[c(2,3,1)], structure_log_files)
               , "Population mismatch between output and logfile.")
})
