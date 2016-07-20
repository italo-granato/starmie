# test-bestK.R
context("checks for bestK")
# load sample dataset
multi_K <- exampleStructure("multiple_runs")
# single runs
single_runs <- structList(multi_K[1:10])
# non-sequential K
Ks <- unlist(lapply(multi_K, getK))
uneven_Ks <- structList(multi_K[Ks %in% c('9', '4', '2')])

test_that("loadStructure I/O checks", {
  # can't pass single struct object to test_that
  expect_error(bestK(multi_K[[1]]), "runs must be structList or admixList object")
  expect_error(bestK(multi_K, method=NA), "method must be one of 'evanno' or 'structure' or 'cverror'")
  expect_error(bestK(multi_K, make_plot="yes"), "plot must be one of TRUE or FALSE")
  # evanno method doesn't work if only single runs
  expect_error(bestK(single_runs, make_plot = FALSE))
  # structure method does output for sd should be all NA in results
  expect_true(all(is.na(bestK(single_runs, method = "structure", make_plot = FALSE)$sd)))
  expect_warning(bestK(uneven_Ks, make_plot = FALSE))
})

