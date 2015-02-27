## function to load and validate the "raw" count data

load_counts <- function() {
  require(testthat)  
  filteredSailfishCounts <- # take a few moments
    read.delim("../data/consolidated-filtered-Sailfish-results.txt")
  test_that("filtered Sailfish data has 58104 rows upon import",
            expect_equal(58104, nrow(filteredSailfishCounts)))
  test_that("Sailfish data has data for exactly 24 samples",
            expect_equal(24, ncol(filteredSailfishCounts)))
  return(filteredSailfishCounts)
}