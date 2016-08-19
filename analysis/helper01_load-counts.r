## function to load and validate the "raw" count data

load_counts <- function() {
  require(testthat)  
  filteredSailfishCounts <- # take a few moments
    read.delim("../data/consolidated-lowExpFiltered-Sailfish-results.15July.txt",
               colClasses = c("character", rep("numeric", 24)), header = TRUE, row.names = 1)
  test_that("filtered Sailfish data has 38197 rows upon import",
            expect_equal(38197, nrow(filteredSailfishCounts)))
  test_that("Sailfish data has data for exactly 24 samples",
            expect_equal(24, ncol(filteredSailfishCounts)))
  return(filteredSailfishCounts)
}