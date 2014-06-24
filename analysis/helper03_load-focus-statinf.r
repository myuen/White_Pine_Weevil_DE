## function to load and validate the statistical inference results for our
## "focus" terms

load_focus_statInf <- function() {
  require(testthat)
  sidf <- read.delim("../results/limma-results-focus-terms.tsv")
  test_that("inference results for our focus terms (still) have 459,263 rows",
            expect_equal(459263, nrow(sidf)))
  
  sidf$focus_term <- 
    factor(sidf$focus_term, levels = 
             c("gTypeH898res_in_control", "gTypeH898res_in_wound",
             "gTypeH898res_in_gallery", "gTypeH898res_all",
             "weevil_in_gTypeQ903susc", "weevil_in_gTypeH898res", "weevil_diff"))

  (t_medians <- aggregate(t ~ focus_term, sidf, median))
  test_that("medians of the t statistics are what we expect",
            expect_equal(t_medians$t,
                         c(0.157361758, 0.01198738, 0.10723831,
                           -0.15538291, -0.09072355, 0.04817940),
                         tolerance = 2 * .Machine$double.eps ^ 0.5))

  return(sidf)
}
