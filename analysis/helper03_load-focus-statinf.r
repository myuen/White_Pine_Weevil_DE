## function to load and validate the statistical inference results for our
## "focus" terms

load_focus_statInf <- function() {
  require(testthat)
  
  sidf <- read.delim("../results/limma-results-focus-terms.15July.tsv", 
                     colClasses = c(rep("character", 3), rep("numeric", 6)))
  
  test_that("inference results for our focus terms have 267379 rows",
            expect_equal(38197 * 7, nrow(sidf)))
  
  sidf$focus_term <- 
    factor(sidf$focus_term, levels = 
             c("constDiff", "woundResp_Q903", "woundResp_H898", "weevilInd_Q903",
               "weevilInd_H898", "weevilCtrl_Q903", "weevilCtrl_H898"))

  options(digits = 15)
  (t_medians <- aggregate(t ~ focus_term, sidf, median))
  
  test_that("medians of the t statistics are what we expect",
            expect_equal(t_medians$t,
                         c(0.11886468313334492, 0.04949321203384519, -0.00856661400850118,
                           -0.09816757866141612, -0.10336498104535216,
                           -0.00798612294370277, -0.10327785045908873),
                         tolerance = 2 * .Machine$double.eps ^ 0.5))
  
  return(sidf)
}
