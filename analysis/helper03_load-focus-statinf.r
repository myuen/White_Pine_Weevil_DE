## function to load and validate the statistical inference results for our
## "focus" terms

load_focus_statInf <- function() {
  require(testthat)
  
  sidf <- read.delim("../results/limma-results-focus-terms.tsv")
  
  test_that("inference results for our focus terms (still) have 755,352 rows",
            expect_equal(58104 * 7, nrow(sidf)))
  
  sidf$focus_term <- 
    factor(sidf$focus_term, levels = 
             c("constDiff", "woundResp_Q903", "woundResp_H898", "weevilInd_Q903",
               "weevilInd_H898", "weevilCtrl_Q903", "weevilCtrl_H898"))


  options(digits = 8)
  (t_medians <- aggregate(t ~ focus_term, sidf, median))
  
  test_that("medians of the t statistics are what we expect",
            expect_equal(t_medians$t,
                         c(0.11699985, 0.0497448213, -0.0021461646, -0.2089632477,
                           -0.0742546845, -0.0917739088, -0.0934183638),
                         tolerance = 2 * .Machine$double.eps ^ 0.5))
  
  return(sidf)
}
