## function to load and validate the statistical inference results for our
## "focus" terms

load_focus_statInf <- function() {
  require(testthat)
  
  sidf <- read.delim("../results/limma-results-focus-terms.tsv")
  
  test_that("inference results for our focus terms (still) have 755,352 rows",
            expect_equal(58104 * 13, nrow(sidf)))
  
  sidf$focus_term <- 
    factor(sidf$focus_term, levels = 
             c("wound_Q903", "wound_H898", "feed_Q903", "feed_H898",
               "combined_effect_Q903", "combined_effect_H898", 
               "ctrl_vs_ctrl", "wound_vs_wound","gallery_vs_gallery", 
               "wounding_diff", "feeding_diff", "combined_diff",
               "induced_vs_const"))

  options(digits = 10)
  (t_medians <- aggregate(t ~ focus_term, sidf, median))
  
  test_that("medians of the t statistics are what we expect",
            expect_equal(t_medians$t,
                         c(0.049744821, -0.002146165, -0.208963248, -0.074254685, 
                           -0.091773909, -0.093418364, 0.116999851, 0.015444845,
                           0.198138740, -0.051140065, 0.103334093, 0.059183495,
                           -0.342118463), 
                         tolerance = 2 * .Machine$double.eps ^ 0.5))
  
  return(sidf)
}
