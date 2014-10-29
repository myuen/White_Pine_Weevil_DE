## function to load and validate the statistical inference results for our
## "focus" terms

load_focus_statInf <- function() {
  require(testthat)
  
  sidf <- read.delim("../results/limma-results-focus-terms.tsv")
  
  test_that("inference results for our focus terms (still) have 825,917 rows",
            expect_equal(65609 * 13, nrow(sidf)))
  
  sidf$focus_term <- 
    factor(sidf$focus_term, levels = 
             c("wounding_in_Q903", "wounding_in_H898", "feeding_in_Q903", "feeding_in_H898",
               "comb_effect_in_Q903", "comb_effect_in_H898", "ctrl_vs_ctrl", "wound_vs_wound",
               "gallery_vs_gallery", "wounding_diff", "feeding_diff", "combined_diff",
               "induced_vs_const"))

  options(digits = 10)
  (t_medians <- aggregate(t ~ focus_term, sidf, median))
  
  test_that("medians of the t statistics are what we expect",
            expect_equal(t_medians$t,
                         c(0.092882336, 0.012555037, -0.155382913, -0.090723545,
                           0.026698663, -0.086665900,  0.157361758, 0.011987382,
                           0.107238311, -0.067348886,  0.048179395,  0.005007678,
                           -0.312069443), tolerance = 1.00293e-07))
#   ,tolerance = 2 * .Machine$double.eps ^ 0.5))
  
  return(sidf)
}
