## function to load and validate the statistical inference results for our
## "focus" terms

load_focus_statInf <- function() {
  require(testthat)
  
  sidf <- read.delim("../results/limma-results-focus-terms.tsv")
  
  test_that("inference results for our focus terms (still) have 852,800 rows",
            expect_equal(65600 * 13, nrow(sidf)))
  
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
                         c(0.093017448, 0.012387609, -0.164613806, -0.093503376, 
                           0.021539496, -0.089969465, 0.162680935, 0.011875199, 
                           0.118430987, -0.067642621, 0.052106222, 0.007929199, -0.321207705), 
                         tolerance = 2 * .Machine$double.eps ^ 0.5))
  
  return(sidf)
}
