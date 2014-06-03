## function to load and validate the experimental design

load_expDes <- function() {
  require(testthat)  
  expDes <- read.delim("../data/White_Pine_Weevil_exp_design.tsv",
                       stringsAsFactors = FALSE)
  expDes <-
    mutate(expDes,
           gType = factor(gType, levels = c("Q903susc", "H898res")),
           txCode = factor(txCode, levels = c('C', 'W', 'G')),
           tx = factor(tx, levels = c("Control", "Wound", "Gallery")))
  test_that("design matrix has 24 rows upon import",
            expect_equal(24, nrow(expDes)))
  return(expDes)
}