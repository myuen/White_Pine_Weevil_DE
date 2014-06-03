#' While doing the main DEA, I ran into problems with the column names of the 
#' design matrix. Initially, I let them arise naturally, from my factor names, 
#' levels, and model form. But then you get a fatal error from 
#' `makeContrasts()`, because at some point they apparently need to be valid R 
#' object names. In some of my informal fussing around, I altered the model 
#' matrix column names enough to satisfy `makeContrasts()`, but then had trouble
#' getting what I wanted from `topTable()`. Specifically, I had trouble 
#' specifying the `coef` or term I wanted to see inference for. Here is where I 
#' document the initial problem and the solution I've settled on.
#' 
#' `spin()` this "by hand" with `precious = TRUE` to ensure markdown is
#' preserved.

library(limma)
library(edgeR)
library(testthat) # facilitate tests that will catch changes on re-analysis

#' ## Differential Expression Analysis on Sitka Spruce Weevil 
#' ## Experiment with limma + voom

#' Load counts from Sailfish. Each row represents a contig.
filteredSailfishCounts <- # take a few moments
  read.delim("../data/consolidated-filtered-Sailfish-results.txt")
str(filteredSailfishCounts) # 'data.frame':  65609 obs. of  24 variables:
test_that("filtered Sailfish data (still) has 65609 rows upon import",
          expect_equal(65609, nrow(filteredSailfishCounts)))
test_that("Sailfish data (still) has data for exactly 24 samples",
          expect_equal(24, ncol(filteredSailfishCounts)))

#' Load experimental design
expDes <- read.delim("../data/White_Pine_Weevil_exp_design.tsv",
                     stringsAsFactors = FALSE)
expDes <-
  transform(expDes,
            gType = factor(gType, levels = c("Q903susc", "H898res")),
            txCode = factor(txCode, levels = c('C', 'W', 'G')),
            tx = factor(tx, levels = c("Control", "Wound", "Gallery")))
expDes$grp <-
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:
test_that("design matrix (still) has 24 rows upon import",
          expect_equal(24, nrow(expDes)))

#' Load counts into DGEList object from edgeR package.
y <- DGEList(counts = filteredSailfishCounts, group = expDes$grp)

#' TMM Normalization by Depth
y <- calcNormFactors(y)

#' make model matrix
#modMat <- model.matrix(~ tx/gType - 1, expDes)
modMat <- model.matrix(~ gType * tx, expDes)

colnames(modMat)

#' voom transformation
#+ first-voom
v <- voom(y, modMat, plot = TRUE) # take a couple moments

#' Linear modelling
fit <- lmFit(v, modMat)

#' create contrast matrix to take difference of interaction terms
#+ error=TRUE
cont_matrix <-
  makeContrasts(weevil = gTypeH898res:txGallery - gTypeH898res:txWound,
                levels = modMat) # here's our first fatal error

#' The `:` is not allowed in R object names. This seems perverse, given that 
#' they will naturally occur in interaction terms formed by model.matrix(). See 
#' [here](https://stat.ethz.ch/pipermail/bioconductor/2010-February/031554.html)
#' for confirmation that I am not just being stupid.
#' 
#' Let's change the column names of the model matrix and redo all the linear 
#' modelling bits. I'm also going to get rid of the parentheses flanking
#' `Intercept`.

colnames(modMat)
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat)

#' voom transformation
#+ second-voom
v <- voom(y, modMat, plot = TRUE) # take a couple moments

#' Linear modelling
fit <- lmFit(v, modMat)

#' create contrast matrix to take difference of interaction terms
cont_matrix <-
  makeContrasts(weevil = gTypeH898res_txGallery - gTypeH898res_txWound,
                levels = modMat)

#' Put the contrasts into action and apply the Empirical Bayes moderation.
fit2 <- contrasts.fit(fit, cont_matrix)
fit3 <- eBayes(fit2)

#' Pull out top hits for our term of interest
tt_diff_interactions <- topTable(fit3, coef = "weevil")
head(tt_diff_interactions)

sessionInfo()