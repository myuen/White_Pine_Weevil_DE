#+ setup, include = FALSE
library(knitr)
opts_chunk$set(fig.path = 'figure/03-dea-with-limma-voom-')

library(edgeR)
library(plyr)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

#' Source purpose-built functions to load and validate the data and the 
#' experimental design. Then call them.
setwd("analysis/")
source("helper01_load-counts.r")
source("helper02_load-exp-des.r")

#' Load counts from Sailfish
x <- load_counts() # takes a few moments
str(x, list.len = 8) # 'data.frame':  58104 obs. of  24 variables

#' Load experimental design
expDes <- load_expDes()
expDes$grp <- # not really sure that we need this?
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:


#' Load counts into DGEList object from edgeR package.
y <- DGEList(counts = x, group = expDes$grp)

#' TMM Normalization by Depth
y <- calcNormFactors(y)

#' make model matrix
modMat <- model.matrix(~ gType * tx, expDes)

# I never needed to fit the model with this parametrization. Instead I got the
# effects of interests from the main parametrization and specific contrasts
# specifed below.
#modMat <- model.matrix(~ tx/gType - 1, expDes)

#' It's hard to believe, but the default names for columns associated with 
#' interaction terms will create fatal errors in makeContrasts below; prevent 
#' that, and a warning about the intercept, by modifying these column names 
#' here; see 90_limma-model-term-name-fiasco.r and .md for more
colnames(modMat)
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat)

#' voom transformation
#+ voom-plot
v <- voom(y, modMat, plot = TRUE) # take a couple moments

#' Linear modelling and forming contrasts
fit <- lmFit(v, modMat)

cont_matrix <-
  makeContrasts(
    # A - Between genotype comparisons
    # A1 - constitutive difference (H898 control - Q903 control)
    constDiff = gTypeH898res,
    # B - Within genotype comparisons
    # B1 - Weevil Induced (Gallery vs. Wound)
    weevilInd_Q903 = txGallery - txWound,
    weevilInd_H898 = txGallery - txWound +
      gTypeH898res_txGallery - gTypeH898res_txWound,
    # B2 - Gallery vs. Control (Weevil Control)
    weevilCtrl_Q903 = txGallery,
    weevilCtrl_H898 = txGallery + gTypeH898res_txGallery,
    levels = modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

summary(decideTests(fit3, p.value = 0.01, lfc = 2))
#    constDiff weevilInd_Q903 weevilInd_H898 weevilCtrl_Q903 weevilCtrl_H898
# -1      3377           4321              2            3618               2
# 0      51092          52302          58102           51825           58102
# 1       3635           1481              0            2661               0

focus_terms <- colnames(cont_matrix)

names(focus_terms) <- focus_terms

#' to get individual t statistics, etc., must use topTable() on each coef
#' separately
statInf_focus_terms <-
  adply(focus_terms, 1,
        function(x) name_rows(topTable(fit3, coef = x,
                                       number = Inf, sort.by = "none")))

statInf_focus_terms <-
  rename(statInf_focus_terms, c(X1 = "focus_term", .rownames = "contig"))

statInf_focus_terms <-
  statInf_focus_terms[ ,c("contig",
                          setdiff(names(statInf_focus_terms),"contig"))]

str(statInf_focus_terms)

test_that("stat inf on the focus terms has 755,352 rows",
          expect_equal(58104 * 5, nrow(statInf_focus_terms)))

(t_medians <- aggregate(t ~ focus_term, statInf_focus_terms, median))

all.equal(t_medians$t, c(0.11699985, -0.20896325, -0.07425468, -0.09177391, -0.09341836))

write.table(statInf_focus_terms,
            "../results/limma-results-focus-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)
