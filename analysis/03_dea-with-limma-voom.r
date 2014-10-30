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
str(x, list.len = 8) # 'data.frame':  65600 obs. of  24 variables:

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
    # A - Within genotype comparisons
    # A1 - Wound vs. Control
    wound_Q903 = txWound,
    wound_H898 = txWound + gTypeH898res_txWound,
    # A2 - Gallery vs. Wound
    feed_Q903 = txGallery - txWound,
    feed_H898 = txGallery - txWound +
      gTypeH898res_txGallery - gTypeH898res_txWound,
    # A3 - Gallery vs. Control (combined effect of wounding and feeding)
    combined_effect_Q903 = txGallery,
    combined_effect_H898 = txGallery + gTypeH898res_txGallery,
    # B - Between genotype comparisons
    # B1 - control vs control
    ctrl_vs_ctrl = gTypeH898res,
    # B2 - wound vs wound
    wound_vs_wound = gTypeH898res + gTypeH898res_txWound,
    # B3 - gallery vs gallery
    gallery_vs_gallery = gTypeH898res + gTypeH898res_txGallery,
    # C - Comparison of comparison
    # C1 -  Mwounding vs. wounding
    wounding_diff = gTypeH898res_txWound,
    # C2 - feeding vs. feeding i.e. (H898G - H898W) - (Q903G - Q903W)
    feeding_diff = gTypeH898res_txGallery - gTypeH898res_txWound,
    # C3 - combined effect vs. combined effect i.e. (H898G - H898C) - (Q903G - Q903C)
    combined_diff = gTypeH898res_txGallery,
    # D - Special case 
    # D1 - Induced vs Constitutive (Q903G - Q903W) - (H898C - Q903C)
    induced_vs_const = txGallery - txWound - gTypeH898res,
    levels = modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

summary(decideTests(fit3, p.value = 0.01, lfc = 2))
#    wound_Q903 wound_H898 feed_Q903 feed_H898 combined_effect_Q903
# -1          1          0      4632         4                 3896
# 0       65599      65600     59054     65596                58117
# 1           0          0      1914         0                 3587
#    combined_effect_H898 ctrl_vs_ctrl wound_vs_wound gallery_vs_gallery
# -1                    2         3679           3769               4390
# 0                 65597        57989          57971              55130
# 1                     1         3932           3860               6080
#    wounding_diff feeding_diff combined_diff induced_vs_const
# -1             0           99           175             5804
# 0          65600        64872         65151            56611
# 1              0          629           274             3185

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


test_that("stat inf on the focus terms has 852,800 rows",
          expect_equal(65600 * 13, nrow(statInf_focus_terms)))

(t_medians <- aggregate(t ~ focus_term, statInf_focus_terms, median))

all.equal(t_medians$t, c(0.093017448, 0.012387609, -0.164613806, -0.093503376, 
            0.021539496, -0.089969465, 0.162680935, 0.011875199, 
            0.118430987, -0.067642621, 0.052106222, 0.007929199, -0.321207705))

write.table(statInf_focus_terms,
            "../results/limma-results-focus-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)
