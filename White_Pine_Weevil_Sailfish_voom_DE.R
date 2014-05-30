library(edgeR)
library(ggplot2)
library(plyr)
library(reshape2)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# Load counts from Sailfish
filteredSailfishCounts <- # take a few moments
  read.delim("consolidated-filtered-Sailfish-results.txt")
str(filteredSailfishCounts, list.len = 8)
test_that("filtered Sailfish data has 65609 rows upon import",
          expect_equal(65609, nrow(filteredSailfishCounts)))
test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(filteredSailfishCounts)))

# Load experimental design
expDes <- read.delim("White_Pine_Weevil_exp_design.tsv",
                     stringsAsFactors = FALSE)
expDes <-
  mutate(expDes,
         gType = factor(gType, levels = c("Q903susc", "H898res")),
         txCode = factor(txCode, levels = c('C', 'W', 'G')),
         tx = factor(tx, levels = c("Control", "Wound", "Gallery")))
expDes$grp <-
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:
test_that("design matrix has 24 rows upon import",
          expect_equal(24, nrow(expDes)))

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = filteredSailfishCounts, group = expDes$grp)

# TMM Normalization by Depth
y <- calcNormFactors(y)

# make model matrix
modMat <- model.matrix(~ gType * tx, expDes)

# I never needed to fit the model with this parametrization. Instead I got the
# effects of interests from the main parametrization and specific contrasts
# specifed below.
#modMat <- model.matrix(~ tx/gType - 1, expDes)

# hard to believe, but the default names for columns associated with interaction
# terms will create fatal errors in makeContrasts below; prevent that, and a
# warning about the intercept, by modifying these column names here; see
# White_Pine_Weevil_Sailfish_limm-model-term-name-fiasco.R and .md for more
colnames(modMat)
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat)

# voom transformation
v <- voom(y, modMat, plot = TRUE) # take a couple moments

# Linear modelling
fit <- lmFit(v, modMat)
cont_matrix <-
  makeContrasts(Intercept, gTypeH898res, txWound, txGallery,
                gTypeH898res_txWound, gTypeH898res_txGallery,
                weevil = gTypeH898res_txGallery - gTypeH898res_txWound,
                wound = gTypeH898res + gTypeH898res_txWound,
                gallery = gTypeH898res + gTypeH898res_txGallery,
                levels = modMat)
fit2 <- contrasts.fit(fit, cont_matrix)
fit3 <- eBayes(fit2)

# get inferential summary for the six terms in the model

# to get individual t statistics, etc., must use topTable() on each coef
# separately
model_terms <- colnames(modMat)
names(model_terms) <- model_terms
statInf_model_terms <-
  adply(model_terms, 1,
        function(x) name_rows(topTable(fit3, coef = x,
                                       number = Inf,sort.by = "none")))
statInf_model_terms <-
  rename(statInf_model_terms, c(X1 = "model_term", .rownames = "contig"))
statInf_model_terms <-
  statInf_model_terms[ ,c("contig", setdiff(names(statInf_model_terms),"contig"))]
head(statInf_model_terms)
str(statInf_model_terms)

## informal tests so we know if things change, differ for Jenny vs Mack, etc.
test_that("stat inf on the model terms has 393654 rows",
          expect_equal(393654, nrow(statInf_model_terms)))
(t_medians <- aggregate(t ~ model_term, statInf_model_terms, median))
all.equal(t_medians$t, c(1.118479141, 0.157361758, 0.092882336,
                         0.026698663, -0.067348886, 0.005007678))

write.table(statInf_model_terms,
            "limma-results-model-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# get inferential summary for the effects we are most interested in

# below, we will rename gTypeH898res to control, to be parallel with would and
# gallery
focus_terms <- c("weevil", "gTypeH898res", "wound", "gallery")
names(focus_terms) <- focus_terms
statInf_focus_terms <-
  adply(focus_terms, 1,
        function(x) name_rows(topTable(fit3, coef = x,
                                       number = Inf, sort.by = "none")))
statInf_focus_terms <-
  rename(statInf_focus_terms, c(X1 = "focus_term", .rownames = "contig"))
statInf_focus_terms$focus_term <-
  revalue(statInf_focus_terms$focus_term, c(gTypeH898res = "control"))

statInf_focus_terms <-
  statInf_focus_terms[ , c("contig",
                           setdiff(names(statInf_focus_terms),"contig"))]
head(statInf_focus_terms)
str(statInf_focus_terms)

## informal tests so we know if things change, differ for Jenny vs Mack, etc.
test_that("stat inf on the focus terms has 262436 rows",
          expect_equal(262436, nrow(statInf_focus_terms)))
(t_medians <- aggregate(t ~ focus_term, statInf_focus_terms, median))
all.equal(t_medians$t, c(0.04817940, 0.157361758, 0.01198738, 0.10723831))
# "Mean relative difference: 2.508257e-08" <-- that's OK!

write.table(statInf_focus_terms,
            "limma-results-focus-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

