library(plyr)
library(testthat) 
library(ggplot2)
library(reshape2)

#' read and validate the statistical inference results
sidf <- read.delim("limma-results-focus-terms.tsv")
str(sidf) # 'data.frame':  262436 obs. of  8 variables
test_that("inference results for our focus terms (still) have 262436 rows",
          expect_equal(262436, nrow(sidf)))
sidf$focus_term <- factor(sidf$focus_term,
                          levels = c("weevil", "control", "wound", "gallery"))
(t_medians <- aggregate(t ~ focus_term, sidf, median))
test_that("medians of the t statistics are what we expect",
          expect_equal(t_medians$t,
                       c(0.04817940, 0.157361758, 0.01198738, 0.10723831),
                       tolerance = 2 * .Machine$double.eps ^ 0.5))

#' read and validate the "raw" count data
filteredSailfishCounts <- # take a few moments
  read.delim("consolidated-filtered-Sailfish-results.txt")
str(filteredSailfishCounts, list.len = 8)
test_that("filtered Sailfish data has 65609 rows upon import",
          expect_equal(65609, nrow(filteredSailfishCounts)))
test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(filteredSailfishCounts)))

#' read and validate the  experimental design
expDes <- read.delim("White_Pine_Weevil_exp_design.tsv",
                     stringsAsFactors = FALSE)
expDes <-
  mutate(expDes,
         gType = factor(gType, levels = c("Q903susc", "H898res")),
         txCode = factor(txCode, levels = c('C', 'W', 'G')),
         tx = factor(tx, levels = c("Control", "Wound", "Gallery")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:
test_that("design matrix has 24 rows upon import",
          expect_equal(24, nrow(expDes)))

#+ weevil-estimates
p <- ggplot(subset(sidf, focus_term == "weevil"), aes(x = logFC))
p + geom_histogram()
p + geom_density()

#+ focus-term-estimates
p <- ggplot(sidf, aes(x = logFC))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#+ focus-term-t-statistics
p <- ggplot(sidf, aes(x = t))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#+ focus-term-p-values
p <- ggplot(sidf, aes(x = P.Value))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#+ focus-term-adjusted-p-values
p <- ggplot(sidf, aes(x = adj.P.Val))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' explore raw data, guided by statistical results
#' 
#' reshape raw count data, i.e. tidy it
jDat <- cbind(expDes, t(filteredSailfishCounts))
names(jDat) <- gsub("WPW_Inoculation_Trinity_C500_", "", names(jDat))
jDat <- melt(jDat,
             id.vars = names(expDes),
             variable.name = 'contig', value.name = 'count')
str(jDat) # 'data.frame':  1574616 obs. of  8 variables:

#' get the contig IDs guided by stat inf results re: weevil effect
weevInf <- droplevels(subset(sidf, focus_term == "weevil", select = -focus_term))
str(weevInf)

weevInf <- arrange(weevInf, P.Value)
head(weevInf)

feat_contigs <- weevInf$contig[1:9]
feat_contigs <- gsub("WPW_Inoculation_Trinity_C500_", "", feat_contigs)

pDat <- subset(jDat, contig %in% feat_contigs)

p <- ggplot(pDat,
            aes(x = count, y = gType, colour = tx))
p + geom_point() + facet_wrap(~ contig, scales="free_x")
