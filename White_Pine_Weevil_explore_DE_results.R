#+ setup, include = FALSE
library(knitr)
library(plyr)
library(testthat) 
library(ggplot2)
library(reshape2)
##opts_chunk$set(fig.path = 'figure/stripplot-')

#' Source purpose-built functions to load and validate the data, the
#' experimental design, and the statistical inference results for our focus
#' terms. Then call them.
source("helper01_load-counts.r")
source("helper02_load-exp-des.r")
source("helper03_load-focus-statinf.r")

x <- load_counts()
str(x, list.len = 8) # 'data.frame':  65609 obs. of  24 variables:
expDes <- load_expDes()
str(expDes) # 'data.frame':  24 obs. of  6 variables:
sidf <- load_focus_statInf()
str(sidf) # 'data.frame':  328045 obs. of  9 variables:

#' Explore the distribution of the "weevil" effect of interest, i.e. the
#' difference between the interaction terms

#+ weevil-estimates, fig.show='hold', out.width='50%'
p <- ggplot(subset(sidf, focus_term == "weevil"), aes(x = logFC))
p + geom_histogram()
p + geom_density()

#' Explore the distribution of estimates for the various effects of interest.
#' Note: one facet will be blank, corresponding to the test whether all terms
#' involving genotype are equal to zero, since there is no single estimate that
#' is appropriate to show.

#+ focus-term-estimates, fig.show='hold', out.width='50%'
p <- ggplot(sidf, aes(x = logFC))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Explore the distribution of t statistics for the various effects of interest.
#' Note: one facet will be blank, corresponding to the test whether all terms 
#' involving genotype are equal to zero, since that is addressed by an F
#' statistic instead.

#+ focus-term-t-statistics, fig.show='hold', out.width='50%'
p <- ggplot(sidf, aes(x = t))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Explore the distribution of p-values for the various effects of interest. 
#' Note: this works same for all effects of interest.

#+ focus-term-p-values, fig.show='hold', out.width='50%'
p <- ggplot(sidf, aes(x = P.Value))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Explore the distribution of Benjamini-Hockberg adjusted p-values for the
#' various effects of interest. Note: this works same for all effects of
#' interest.

#+ focus-term-adjusted-p-values, fig.show='hold', out.width='50%'
p <- ggplot(sidf, aes(x = adj.P.Val))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

## I AM HERE

#' Find the top hit w/r/t genotype effects
hit_sidf_row <- with(sidf,
                     which(focus_term == "gTypeH898res_all" & rank(P.Value) < 5))
sidf[hit_sidf_row, ]
hit_contig <- sidf$contig[hit_sidf_row]

hit_data_row <- which(rownames(x) %in% hit_contig)
x[hit_data_row, ]

jDat <- cbind(expDes, t(x[hit_data_row, ]))
names(jDat) <- gsub("WPW_Inoculation_Trinity_C500_", "", names(jDat))
jDat <- melt(jDat,
             id.vars = names(expDes),
             variable.name = 'contig', value.name = 'count')
str(jDat) # 'data.frame':  1574616 obs. of  8 variables:

p <- ggplot(jDat, aes(x = count, y = gType, colour = tx))
p + geom_point() + facet_wrap(~ contig, scales="free_x")



#' explore raw data, guided by statistical results
#' 
#' reshape raw count data, i.e. tidy it
# jDat <- cbind(expDes, t(filteredSailfishCounts))
# names(jDat) <- gsub("WPW_Inoculation_Trinity_C500_", "", names(jDat))
# jDat <- melt(jDat,
#              id.vars = names(expDes),
#              variable.name = 'contig', value.name = 'count')
# str(jDat) # 'data.frame':  1574616 obs. of  8 variables:
# 
# #' get the contig IDs guided by stat inf results re: weevil effect
# weevInf <- droplevels(subset(sidf, focus_term == "weevil", select = -focus_term))
# str(weevInf)
# 
# weevInf <- arrange(weevInf, P.Value)
# head(weevInf)
# 
# feat_contigs <- weevInf$contig[1:9]
# feat_contigs <- gsub("WPW_Inoculation_Trinity_C500_", "", feat_contigs)
# 
# pDat <- subset(jDat, contig %in% feat_contigs)
# 
# p <- ggplot(pDat,
#             aes(x = count, y = gType, colour = tx))
# p + geom_point() + facet_wrap(~ contig, scales="free_x")
