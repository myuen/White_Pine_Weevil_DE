#+ setup, include = FALSE
library(knitr)
library(plyr)
library(testthat) 
library(ggplot2)
library(reshape2)
library(fdrtool) # fdrtool() to compute pi0 (or eta0, as they say)
opts_chunk$set(fig.path = 'figure/04_explore-dea-results-')

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

#' Explore the distribution of the estimated "weevil" effect of interest, i.e.
#' the difference between the interaction terms

#+ weevil-estimates, fig.show='hold', out.width='50%', message = FALSE
p <- ggplot(subset(sidf, focus_term == "weevil"), aes(x = logFC))
p + geom_histogram()
p + geom_density()

#' Explore the distribution of estimates for the various effects of interest.
#' Note: one facet will be blank, corresponding to the test whether all terms
#' involving genotype are equal to zero, since there is no single estimate that
#' is appropriate to show.

#+ focus-term-estimates, fig.show='hold', out.width='50%', message = FALSE
p <- ggplot(sidf, aes(x = logFC))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Explore the distribution of t statistics for the various effects of interest.
#' Note: one facet will be blank, corresponding to the test whether all terms 
#' involving genotype are equal to zero, since that is addressed by an F
#' statistic instead.

#+ focus-term-t-statistics, fig.show='hold', out.width='50%', message = FALSE
p <- ggplot(sidf, aes(x = t))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Explore the distribution of p-values for the various effects of interest. 
#' Note: this works same for all effects of interest.

#' For each focus term, estimate pi0 = the proportion of null contigs (note: the
#' authors of the fdrtool package use the nonstandard term eta0)
(pi0 <- ddply(sidf, ~ focus_term, function(z) {
  fdr_result <- fdrtool(z$P.Value, statistic = "pvalue",
                        plot = FALSE, verbose = FALSE)
  return(fdr_result$param[1, "eta0"])  
}))

#+ focus-term-p-values, fig.show='hold', out.width='50%', message = FALSE
p <- ggplot(sidf, aes(x = P.Value))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Explore the distribution of Benjamini-Hockberg adjusted p-values for the
#' various effects of interest. Note: this works same for all effects of
#' interest.

#+ focus-term-adjusted-p-values, fig.show='hold', out.width='50%', message = FALSE
p <- ggplot(sidf, aes(x = adj.P.Val))
p + geom_histogram() + facet_wrap( ~ focus_term)
p + geom_density(aes(colour = focus_term))

#' Source function to extract and tidy count data
source("helper04_extract-and-tidy.r")

#' Find the top hits w/r/t genotype effects
y <- subset(sidf, focus_term == "gTypeH898res_all")
hit_sidf_row <- with(y, which(rank(P.Value) < 5))
y[hit_sidf_row, ]
hit_contig <- y$contig[hit_sidf_row]

jDat <- extract_and_tidy(hit_contig, x, expDes)
str(jDat)

p <- ggplot(jDat, aes(x = count, y = gType, colour = tx))
p + geom_jitter(position = position_jitter(height = .15)) +
  facet_wrap(~ contig, scales="free_x")

#' Find the top hits w/r/t the weevil effect
y <- subset(sidf, focus_term == "weevil")
hit_sidf_row <- with(y, which(rank(P.Value) < 7))
y[hit_sidf_row, ]
hit_contig <- y$contig[hit_sidf_row]

jDat <- extract_and_tidy(hit_contig, x, expDes)
str(jDat)

p <- ggplot(jDat, aes(x = count, y = gType, colour = tx))
p + geom_jitter(position = position_jitter(height = .15)) +
  facet_wrap(~ contig, scales="free_x")

#' Find the top hits w/r/t genotype in the control condition specifically
y <- subset(sidf, focus_term == "gTypeH898res_in_control")
hit_sidf_row <- with(y, which(rank(P.Value) < 10))
y[hit_sidf_row, ]
hit_contig <- y$contig[hit_sidf_row]

jDat <- extract_and_tidy(hit_contig, x, expDes)
str(jDat)

p <- ggplot(subset(jDat, tx == "Control"),
            aes(x = count, y = gType))
p + geom_jitter(position = position_jitter(height = .15)) +
  facet_wrap(~ contig, scales="free_x")
