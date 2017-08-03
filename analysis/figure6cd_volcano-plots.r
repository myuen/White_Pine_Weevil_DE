### figure 6c + 6d

library(ggplot2)

setwd("analysis/")

source("helper03_load-focus-statinf.r")
source("helper05_figures-maker.R")

sidf <- load_focus_statInf()
str(sidf)

# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfcCutoff) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05


fig6c <- volcanoPlot(subset(sidf, sidf$focus_term == "constDiff"), lfcCutoff, pCutoff, "(c)")
ggsave("../results/figures/fig6c.3Aug2017.tiff", plot = fig6c, height = 3.5, width = 4)

fig6d <- volcanoPlot(subset(sidf, sidf$focus_term == "weevilInd_Q903"), lfcCutoff, pCutoff, "(d)")
ggsave("../results/figures/fig6d.3Aug2017.tiff", plot = fig6d, height = 3.5, width = 4)
