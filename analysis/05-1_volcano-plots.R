library(ggplot2)
library(purrr)

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


# You should be able to pass focus_terms <- levels(sidf$focus_term) to lapply
# but one of our focus_term do not produce any sig. DE contigs and will crash ggplot

#        focus_term   up down  sum
# 1       constDiff 2279 2242 4521
# 2  woundResp_Q903    0    1    1
# 3  woundResp_H898    0    0    0
# 4  weevilInd_Q903 2000 3763 5763
# 5  weevilInd_H898    0    0    0
# 6 weevilCtrl_Q903 3319 3925 7244
# 7 weevilCtrl_H898    2    1    3

focus_terms <- c("constDiff", "weevilInd_Q903", "weevilCtrl_Q903")

# lapply(focus_terms, function(x) {
#   p <- volcanoPlot(subset(sidf, sidf$focus_term == x), lfcCutoff, pCutoff, x);
#   
#   ggsave(paste0("../results/figures/", x, "_vPlot.31May2016.png"), 
#          plot = p, height = 8.5, width = 11)
#   }
# )

map(focus_terms, function(x) {
  p <- volcanoPlot(subset(sidf, sidf$focus_term == x), lfcCutoff, pCutoff, x);
  
  ggsave(paste0("../results/figures/", x, "_vPlot.24Aug2016.png"), 
         plot = p, height = 8.5, width = 11)
}
)
