library(ggplot2)

setwd("analysis/")

source("helper03_load-focus-statinf.r")
source("helper05_figures-maker.R")


sidf <- load_focus_statInf()
str(sidf)

sidf_wide  <- reshape(sidf, direction = "wide", 
                      timevar = "focus_term", idvar = c("contig"), 
                      drop = c("F", "P.Value", "B", "AveExpr", "t"))
str(sidf_wide)


topN <- 200
lfc <- 2
pCutoff <- 0.01


# Volcano plot for (Q903G - Q903C) - (H898C - Q903C)

# needs a better title for plot
ivc_vp <- volcanoPlot(
  df = subset(sidf, sidf$focus == "induced_vs_const"), 
  lfc, pCutoff, "Volcano Plot for Q903C Inducible Expression to Control Constitutive Expression")

ggsave("../results/figures/inducible_vs_constitutive_vPlot.png", plot = ivc_vp,
       height = 9, width = 16)


# Log 2 Fold Change against Average Expression Counts
ivc_logFC2avgExpr <- plotFC2AvgExpr(sidf, "induced_vs_const", lfc, pCutoff)
ggsave("../results/figures/inducible_vs_constitutive_logFC2avgExpr.png", 
       plot = ivc_logFC2avgExpr, height = 9, width = 16)
