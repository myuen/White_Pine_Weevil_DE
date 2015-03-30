library(ggplot2)

setwd("analysis/")

source("helper03_load-focus-statinf.r")
source("helper05_figures-maker.R")

sidf <- load_focus_statInf()
str(sidf)

lfc <- 2
pCutoff <- 0.01

# woundResp_H898 produce no diff exp and exluded from list
focus_terms <- c("constDiff", "woundResp_Q903", 
                  "weevilInd_Q903", "weevilInd_H898", 
                  "weevilCtrl_Q903", "weevilCtrl_H898")

lapply(focus_terms, function(x) {
  p <- volcanoPlot(subset(sidf, sidf$focus_term == x), lfc, pCutoff, x);
  
  ggsave(paste0("../results/figures/", x, "_vPlot.png"), 
         plot = p, height = 8.5, width = 11)
  
  }
)
