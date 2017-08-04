### figure 6c + 6d

library(ggplot2)

setwd("analysis/")

source("helper03_load-focus-statinf.r")

sidf <- load_focus_statInf()
str(sidf)

# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfcCutoff) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05

cd <- subset(sidf, sidf$focus_term == "constDiff")
wiq <- subset(sidf, sidf$focus_term == "weevilInd_Q903")

minX <- min(floor(c(cd$logFC, wiq$logFC)))
maxX <- max(ceiling(c(cd$logFC, wiq$logFC)))

maxY <- ceiling(max(-log10(c(cd$adj.P.Val, wiq$adj.P.Val))))

# Function to make Volcano plot
volcanoPlot <- function(df, lfcCutoff, pCutoff, title) {
  
  # Plot insignificant points in black (i.e. lower than logFC cutoff 
  # and/or insignificant p-value)
  insig = subset(df, (abs(df$logFC) < lfcCutoff | df$adj.P.Val > pCutoff))
  sig = subset(df, (abs(logFC) >= lfcCutoff & adj.P.Val <= pCutoff))
  
  g <- ggplot(data = insig, aes(x = logFC, y = -log10(adj.P.Val))) + 
    geom_point(colour = "black", size = 0.7) + 
    
    geom_abline(aes(intercept = -log10(pCutoff), slope = 0), colour = "blue") + 
    geom_vline(xintercept = lfcCutoff, colour = "blue") + 
    geom_vline(xintercept = -(lfcCutoff), colour = "blue") +
    labs(title = title,
         x = expression(paste('Log'[2], " Fold Change")), 
         y = expression(paste('-Log'[10], " (Adjusted p-Value)"))) + 
    scale_x_continuous(limits = c(minX, maxX)) + 
    scale_y_continuous(limits = c(0, maxY)) + 
    theme_bw() + 
    theme(plot.title = element_text(size = rel(2)), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (nrow(sig) > 0) {
    # Add second layer with significant points in red (i.e. higher than logFC cutoff
    # with stat significant p-value)
    g <- g + geom_point(data = sig, aes(logFC, -log10(adj.P.Val), color = "red"), size = 0.7)
  }
  
  return (g)
}


fig6c <- volcanoPlot(cd, lfcCutoff, pCutoff, "(c)")
ggsave("../results/figures/fig6c.4Aug2017.tiff", plot = fig6c, height = 3.5, width = 4)

fig6d <- volcanoPlot(wiq, lfcCutoff, pCutoff, "(d)")
ggsave("../results/figures/fig6d.4Aug2017.tiff", plot = fig6d, height = 3.5, width = 4)
