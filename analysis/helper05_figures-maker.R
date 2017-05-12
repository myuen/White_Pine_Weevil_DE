require (ggplot2)

## function to make Volcano plot
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
    labs(title = title, x = "Log 2 Fold Change", y = "-log 10 (Adjusted P Value)") + 
    scale_x_continuous(limits = c(-20, 20)) + theme_bw() + 
    theme(plot.title = element_text(size = rel(2)), legend.position = "none")
  
  if (nrow(sig) > 0) {
    # Add second layer with significant points in red (i.e. higher than logFC cutoff
    # with stat significant p-value)
    g <- g + geom_point(data = sig, aes(logFC, -log10(adj.P.Val), color = "red"), size = 0.7)
  }

  return (g)
}


## function to plot logFC against Average Expression
plotFC2AvgExpr <- function(df, focus, lfcCutoff, pCutoff) {
  g <- ggplot(subset(df, focus_term == focus & adj.P.Val > pCutoff), 
              aes(x = AveExpr, y = logFC)) + 
    # stat insignificant data points
    geom_point(shape = 20, color = "#666666", size = 1) + 
    # stat. sig. and up regulated data points
    geom_point(data = subset(df, focus_term == focus & 
                               adj.P.Val <= pCutoff & logFC >= lfcCutoff),
               aes(x = AveExpr, y = logFC),
               colour = "red", fill = "red", shape = 24, size = 0.8) + 
    # stat. sig. and down regulated data points
    geom_point(data = subset(df, df$focus_term == focus & 
                               adj.P.Val <= pCutoff & logFC <= (-1 * lfcCutoff)),
               aes(x = AveExpr, y = logFC), 
               colour = "blue", fill= "blue", shape = 25, size = 0.8) + 
    geom_abline(aes(intercept = 0, slope = 0), colour = "black") + 
    scale_y_continuous(limits = c(-20, 20)) + theme_bw() +
    labs(title = paste0("logFC against Expression in ", focus), 
         x = "Log 2 Average Expression Counts", y = "Log 2 Fold Changes") +
    theme(plot.title = element_text(size = rel(2)))
  return (g)
}
