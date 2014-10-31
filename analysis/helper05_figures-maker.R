require (ggplot2)

## function to make Volcano plot
volcanoPlot <- function(df, lfc, pCutoff, title) {
  require (ggplot2)
  g <- ggplot(data = subset(df, (abs(df$logFC) < lfc | df$adj.P.Val > pCutoff)), 
              aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(colour = "black", size = 1) + 
    geom_point(data = subset(df, (abs(logFC) >= lfc & adj.P.Val <= pCutoff)), 
               aes(logFC, -log10(adj.P.Val), color = "red"), size = 1) + 
    geom_abline(aes(intercept = -log10(pCutoff), slope = 0), colour = "blue") + 
    geom_vline(xintercept = lfc, colour = "blue") + 
    geom_vline(xintercept = -(lfc), colour = "blue") +
    labs(title = title, x = "Log 2 Fold Change", y = "-log 10 (Adjusted P Value)") + 
    scale_x_continuous(limits = c(-20, 20)) + theme_bw() + 
    theme(plot.title = element_text(size = rel(2)), legend.position = "none")
  return (g)
}


## function to plot logFC against Average Expression
plotFC2AvgExpr <- function(df, focus, lfc, pCutoff) {
  g <- ggplot(subset(df, focus_term == focus & adj.P.Val > pCutoff), 
              aes(x = 2 ^ AveExpr, y = logFC)) + geom_point(shape = 20, color = "#666666") + 
    # up regulated data points
    geom_point(data = subset(df, focus_term == focus & 
                               adj.P.Val <= pCutoff & logFC > lfc),
               aes(x = 2 ^ AveExpr, y = logFC), colour = "red", shape = 17) + 
    # down regulated data points
    geom_point(data = subset(df, df$focus_term == focus & 
                               adj.P.Val <= pCutoff & logFC < lfc),
               aes(x = 2 ^ AveExpr, y = logFC), colour = "red", fill= "red", shape = 25) + 
    geom_abline(aes(intercept = 0, slope = 0), colour = "black") + 
    scale_x_log10(breaks=c(0.1, 1, 10, 100, 1000), labels=c(0.1, 1, 10, 100, 1000)) + 
    scale_y_continuous(limits = c(-20, 20)) + theme_bw() +
    labs(title = paste0("logFC against Expression in ", focus), 
         x = "Average Expression Counts", y = "Log 2 Fold Changes") +
    theme(plot.title = element_text(size = rel(2)))
  return (g)
}
