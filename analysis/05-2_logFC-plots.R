# Plot logFC of multiple contrasts

library(ggplot2)
library(grid)

setwd("analysis/")

source("helper03_load-focus-statinf.r")
source("helper06_classifyCtgs.R")

sidf <- load_focus_statInf()
str(sidf)

# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05

# List contrast to compare
cc <- c("constDiff", "weevilInd_Q903")

# subset data
mDat <- subset(sidf, sidf$focus_term %in% cc)
summary(mDat$focus_term)


mDat_wide <- reshape(mDat, direction = "wide", 
                     timevar = "focus_term", idvar = "cds",
                     drop = c("F", "P.Value", "B", "AveExpr", "t", "contig"))
mDat_wide <- data.frame(mDat_wide, row.names = 1)


# Call function to classify contig based on stat. sig and log fold-change
# mDat_wide$type <- apply(mDat_wide, 1, classifyCtgs)
mDat_wide$type <- apply(mDat_wide, 1, 
                        function(x) classifyCtgs(x, lfcCutoff, pCutoff))


# Re-order levels for presentation in graph
mDat_wide$type <- 
  factor(mDat_wide$type, levels = c('UpUp', 'UpDown', 'DownUp', 'DownDown', 
                                    'SigBoring', 'NS','NSNS'))

summary(mDat_wide$type)
#      UpUp    UpDown    DownUp  DownDown SigBoring        NS      NSNS 
#        87       106        65       641      1385     13338     22985 


# Plot
q <- ggplot(subset(mDat_wide, mDat_wide$type == "NS" | mDat_wide$type == "NSNS"),
            aes_string(x = colnames(mDat_wide)[1], y = colnames(mDat_wide)[3])) + 
  geom_point(color = "gray90", size = 1, alpha = 0.75) + 
  geom_point(data = subset(mDat_wide, mDat_wide$type != "NS" & mDat_wide$type != "NSNS"),
             aes_string(x = colnames(mDat_wide)[1], y = colnames(mDat_wide)[3], 
                        color = "type"), size = 1, alpha = 0.9) + 
  scale_colour_manual(
    # Change color of label
    values = c("#b2182b", "#ef8a62", "#67a9cf", "#2166ac", "gray40"), 
    # Change legend label
    labels = 
      c("Stat. Sig. and Up Regulation in Both Contrasts",
        paste0("Stat. Sig. Up Regulation in ", cc[1] , " and \nDown Regulation in ", cc[2]), 
        paste0("Stat. Sig. Down Regulation in ", cc[1], " and \nUp Regulation in ", cc[2]), 
        "Stat. Sig. and Down Regulation in Both Contrasts",
        "Stat. Sig. but not Sig. Diff. Expressed"), 
    guide = guide_legend(title = NULL)) +
  labs(title = paste0("Expression Fold Change for ", cc[1], " Compare against ", cc[2]), 
       x = paste0("Log 2 Fold Change for ", cc[1]), 
       y = paste0("Log 2 Fold Change for ", cc[2])) +  
  theme_bw() + 
  theme(plot.title = element_text(size = rel(2)),
        legend.position = c(0.88, 0.15),
        legend.key.height = unit(0.85, "cm")) + 
  guides(colour = guide_legend(title = NULL, override.aes = list(size = 2))) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_abline(aes(intercept = 0, slope = 0), colour = "black")

ggsave(paste0("../results/figures/", cc[1], "_To_", cc[2], ".31May2016.png"), 
       plot = q, width = 16, height = 9)
