library(ggplot2)
library(grid)

setwd("analysis/")

source("helper03_load-focus-statinf.r")
source("helper06_classifyCtgs.R")

sidf <- load_focus_statInf()
str(sidf)

lfc <- 2
pCutoff <- 0.01

## log FC plot of ctrl_vs_ctrl(H898C - Q903C) against combined_effect_Q903 (Q903G - Q903W)

# subset data, in the order specificed by the user
mDat <- subset(sidf, sidf$focus_term == "constDiff")
mDat <- rbind(mDat, subset(sidf, sidf$focus_term == "weevilInd_Q903"))
summary(mDat$focus_term)

mDat_wide <- reshape(mDat, direction = "wide", 
                     timevar = "focus_term", idvar = "contig",
                     drop = c("F", "P.Value", "B", "AveExpr", "t"))
mDat_wide <- data.frame(mDat_wide, row.names = 1)


# mDat_wide col order : 
# 1) logFC of 1st focus, 2) adj.P.Val of 1st focus, 
# 3) logFC of 2nd focus, 4) adj.P.Val of 2nd focus

# Call function to classify contig based on stat. sig and log fold-change
mDat_wide$type <- apply(mDat_wide, 1, classifyCtgs)

# Re-order levels for presentation in graph
mDat_wide$type <- 
  factor(mDat_wide$type, levels = c('UpUp', 'UpDown', 'DownUp', 'DownDown', 
                                    'SigBoring', 'NS','NSNS'))

summary(mDat_wide$type)
# UpUp    UpDown    DownUp  DownDown SigBoring        NS      NSNS 
#   40        86        59       825       666     14490     41938 


# Plot
q <- ggplot(subset(mDat_wide, mDat_wide$type == "NS" | mDat_wide$type == "NSNS"),
            aes(x = logFC.constDiff, y = logFC.weevilInd_Q903)) + 
  geom_point(color = "gray90", size = 1, alpha = 0.75) + 
  geom_point(data = subset(mDat_wide, mDat_wide$type != "NS" & mDat_wide$type != "NSNS"),
             aes(x = logFC.constDiff, y = logFC.weevilInd_Q903, color = type),
             size = 1, alpha = 0.9) + 
  scale_colour_manual(
    # Change color of label
    values = c("#b2182b", "#ef8a62", "#67a9cf", "#2166ac", "gray40"), 
    # Change legend label
    labels = 
      c("Stat. Sig. and Up Reg. in Both Contrasts",
        "Stat. Sig. Up Reg. in Constitutive Diff. and \nDown Reg. in Q903 Weevil Induced", 
        "Stat. Sig. Down Reg. in Constitutive Diff. and \nUp Reg. in Q903 Weevil Induced", 
        "Stat. Sig. and Down Reg. in Both Contrasts",
        "Stat. Sig but not Sig. Diff. Expressed"), 
    guide = guide_legend(title = NULL)) +
  labs(title = "Expression Fold Change for Control to Control Comparison against 
       Combined Wounding and Feeding Effect in Q903", 
       x = "Log 2 Fold Change for Constitutive Difference", 
       y = "Log 2 Fold Change for Weevil Induced in Q903") + 
  theme_bw() + 
  theme(plot.title = element_text(size = rel(2)),
        legend.position = c(0.88, 0.15),
        legend.key.height = unit(0.85, "cm")) + 
  guides(colour = guide_legend(title = NULL, override.aes = list(size = 2))) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_abline(aes(intercept = 0, slope = 0), colour = "black")


ggsave("../results/figures/constDiff_To_Q903WeevilInd.png", 
       plot = q, width = 16, height = 9)
