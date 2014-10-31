library(ggplot2)

setwd("analysis/")

# source("helper01_load-counts.r")
# source("helper02_load-exp-des.r")
source("helper03_load-focus-statinf.r")
# source("helper04_extract-and-tidy.r")
# source("helper05_figures-maker.R")
source("helper06_classifyCtgs.R")


sidf <- load_focus_statInf()
str(sidf)


topN <- 200
lfc <- 2
pCutoff <- 0.01


## log FC plot of ctrl_vs_ctrl(H898C - Q903C) against combined_effect_Q903 (Q903G - Q903C)

# subset data, in the order specificed by the user
mDat <- subset(sidf, sidf$focus_term == "ctrl_vs_ctrl")
mDat <- rbind(mDat, subset(sidf, sidf$focus_term == "combined_effect_Q903"))
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
# 120        50        86       805       673     16538     47337


# Plot
q <- ggplot(subset(mDat_wide, mDat_wide$type == "NS" | mDat_wide$type == "NSNS"),
            aes(x = logFC.ctrl_vs_ctrl, y = logFC.combined_effect_Q903)) + 
  geom_point(color = "gray80", size = 2, alpha = 0.75) + 
  geom_point(data = subset(mDat_wide, mDat_wide$type != "NS" & mDat_wide$type != "NSNS"),
             aes(x = logFC.ctrl_vs_ctrl, y = logFC.combined_effect_Q903, color = type),
             size = 2, alpha = 0.9) + 
  scale_colour_manual(
    # Change color of label
    values = c("#d7191c", "#fdae61", "#2c7bb6", "#a6d96a", "gray50"), 
    # Change legend label
    labels = c("Stat. Sig. and Up Reg. in Both Contrasts",
               "Stat. Sig. Up Reg. in Ctrl vs. Ctrl and Down Reg. in Q903 Feeding", 
               "Stat. Sig. Down Reg. in Ctrl vs. Ctrl and Up Reg. in Q903 Feeding", 
               "Stat. Sig. and Down Reg. in Both Contrasts",
               "Stat. Sig but not Sig. Diff. Expressed"), 
    guide = guide_legend(title = NULL)) +
  labs(title = "Expression Fold Change for Control to Control Comparison against 
       Combined Wounding and Feeding Effect in Q903", 
       x = "Log 2 Fold Change for Control vs Control", 
       y = "Log 2 Fold Change for Feeding and Wounding Effect in Q903") + 
  theme_bw() + theme(plot.title = element_text(size = rel(2)), 
                     legend.position = c(0.8, 0.2)) + 
  geom_vline(xintercept = 0, colour = "black") +
  geom_abline(aes(intercept = 0, slope = 0), colour = "black")


ggsave("../results/figures/ctrlVsCtrl_To_Q903FeedAndWound.png", 
       plot = q, width = 16, height = 9)
