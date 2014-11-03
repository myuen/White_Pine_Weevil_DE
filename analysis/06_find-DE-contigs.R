library(plyr)

setwd("analysis/")
source("helper03_load-focus-statinf.r")


topN <- 200
lfc <- 2
pCutoff <- 0.01


sidf <- load_focus_statInf()
str(sidf)

sidf_wide <- reshape(sidf, direction = "wide", timevar = "focus_term", idvar = "contig",
                     drop = c("t", "F", "P.Value", "B"))
str(sidf_wide)


# function to contigs that are common between two given focus_terms
findCommonContigs <- function(df, pCutoff, lfc, x, y) {
  subset(df, eval(parse(
    text = paste0("df$adj.P.Val.", x, " <= ", pCutoff, " & ",
                 "df$adj.P.Val.", y, " <= ", pCutoff, " & ",
                 "abs(df$logFC.", x, ") > ", lfc, " & ",
                 "abs(df$logFC.", y, ") > ", lfc)))) 
}


# Common contigs between feeding in Q903 and combined in Q903
dim(findCommonContigs(sidf_wide, pCutoff, lfc, "feed_Q903", "combined_effect_Q903"))
# [1] 5132   40


# Common contigs between ctrl_vs_ctrl and cmpd in Q903
dim(findCommonContigs(sidf_wide, pCutoff, lfc, "ctrl_vs_ctrl", "combined_effect_Q903"))
# [1] 1067   40


# Count all differentially expressed contigs for each focus_term
sig_DE_counts <- ddply(sidf, ~ focus_term, function(x) {
  up <- length(subset(x, x$adj.P.Val <= 0.01 & x$logFC > 0)[,"contig"]);
  down <- length(subset(x, x$adj.P.Val <= 0.01 & x$logFC < 0)[,"contig"]);
  sum <- up + down
  return (c("up" = up, "down" = down, "sum" = sum));
})

#              focus_term   up down   sum
# 1            wound_Q903    0    1     1
# 2            wound_H898    0    0     0
# 3             feed_Q903 3244 6272  9516
# 4             feed_H898    0    4     4
# 5  combined_effect_Q903 5366 5332 10698
# 6  combined_effect_H898    1    2     3
# 7          ctrl_vs_ctrl 4897 4496  9393
# 8        wound_vs_wound 4660 4638  9298
# 9    gallery_vs_gallery 7788 5956 13744
# 10        wounding_diff    0    0     0
# 11         feeding_diff  693  144   837
# 12        combined_diff  293  208   501
# 13     induced_vs_const 3950 6705 10655
