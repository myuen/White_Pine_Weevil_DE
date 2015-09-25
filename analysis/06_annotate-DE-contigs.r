library(plyr)

setwd("analysis/")
source("helper03_load-focus-statinf.r")

sidf <- load_focus_statInf()
str(sidf)

lfc <- 2
pCutoff <- 0.01


# Count all differentially expressed contigs for each focus_term
ddply(sidf, ~ focus_term, function(x) {
  up <- length(subset(x, x$adj.P.Val <= pCutoff & x$logFC >= lfc)[,"contig"]);
  down <- length(subset(x, x$adj.P.Val <= pCutoff & x$logFC <= -1 * lfc)[,"contig"]);
  sum <- up + down;
  return (c("up" = up, "down" = down, "sum" = sum));
})

#        focus_term   up down  sum
# 1       constDiff 3635 3377 7012
# 2  woundResp_Q903    0    1    1
# 3  woundResp_H898    0    0    0
# 4  weevilInd_Q903 1481 4321 5802
# 5  weevilInd_H898    0    2    2
# 6 weevilCtrl_Q903 2661 3618 6279
# 7 weevilCtrl_H898    0    2    2


# BLAST results on tab-delimited format.  We get at most 10 hits per query sequences.
annots <- read.delim("../results/WPW_Inoculation_Trinity_C500.diffExp.lfc2.blastxNr.txt", 
                     sep = "\t", header = TRUE)
dim(annots) # [1] 129980     15


# Function to concatenate multiple blast annontations for each query sequence
concatAnnot <- function(x) {
  annot <- x$salltitles
  annot <- gsub("\\<\\>", " ", annot, perl = TRUE)
  return (paste(annot, collapse = " ; "))
}


annotsCollapsed <-
  ddply(annots, ~ qseqid, concatAnnot)
colnames(annotsCollapsed) <- c("contig", "annot")


write.table(annotsCollapsed,
            "../results/WPW_Inoculation_Trinity_C500.diffExp.lfc2.blastxNr.collapsed.txt",
  sep = "\t", quote = FALSE, row.names = FALSE)


# Find all DE contigs
sigDE <- 
  ddply(sidf, ~ focus_term, 
        function(x) {
          tmp_sigDE <- subset(x, (abs(logFC)) >= lfc & adj.P.Val <= pCutoff);
          tmp_sigDE <- tmp_sigDE[order(tmp_sigDE$logFC, decreasing = TRUE), ];
        })


focus_terms <- c("constDiff", "woundResp_Q903", "woundResp_H898",
                 "weevilInd_Q903", "weevilInd_H898", 
                 "weevilCtrl_Q903", "weevilCtrl_H898")


lapply(focus_terms, function(x) {
  write.table(
    merge(subset(sigDE, sigDE$focus_term == x), annotsCollapsed, all.x = TRUE),
    paste0("../results/", x, ".p", pCutoff, "_lfc", lfc, ".annot.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE)
})


#####


# function to contigs that are common between two given focus_terms
findCommonContigs <- function(df, pCutoff, lfc, x, y) {
  subset(df, eval(parse(
    text = paste0("df$adj.P.Val.", x, " <= ", pCutoff, " & ",
                  "df$adj.P.Val.", y, " <= ", pCutoff, " & ",
                  "abs(df$logFC.", x, ") >= ", lfc, " & ",
                  "abs(df$logFC.", y, ") >= ", lfc)))) 
}


# Common contigs between feeding in Q903 and combined in Q903
weevilIndCommon <- findCommonContigs(sidf_wide, pCutoff, lfc, 
                                     "weevilInd_H898", "weevilInd_Q903")
str(weevilIndCommon)
# 2 obs. of  22 variables:

write.table(
  merge(weevilIndCommon, annotsCollapsed, all.x = TRUE),
  paste0("../results/weevilIndQ903_vs_weevilIndH898.p", pCutoff, "_lfc", lfc, ".annot.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)


# Common contigs between ctrl_vs_ctrl and cmpd in Q903
weevilInd_Q903_vs_constDiff <- 
  findCommonContigs(sidf_wide, pCutoff, lfc, "weevilInd_Q903", "constDiff")
str(weevilInd_Q903_vs_constDiff)
# 1010 obs. of  22 variables

write.table(
  merge(weevilInd_Q903_vs_constDiff, annotsCollapsed, all.x = TRUE),
  paste0("../results/weevilInd_Q903_vs_constDiff.p", pCutoff, "_lfc", lfc, ".annot.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)
