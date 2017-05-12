library(plyr)

setwd("analysis/")
source("helper03_load-focus-statinf.r")

sidf <- load_focus_statInf()
str(sidf)


# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05



# Count all differentially expressed contigs for each focus_term
ddply(sidf, ~ focus_term, function(x) {
  up <- length(subset(x, x$adj.P.Val <= pCutoff & x$logFC >= lfcCutoff)[,"cds"]);
  down <- length(subset(x, x$adj.P.Val <= pCutoff & x$logFC <= -1 * lfcCutoff)[,"cds"]);
  sum <- up + down;
  return (c("up" = up, "down" = down, "sum" = sum));
})

#        focus_term   up down  sum
# 1       constDiff 2279 2242 4521
# 2  woundResp_Q903    0    1    1
# 3  woundResp_H898    0    0    0
# 4  weevilInd_Q903 2000 3763 5763
# 5  weevilInd_H898    0    0    0
# 6 weevilCtrl_Q903 3319 3925 7244
# 7 weevilCtrl_H898    2    1    3


# BLAST results output on tab-delimited format (i.e. run with -outfmt '6' 
# option on comamnd line BLAST).  We ran with 10 max target hits returned 
# per query sequences.
annots <- read.delim(
    "../results/WPW_Inoculation_Trinity_BBMerged_31May2016.transdecoder.cds.nr.diffExp.blastxNR.txt", 
    sep = "\t", header = TRUE)
dim(annots) # [1] 115033     15


# Function to concatenate multiple blast annontations for each query sequence
concatAnnot <- function(x) {
  annot <- x$salltitles
  annot <- gsub("\\<\\>", " ", annot, perl = TRUE)
  return (paste(annot, collapse = " ; "))
}

annotsCollapsed <-
  ddply(annots, ~ qseqid, concatAnnot)
colnames(annotsCollapsed) <- c("contig", "annot")
str(annotsCollapsed) # 'data.frame':	10301 obs. of  2 variables:

write.table(annotsCollapsed,
            "../results/WPW_Inoculation_Trinity_BBMerged_31May2016.transdecoder.cds.nr.diffExp.blastxNR.collapsed.txt",
  sep = "\t", quote = FALSE, row.names = FALSE)


# Find all DE contigs
sigDE <- 
  ddply(sidf, ~ focus_term, 
        function(x) {
          tmp_sigDE <- subset(x, (abs(logFC)) >= lfcCutoff & adj.P.Val <= pCutoff);
          tmp_sigDE <- tmp_sigDE[order(tmp_sigDE$logFC, decreasing = TRUE), ];
        })

focus_terms <- as.vector(levels(sidf$focus_term))

lapply(focus_terms, function(x) {
  write.table(
    merge(subset(sigDE, sigDE$focus_term == x), annotsCollapsed, 
          by.x = "cds", by.y = "cds", all.x = TRUE),
    paste0("../results/", x, ".p", pCutoff, "_lfc", lfcCutoff, ".annotated.31May.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE)
})


# function to contigs that are common between two given focus_terms
findCommonContigs <- function(df, pCutoff, lfc, x, y) {
  subset(df, eval(parse(
    text = paste0("df$adj.P.Val.", x, " <= ", pCutoff, " & ",
                  "df$adj.P.Val.", y, " <= ", pCutoff, " & ",
                  "abs(df$logFC.", x, ") >= ", lfcCutoff, " & ",
                  "abs(df$logFC.", y, ") >= ", lfcCutoff)))) 
}


# Common contigs between ctrl_vs_ctrl and cmpd in Q903
weevilInd_Q903_vs_constDiff <- 
  findCommonContigs(sidf_wide, pCutoff, lfc, "weevilInd_Q903", "constDiff")
str(weevilInd_Q903_vs_constDiff)
# 766 obs. of  15 variables


write.table(
  merge(weevilInd_Q903_vs_constDiff, annotsCollapsed, all.x = TRUE),
  paste0("../results/weevilInd_Q903_vs_constDiff.p", pCutoff, "_lfc", lfcCutoff, ".annot.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)
