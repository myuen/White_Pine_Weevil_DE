library(testthat) # facilitate tests that will catch changes on re-analysis
library(plyr)     ## **ply() and revalue()

# Marshal output files from Sailfish
jFiles <- list.files("Sailfish-results/", full.names = TRUE)
test_that("Exactly 24 Sailfish output files are found",
          expect_equal(24, length(jFiles)))

## Extract, e.g. 898C1 from Sailfish-results/898C1_quant.sf
tmp <- gsub("([0-9]+[CGW][1234])_quant.sf", "\\1", basename(jFiles))
## Prepend "H" to 898 and "Q" to 903
tmp <- gsub("898", "H898", tmp)
tmp <- gsub("903", "Q903", tmp)
## Use as names for good side effects later
names(jFiles) <- tmp
  
## Read one file to learn how many rows we expect and to grab rownames
tmp <- read.table(jFiles[1], row.names = 1,
                  ## specifying colClasses speeds this up 2x
                  colClasses = rep(c("character", "numeric"), c(1, 6)))
(n <- nrow(tmp)) # 492317
jRowNames <- rownames(tmp)
test_that("First Sailfish output files has 492317 rows",
          expect_equal(492317, n))

## Read in all Sailfish data
system.time(
rawSailfishCounts <-
  aaply(jFiles, 1, function(x) {
    jDat <-
      read.table(x, row.names = 1, nrows = n * 1.1,
                 ## specifying colClasses speeds this up 2x
                 colClasses = rep(c("character", "numeric"), c(1, 6)))
    test_that("Sailfish output files have expected number of rows",
              expect_equal(n, nrow(jDat)))
    return(jDat$V7)
  })
) ## ~90 seconds for JB

## NOTE: rawSailfishCounts is transposed relative to what we expect / want at
## this point! Will put up with this here and resolve upon export / import.
colnames(rawSailfishCounts) <- jRowNames
str(rawSailfishCounts)
# num [1:24, 1:492317] 0 0 1.87 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ X1: chr [1:24] "H898C1" "H898C2" "H898C3" "H898C4" ...
# ..$   : chr [1:492317] "WPW_Inoculation_Trinity_C500_comp100002_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq2" "WPW_Inoculation_Trinity_C500_comp100011_c0_seq1" ...

# We discovered more ribosomal RNA post-assembly using BLAST.  The following 
# line removes putative ribosomal RNA from the table.
rRNA <- scan("putativeRibosomalRNA.id", what = "")
str(rRNA) # chr [1:389] "WPW_Inoculation_Trinity_C500_comp27782_c0_seq1" ...
summary(colnames(rawSailfishCounts) %in% rRNA)
#    Mode   FALSE    TRUE    NA's 
# logical  491928     389       0 
rawSailfishCounts <- rawSailfishCounts[, !(colnames(rawSailfishCounts) %in% rRNA)]
str(rawSailfishCounts) # num [1:24, 1:491928]

(n <- ncol(rawSailfishCounts)) # 491928
test_that("Sailfish output has 491928 rows after rRNA filtering",
          expect_equal(491928, n))

## enact the row / column transposition now
write.table(t(rawSailfishCounts), "consolidated-Sailfish-results.txt",
            sep = "\t", quote = FALSE)

## from JB checking she got same data as original
# library(tools) # md5sum()
# all.equal(md5sum("bak/consolidated-Sailfish-results-BAK.txt"),
#           md5sum("consolidated-Sailfish-results.txt"))
# [1] "Names: 1 string mismatch"
# presumably this is just the file name mismatch
# diff'ing in shell turns up no differences

## write design matrix

desMat <- data.frame(sample = rownames(rawSailfishCounts),
                     gTypeCode = factor(substr(rownames(rawSailfishCounts), 1, 4),
                                        levels = c("Q903", "H898")),
                     # H898 = Resistance Genotype; Q903 = Susceptible Genotype
                     gType = factor(NA),
                     txCode = substr(rownames(rawSailfishCounts), 5, 5),
                     # C = Control, G = Gallery, W = Wounding
                     tx = factor(NA),
                     bioRep = as.numeric(substr(rownames(rawSailfishCounts), 6, 6)))
desMat$gType <-
  revalue(desMat$gTypeCode, c('H898' = "res", 'Q903' = "susc"))
desMat$tx <-
  revalue(desMat$txCode, c('C' = "Control", 'G' = "Gallery", 'W' = "Wound"))
desMat$grp <- with(desMat, interaction(gTypeCode, tx))
desMat
str(desMat)

write.table(desMat, "White_Pine_Weevil_design_matrix.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
