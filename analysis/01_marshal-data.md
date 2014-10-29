



This report was automatically generated with the R package **knitr**
(version 1.7).


```r
library(testthat) # facilitate tests that will catch changes on re-analysis
library(plyr)     ## **ply() and revalue()

# Marshal output files from Sailfish
jFiles <- list.files("../data/Sailfish-results/", full.names = TRUE)
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
```

```
## [1] 492317
```

```r
jRowNames <- rownames(tmp)
test_that("First Sailfish output file has 492317 rows",
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
```

```
##    user  system elapsed 
##  72.610   0.608  73.609
```

```r
## NOTE: rawSailfishCounts is transposed relative to what we expect / want at
## this point! Will put up with this here and resolve upon export / import.
colnames(rawSailfishCounts) <- jRowNames
str(rawSailfishCounts)
```

```
##  num [1:24, 1:492317] 0 0 1.87 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ X1: chr [1:24] "H898C1" "H898C2" "H898C3" "H898C4" ...
##   ..$   : chr [1:492317] "WPW_Inoculation_Trinity_C500_comp100002_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq2" "WPW_Inoculation_Trinity_C500_comp100011_c0_seq1" ...
```

```r
# num [1:24, 1:492317] 0 0 1.87 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ X1: chr [1:24] "H898C1" "H898C2" "H898C3" "H898C4" ...
# ..$   : chr [1:492317] "WPW_Inoculation_Trinity_C500_comp100002_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq2" "WPW_Inoculation_Trinity_C500_comp100011_c0_seq1" ...


# We discovered more ribosomal RNA post-assembly using BLAST.  The following 
# line removes putative ribosomal RNA from the table.
rRNA <- scan("../data/putativeRibosomalRNA.id", what = "")
str(rRNA) # chr [1:389] "WPW_Inoculation_Trinity_C500_comp27782_c0_seq1" ...
```

```
##  chr [1:389] "WPW_Inoculation_Trinity_C500_comp27782_c0_seq1" ...
```

```r
summary(colnames(rawSailfishCounts) %in% rRNA)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical  491928     389       0
```

```r
#    Mode   FALSE    TRUE    NA's 
# logical  491928     389       0 

microbial <- scan("../data/microbialContamination.id", what = "")
str(microbial) #  chr [1:181] "WPW_Inoculation_Trinity_C500_comp150393_c0_seq1" ...
```

```
##  chr [1:181] "WPW_Inoculation_Trinity_C500_comp150393_c0_seq1" ...
```

```r
summary(colnames(rawSailfishCounts) %in% microbial)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical  492136     181       0
```

```r
# Mode   FALSE    TRUE    NA's 
# logical  492136     181       0

human <- scan("../data/humanContamination.id", what = "")
str(human) # chr [1:224] "WPW_Inoculation_Trinity_C500_comp107462_c0_seq1" ...
```

```
##  chr [1:224] "WPW_Inoculation_Trinity_C500_comp107462_c0_seq1" ...
```

```r
summary(colnames(rawSailfishCounts) %in% human)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical  492093     224       0
```

```r
# Mode   FALSE    TRUE    NA's 
# logical  492093     224       0

weevil <- scan("../data/weevilContamination.id", what = "")
str(weevil) # chr [1:403] "WPW_Inoculation_Trinity_C500_comp100267_c0_seq1" ...
```

```
##  chr [1:403] "WPW_Inoculation_Trinity_C500_comp100267_c0_seq1" ...
```

```r
summary(colnames(rawSailfishCounts) %in% weevil)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical  491914     403       0
```

```r
# Mode   FALSE    TRUE    NA's 
# logical  491914     403       0

fungal <- scan("../data/fungalContamination.id", what = "")
str(fungal) # chr [1:8361] "WPW_Inoculation_Trinity_C500_comp100413_c0_seq1" ...
```

```
##  chr [1:8361] "WPW_Inoculation_Trinity_C500_comp100413_c0_seq1" ...
```

```r
summary(colnames(rawSailfishCounts) %in% fungal)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical  483956    8361       0
```

```r
#    Mode   FALSE    TRUE    NA's 
# logical  483956    8361       0

contaminants <- c(rRNA, human, microbial, fungal, weevil)
contaminants <- unique(sort(contaminants))
str(contaminants) # chr [1:9270] "WPW_Inoculation_Trinity_C500_comp100267_c0_seq1" ...
```

```
##  chr [1:9270] "WPW_Inoculation_Trinity_C500_comp100267_c0_seq1" ...
```

```r
summary(colnames(rawSailfishCounts) %in% contaminants)
```

```
##    Mode   FALSE    TRUE    NA's 
## logical  483047    9270       0
```

```r
# Mode   FALSE    TRUE    NA's 
# logical  483047    9270       0 

rawSailfishCounts <- rawSailfishCounts[, !(colnames(rawSailfishCounts) %in% contaminants)]
str(rawSailfishCounts) # num [1:24, 1:483047]
```

```
##  num [1:24, 1:483047] 0 0 1.87 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ X1: chr [1:24] "H898C1" "H898C2" "H898C3" "H898C4" ...
##   ..$   : chr [1:483047] "WPW_Inoculation_Trinity_C500_comp100002_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100009_c0_seq2" "WPW_Inoculation_Trinity_C500_comp100011_c0_seq1" ...
```

```r
(n <- ncol(rawSailfishCounts)) # 483047
```

```
## [1] 483047
```

```r
test_that("Sailfish output has 483047 rows after rRNA filtering",
          expect_equal(483047, n))


## enact the row / column transposition now
write.table(t(rawSailfishCounts), "../data/consolidated-Sailfish-results.txt",
            sep = "\t", quote = FALSE)


## from JB checking she got same data as original
# library(tools) # md5sum()
# all.equal(md5sum("bak/consolidated-Sailfish-results-BAK.txt"),
#           md5sum("consolidated-Sailfish-results.txt"))
# [1] "Names: 1 string mismatch"
# presumably this is just the file name mismatch
# diff'ing in shell turns up no differences

## write design matrix

expDes <- data.frame(sample = rownames(rawSailfishCounts),
                     gType = factor(substr(rownames(rawSailfishCounts), 1, 4),
                                        levels = c("Q903", "H898")),
                     txCode = factor(substr(rownames(rawSailfishCounts), 5, 5),
                                     levels = c('C', 'W', 'G')),
                     # C = Control, W = Wounding, G = Gallery 
                     tx = factor(NA),
                     bioRep = as.numeric(substr(rownames(rawSailfishCounts), 6, 6)))
expDes$gType <- # H898 = Resistance Genotype; Q903 = Susceptible Genotype
  revalue(expDes$gType, c('H898' = "H898res", 'Q903' = "Q903susc"))
expDes$tx <-
  revalue(expDes$txCode, c('C' = "Control", 'W' = "Wound", 'G' = "Gallery"))
expDes$grp <- with(expDes, interaction(gType, tx))
expDes
```

```
##    sample    gType txCode      tx bioRep              grp
## 1  H898C1  H898res      C Control      1  H898res.Control
## 2  H898C2  H898res      C Control      2  H898res.Control
## 3  H898C3  H898res      C Control      3  H898res.Control
## 4  H898C4  H898res      C Control      4  H898res.Control
## 5  H898G1  H898res      G Gallery      1  H898res.Gallery
## 6  H898G2  H898res      G Gallery      2  H898res.Gallery
## 7  H898G3  H898res      G Gallery      3  H898res.Gallery
## 8  H898G4  H898res      G Gallery      4  H898res.Gallery
## 9  H898W1  H898res      W   Wound      1    H898res.Wound
## 10 H898W2  H898res      W   Wound      2    H898res.Wound
## 11 H898W3  H898res      W   Wound      3    H898res.Wound
## 12 H898W4  H898res      W   Wound      4    H898res.Wound
## 13 Q903C1 Q903susc      C Control      1 Q903susc.Control
## 14 Q903C2 Q903susc      C Control      2 Q903susc.Control
## 15 Q903C3 Q903susc      C Control      3 Q903susc.Control
## 16 Q903C4 Q903susc      C Control      4 Q903susc.Control
## 17 Q903G1 Q903susc      G Gallery      1 Q903susc.Gallery
## 18 Q903G2 Q903susc      G Gallery      2 Q903susc.Gallery
## 19 Q903G3 Q903susc      G Gallery      3 Q903susc.Gallery
## 20 Q903G4 Q903susc      G Gallery      4 Q903susc.Gallery
## 21 Q903W1 Q903susc      W   Wound      1   Q903susc.Wound
## 22 Q903W2 Q903susc      W   Wound      2   Q903susc.Wound
## 23 Q903W3 Q903susc      W   Wound      3   Q903susc.Wound
## 24 Q903W4 Q903susc      W   Wound      4   Q903susc.Wound
```

```r
str(expDes)
```

```
## 'data.frame':	24 obs. of  6 variables:
##  $ sample: Factor w/ 24 levels "H898C1","H898C2",..: 1 2 3 4 5 6 7 8 9 10 ...
##  $ gType : Factor w/ 2 levels "Q903susc","H898res": 2 2 2 2 2 2 2 2 2 2 ...
##  $ txCode: Factor w/ 3 levels "C","W","G": 1 1 1 1 3 3 3 3 2 2 ...
##  $ tx    : Factor w/ 3 levels "Control","Wound",..: 1 1 1 1 3 3 3 3 2 2 ...
##  $ bioRep: num  1 2 3 4 1 2 3 4 1 2 ...
##  $ grp   : Factor w/ 6 levels "Q903susc.Control",..: 2 2 2 2 6 6 6 6 4 4 ...
```

```r
write.table(expDes, "../data/White_Pine_Weevil_exp_design.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

The R session information (including the OS info, R version and all
packages used):


```r
sessionInfo()
```

```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
##  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
##  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] plyr_1.8.1     testthat_0.9.1
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_1.0    knitr_1.7      Rcpp_0.11.3   
## [5] stringr_0.6.2  tools_3.1.1
```

```r
Sys.time()
```

```
## [1] "2014-10-29 12:44:25 PDT"
```

