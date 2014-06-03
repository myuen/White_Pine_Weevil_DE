



This report was automatically generated with the R package **knitr**
(version 1.5.33).


```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(ggplot2)
library(plyr)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# Load counts from Sailfish
rawSailfishCounts <- read.delim("../data/consolidated-Sailfish-results.txt")
str(rawSailfishCounts) # 'data.frame':  491928 obs. of  24 variables:
```

```
## 'data.frame':	491928 obs. of  24 variables:
##  $ H898C1: num  0 0 0 0 3.68 ...
##  $ H898C2: num  0 0 0 0 13.8 ...
##  $ H898C3: num  1.87 0 0 0 2.19 ...
##  $ H898C4: num  0 0 0 0 3.83 ...
##  $ H898G1: num  0 0 0 0 0 ...
##  $ H898G2: num  1.25 0 0 0 0 ...
##  $ H898G3: num  0 0 0 0 0 ...
##  $ H898G4: num  0 0 0 0 0 0 0 0 0 0 ...
##  $ H898W1: num  0 0 0 0 2.79 ...
##  $ H898W2: num  0 0 0 0 0 ...
##  $ H898W3: num  0 0 0 0 0 ...
##  $ H898W4: num  0 0 0 0 0 ...
##  $ Q903C1: num  5.18 15.56 4.56 0 0 ...
##  $ Q903C2: num  0 0 7.06 0 0 ...
##  $ Q903C3: num  1.69 0 0 0 5.99 ...
##  $ Q903C4: num  0 0.956 13.083 0 0 ...
##  $ Q903G1: num  0 0 0 0 0 0 0 0 0 0 ...
##  $ Q903G2: num  0 0 0 10.45 5.21 ...
##  $ Q903G3: num  0 0 0 0 0 ...
##  $ Q903G4: num  0.714 0 0 0 0 ...
##  $ Q903W1: num  0 0 0 0 0 ...
##  $ Q903W2: num  1.56 0 0 0 0 ...
##  $ Q903W3: num  0 0 0 0 0.948 ...
##  $ Q903W4: num  0 0 0 0 0 ...
```

```r
test_that("Sailfish data has 491928 rows upon import",
          expect_equal(491928, nrow(rawSailfishCounts)))
test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(rawSailfishCounts)))

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = rawSailfishCounts)
lib_size <- data.frame(raw = y$samples$lib.size)

# exploring the phenomenon of low expression:
# convert to counts per million
# is cpm > 1? <-- our effective definition of "present"
# sum across all 24 samples
# tabulate that frequency
present_sample_count <- as.data.frame(with(y, table(rowSums(cpm(y) > 1))))
names(present_sample_count) <- c("num.present", "freq")
present_sample_count$num.present <-
  with(present_sample_count,
       as.numeric(levels(num.present)[num.present]))
p <- ggplot(subset(present_sample_count, num.present > 0),
            aes(x = as.factor(num.present), y = freq)) +
  geom_bar(stat = "identity") + coord_flip() +
  xlab("frequency of: number of samples contig is present in")
p <- p + annotate("text", y = Inf, x = 12, hjust = 1.1,
                  label = 'contig called "present" in a sample if cpm > 1') +
  annotate("text", y = Inf, x = 3, hjust = 1.1,
           label = paste(with(present_sample_count, freq[num.present == 0]),
                         'contigs called "present" in NO samples')) +
  annotate("text", y = Inf, x = 22.5, hjust = 1.1,
           label = paste(with(present_sample_count, freq[num.present == 24]),
                         'contigs called "present" in all samples'))
ggsave("figure/02_pre-dea-filtering-preDE-filtering.png", plot = p,
       height = 7, width = 7)

## specifying height and width prevents an empty Rplots.pdf file from being left
## behind; see
## http://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript

# Filtering low expression genes
# We are setting an arbitary threshold and only keeping contigs with more than
# 1 count-per-million (cpm) in at least 2 samples
y <- y[(rowSums(cpm(y) > 1) > 1), ]
test_that("After low expression filter, we have 65609 rows",
          expect_equal(65609, nrow(y)))
# 65609 (down from 491928) ~= we have about 13% of original rows

## write to file
write.table(y$counts, "../data/consolidated-filtered-Sailfish-results.txt",
            sep = "\t", quote = FALSE)
```

The R session information (including the OS info, R version and all
packages used):


```r
sessionInfo()
```

```
## R version 3.1.0 (2014-04-10)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] C
## 
## attached base packages:
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] testthat_0.8.1  plyr_1.8.1      ggplot2_0.9.3.1 edgeR_3.4.2    
## [5] limma_3.18.13  
## 
## loaded via a namespace (and not attached):
##  [1] MASS_7.3-33      Rcpp_0.11.1      colorspace_1.2-4 digest_0.6.4    
##  [5] evaluate_0.5.5   formatR_0.10     grid_3.1.0       gtable_0.1.2    
##  [9] knitr_1.5.33     labeling_0.2     munsell_0.4.2    proto_0.3-10    
## [13] reshape2_1.4     scales_0.2.4     stringr_0.6.2    tools_3.1.0
```

```r
Sys.time()
```

```
## [1] "2014-06-03 15:38:13 PDT"
```

