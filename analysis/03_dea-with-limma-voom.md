



This report was automatically generated with the R package **knitr**
(version 1.5.33).


```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(plyr)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

#' Source purpose-built functions to load and validate the data, the
#' experimental design, and the statistical inference results for our focus
#' terms. Then call them.
source("helper01_load-counts.r")
source("helper02_load-exp-des.r")

# Load counts from Sailfish
x <- load_counts() # takes a few moments
str(x, list.len = 8) # 'data.frame':  65609 obs. of  24 variables:
```

```
## 'data.frame':	65609 obs. of  24 variables:
##  $ H898C1: num  9450.8 85.7 136.4 639.5 57.1 ...
##  $ H898C2: num  7777.8 70.8 98.7 550.2 49.9 ...
##  $ H898C3: num  7718.5 132.9 144 472.2 51.5 ...
##  $ H898C4: num  8283.8 122.6 134.5 488.4 50.9 ...
##  $ H898G1: num  8717.3 99.9 121.4 422.3 60.1 ...
##  $ H898G2: num  8209.8 99.9 167.9 420 49.8 ...
##  $ H898G3: num  6860.8 68.1 123 411.6 36.8 ...
##  $ H898G4: num  5362.4 98.3 62.2 291.4 46.6 ...
##   [list output truncated]
```

```r
# Load experimental design
expDes <- load_expDes()
expDes$grp <-
  with(expDes, factor(grp,
                      levels = paste(levels(gType),
                                     rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables:
```

```
## 'data.frame':	24 obs. of  6 variables:
##  $ sample: chr  "H898C1" "H898C2" "H898C3" "H898C4" ...
##  $ gType : Factor w/ 2 levels "Q903susc","H898res": 2 2 2 2 2 2 2 2 2 2 ...
##  $ txCode: Factor w/ 3 levels "C","W","G": 1 1 1 1 3 3 3 3 2 2 ...
##  $ tx    : Factor w/ 3 levels "Control","Wound",..: 1 1 1 1 3 3 3 3 2 2 ...
##  $ bioRep: int  1 2 3 4 1 2 3 4 1 2 ...
##  $ grp   : Factor w/ 6 levels "Q903susc.Control",..: 2 2 2 2 6 6 6 6 4 4 ...
```

```r
# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = x, group = expDes$grp)

# TMM Normalization by Depth
y <- calcNormFactors(y)

# make model matrix
modMat <- model.matrix(~ gType * tx, expDes)

# I never needed to fit the model with this parametrization. Instead I got the
# effects of interests from the main parametrization and specific contrasts
# specifed below.
#modMat <- model.matrix(~ tx/gType - 1, expDes)

# hard to believe, but the default names for columns associated with interaction
# terms will create fatal errors in makeContrasts below; prevent that, and a
# warning about the intercept, by modifying these column names here; see
# White_Pine_Weevil_Sailfish_limm-model-term-name-fiasco.R and .md for more
colnames(modMat)
```

```
## [1] "(Intercept)"            "gTypeH898res"          
## [3] "txWound"                "txGallery"             
## [5] "gTypeH898res:txWound"   "gTypeH898res:txGallery"
```

```r
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat)
```

```
## [1] "Intercept"              "gTypeH898res"          
## [3] "txWound"                "txGallery"             
## [5] "gTypeH898res_txWound"   "gTypeH898res_txGallery"
```

```r
# voom transformation
pdf("figure/03_dea-with-limma-voom_voom-plot.pdf")
v <- voom(y, modMat, plot = TRUE) # take a couple moments
dev.off()
```

```
## pdf 
##   2
```

```r
# Linear modelling
fit <- lmFit(v, modMat)
cont_matrix <-
  makeContrasts(Intercept, gTypeH898res, txWound, txGallery,
                gTypeH898res_txWound, gTypeH898res_txGallery,
                weevil = gTypeH898res_txGallery - gTypeH898res_txWound,
                wound = gTypeH898res + gTypeH898res_txWound,
                gallery = gTypeH898res + gTypeH898res_txGallery,
                levels = modMat)
fit2 <- contrasts.fit(fit, cont_matrix)
fit3 <- eBayes(fit2)

# get inferential summary for the six terms in the model

# to get individual t statistics, etc., must use topTable() on each coef
# separately
model_terms <- colnames(modMat)
names(model_terms) <- model_terms
statInf_model_terms <-
  adply(model_terms, 1,
        function(x) name_rows(topTable(fit3, coef = x,
                                       number = Inf,sort.by = "none")))
statInf_model_terms <-
  rename(statInf_model_terms, c(X1 = "model_term", .rownames = "contig"))
statInf_model_terms <-
  statInf_model_terms[ ,c("contig", setdiff(names(statInf_model_terms),"contig"))]
head(statInf_model_terms)
```

```
##                                            contig model_term     logFC
## 1 WPW_Inoculation_Trinity_C500_comp100267_c0_seq1  Intercept  6.513861
## 2 WPW_Inoculation_Trinity_C500_comp100847_c0_seq1  Intercept  0.172678
## 3 WPW_Inoculation_Trinity_C500_comp100915_c0_seq1  Intercept  0.353133
## 4 WPW_Inoculation_Trinity_C500_comp101147_c0_seq1  Intercept  2.853620
## 5 WPW_Inoculation_Trinity_C500_comp101332_c0_seq1  Intercept -0.004685
## 6 WPW_Inoculation_Trinity_C500_comp103014_c0_seq2  Intercept  1.037792
##   AveExpr         t   P.Value adj.P.Val       B
## 1  6.5829 120.47702 4.301e-31 1.411e-27 57.8451
## 2  0.4102   0.74627 4.638e-01 5.138e-01 -6.5607
## 3  0.3805   1.64358 1.153e-01 1.493e-01 -5.5743
## 4  2.7002  24.73168 6.475e-17 5.337e-16 28.7967
## 5 -0.2426  -0.02072 9.837e-01 9.863e-01 -6.8046
## 6  1.6571   4.16973 4.398e-04 9.124e-04 -0.4998
```

```r
str(statInf_model_terms)
```

```
## 'data.frame':	393654 obs. of  8 variables:
##  $ contig    : chr  "WPW_Inoculation_Trinity_C500_comp100267_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100847_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100915_c0_seq1" "WPW_Inoculation_Trinity_C500_comp101147_c0_seq1" ...
##  $ model_term: Factor w/ 6 levels "Intercept","gTypeH898res",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ logFC     : num  6.51386 0.17268 0.35313 2.85362 -0.00468 ...
##  $ AveExpr   : num  6.583 0.41 0.381 2.7 -0.243 ...
##  $ t         : num  120.477 0.7463 1.6436 24.7317 -0.0207 ...
##  $ P.Value   : num  4.30e-31 4.64e-01 1.15e-01 6.47e-17 9.84e-01 ...
##  $ adj.P.Val : num  1.41e-27 5.14e-01 1.49e-01 5.34e-16 9.86e-01 ...
##  $ B         : num  57.85 -6.56 -5.57 28.8 -6.8 ...
```

```r
## informal tests so we know if things change, differ for Jenny vs Mack, etc.
test_that("stat inf on the model terms has 393654 rows",
          expect_equal(393654, nrow(statInf_model_terms)))
(t_medians <- aggregate(t ~ model_term, statInf_model_terms, median))
```

```
##               model_term         t
## 1              Intercept  1.118479
## 2           gTypeH898res  0.157362
## 3                txWound  0.092882
## 4              txGallery  0.026699
## 5   gTypeH898res_txWound -0.067349
## 6 gTypeH898res_txGallery  0.005008
```

```r
all.equal(t_medians$t, c(1.118479141, 0.157361758, 0.092882336,
                         0.026698663, -0.067348886, 0.005007678))
```

```
## [1] TRUE
```

```r
write.table(statInf_model_terms,
            "../results/limma-results-model-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

# get inferential summary for the effects we are most interested in
focus_patterns <-
  c("weevil", "gTypeH898res", "^gTypeH898res$", "wound", "gallery")
focus_terms <-
  c("weevil", "gTypeH898res_all", "gTypeH898res_in_control",
    "gTypeH898res_in_wound", "gTypeH898res_in_gallery")
names(focus_patterns) <- focus_terms
statInf_focus_terms <-
  alply(focus_patterns, 1,
        function(x) name_rows(topTable(fit3,
                                       coef = grep(x, colnames(coef(fit3))),
                                       number = Inf, sort.by = "none")))
names(statInf_focus_terms) <- focus_terms

#' I must massage these results before I can rbind them
#' Must address fact that they don't have exactly the same variables
var_names <- unique(unlist(llply(statInf_focus_terms, names)))
ldply(statInf_focus_terms, function(x) {
  y <- var_names %in% names(x)
  names(y) <- var_names
  y
})
```

```
##                        X1 logFC AveExpr     t P.Value adj.P.Val     B
## 1                  weevil  TRUE    TRUE  TRUE    TRUE      TRUE  TRUE
## 2        gTypeH898res_all FALSE    TRUE FALSE    TRUE      TRUE FALSE
## 3 gTypeH898res_in_control  TRUE    TRUE  TRUE    TRUE      TRUE  TRUE
## 4   gTypeH898res_in_wound  TRUE    TRUE  TRUE    TRUE      TRUE  TRUE
## 5 gTypeH898res_in_gallery  TRUE    TRUE  TRUE    TRUE      TRUE  TRUE
##   .rownames gTypeH898res gTypeH898res_txWound gTypeH898res_txGallery     F
## 1      TRUE        FALSE                FALSE                  FALSE FALSE
## 2      TRUE         TRUE                 TRUE                   TRUE  TRUE
## 3      TRUE        FALSE                FALSE                  FALSE FALSE
## 4      TRUE        FALSE                FALSE                  FALSE FALSE
## 5      TRUE        FALSE                FALSE                  FALSE FALSE
```

```r
#' get rid of the three separate estimates in the case where we are testing for
#' equality with zero for three terms at once
statInf_focus_terms$gTypeH898res_all <-
  subset(statInf_focus_terms$gTypeH898res_all,
         select = -c(gTypeH898res, gTypeH898res_txWound, gTypeH898res_txGallery))

#' rbind them
statInf_focus_terms <- ldply(statInf_focus_terms, function(x) x)
statInf_focus_terms <-
  rename(statInf_focus_terms, c(X1 = "focus_term", .rownames = "contig"))
summary(statInf_focus_terms)
```

```
##                    focus_term        logFC          AveExpr      
##  weevil                 :65609   Min.   :-22     Min.   :-6.525  
##  gTypeH898res_all       :65609   1st Qu.: -1     1st Qu.:-2.139  
##  gTypeH898res_in_control:65609   Median :  0     Median : 0.040  
##  gTypeH898res_in_wound  :65609   Mean   :  0     Mean   : 0.316  
##  gTypeH898res_in_gallery:65609   3rd Qu.:  1     3rd Qu.: 2.707  
##                                  Max.   : 15     Max.   :12.904  
##                                  NA's   :65609                   
##        t            P.Value         adj.P.Val            B        
##  Min.   :-27     Min.   :0.0000   Min.   :0.0000   Min.   :-9     
##  1st Qu.: -1     1st Qu.:0.0157   1st Qu.:0.0634   1st Qu.:-6     
##  Median :  0     Median :0.1852   Median :0.3857   Median :-5     
##  Mean   :  0     Mean   :0.3001   Mean   :0.4196   Mean   :-4     
##  3rd Qu.:  1     3rd Qu.:0.5388   3rd Qu.:0.7428   3rd Qu.:-4     
##  Max.   : 31     Max.   :1.0000   Max.   :1.0000   Max.   :27     
##  NA's   :65609                                     NA's   :65609  
##     contig                F         
##  Length:328045      Min.   :  0     
##  Class :character   1st Qu.:  1     
##  Mode  :character   Median :  3     
##                     Mean   : 11     
##                     3rd Qu.:  7     
##                     Max.   :913     
##                     NA's   :262436
```

```r
#' rearrange the variables
vars_in_order <- c("contig", "focus_term", "logFC", "AveExpr",
                   "t", "F", "P.Value", "adj.P.Val", "B")
statInf_focus_terms <- statInf_focus_terms[vars_in_order]
head(statInf_focus_terms)
```

```
##                                            contig focus_term    logFC
## 1 WPW_Inoculation_Trinity_C500_comp100267_c0_seq1     weevil  0.27711
## 2 WPW_Inoculation_Trinity_C500_comp100847_c0_seq1     weevil -0.09150
## 3 WPW_Inoculation_Trinity_C500_comp100915_c0_seq1     weevil  0.97003
## 4 WPW_Inoculation_Trinity_C500_comp101147_c0_seq1     weevil -0.02553
## 5 WPW_Inoculation_Trinity_C500_comp101332_c0_seq1     weevil  0.05371
## 6 WPW_Inoculation_Trinity_C500_comp103014_c0_seq2     weevil  0.85351
##   AveExpr       t  F P.Value adj.P.Val      B
## 1  6.5829  2.5788 NA 0.01758    0.1449 -4.240
## 2  0.4102 -0.2072 NA 0.83788    0.9433 -6.053
## 3  0.3805  2.2200 NA 0.03766    0.2200 -3.877
## 4  2.7002 -0.1069 NA 0.91586    0.9721 -6.598
## 5 -0.2426  0.1095 NA 0.91382    0.9714 -5.910
## 6  1.6571  1.9201 NA 0.06865    0.3031 -4.622
```

```r
str(statInf_focus_terms)
```

```
## 'data.frame':	328045 obs. of  9 variables:
##  $ contig    : chr  "WPW_Inoculation_Trinity_C500_comp100267_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100847_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100915_c0_seq1" "WPW_Inoculation_Trinity_C500_comp101147_c0_seq1" ...
##  $ focus_term: Factor w/ 5 levels "weevil","gTypeH898res_all",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ logFC     : num  0.2771 -0.0915 0.97 -0.0255 0.0537 ...
##  $ AveExpr   : num  6.583 0.41 0.381 2.7 -0.243 ...
##  $ t         : num  2.579 -0.207 2.22 -0.107 0.11 ...
##  $ F         : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ P.Value   : num  0.0176 0.8379 0.0377 0.9159 0.9138 ...
##  $ adj.P.Val : num  0.145 0.943 0.22 0.972 0.971 ...
##  $ B         : num  -4.24 -6.05 -3.88 -6.6 -5.91 ...
```

```r
## informal tests so we know if things change, differ for Jenny vs Mack, etc.
test_that("stat inf on the focus terms has 328045 rows",
          expect_equal(328045, nrow(statInf_focus_terms)))
(t_medians <- aggregate(t ~ focus_term, statInf_focus_terms, median))
```

```
##                focus_term       t
## 1                  weevil 0.04818
## 2 gTypeH898res_in_control 0.15736
## 3   gTypeH898res_in_wound 0.01199
## 4 gTypeH898res_in_gallery 0.10724
```

```r
all.equal(t_medians$t, c(0.04817940, 0.157361758, 0.01198738, 0.10723831))
```

```
## [1] "Mean relative difference: 2.508e-08"
```

```r
# "Mean relative difference: 2.508257e-08" <-- that's OK!

write.table(statInf_focus_terms,
            "../results/limma-results-focus-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)
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
## [1] testthat_0.8.1 plyr_1.8.1     edgeR_3.4.2    limma_3.18.13 
## 
## loaded via a namespace (and not attached):
## [1] Rcpp_0.11.1    digest_0.6.4   evaluate_0.5.5 formatR_0.10  
## [5] knitr_1.5.33   stringr_0.6.2  tools_3.1.0
```

```r
Sys.time()
```

```
## [1] "2014-06-03 14:30:49 PDT"
```

