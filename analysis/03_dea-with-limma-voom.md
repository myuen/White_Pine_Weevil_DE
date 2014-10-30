


Source purpose-built functions to load and validate the data and the 
experimental design. Then call them.


```r
source("helper01_load-counts.r")
source("helper02_load-exp-des.r")
```

Load counts from Sailfish


```r
x <- load_counts() # takes a few moments
str(x, list.len = 8) # 'data.frame':  65600 obs. of  24 variables:
```

```
## 'data.frame':	65600 obs. of  24 variables:
##  $ H898C1: num  85.7 136.4 639.5 57.1 411.4 ...
##  $ H898C2: num  70.8 98.7 550.2 49.9 227.2 ...
##  $ H898C3: num  132.9 144 472.2 51.5 398.2 ...
##  $ H898C4: num  122.6 134.5 488.4 50.9 268.1 ...
##  $ H898G1: num  99.9 121.4 422.3 60.1 423.1 ...
##  $ H898G2: num  99.9 167.9 420 49.8 339.4 ...
##  $ H898G3: num  68.1 123 411.6 36.8 285 ...
##  $ H898G4: num  98.3 62.2 291.4 46.6 208.6 ...
##   [list output truncated]
```

Load experimental design


```r
expDes <- load_expDes()
expDes$grp <- # not really sure that we need this?
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

Load counts into DGEList object from edgeR package.


```r
y <- DGEList(counts = x, group = expDes$grp)
```

TMM Normalization by Depth


```r
y <- calcNormFactors(y)
```

make model matrix


```r
modMat <- model.matrix(~ gType * tx, expDes)

# I never needed to fit the model with this parametrization. Instead I got the
# effects of interests from the main parametrization and specific contrasts
# specifed below.
#modMat <- model.matrix(~ tx/gType - 1, expDes)
```

It's hard to believe, but the default names for columns associated with 
interaction terms will create fatal errors in makeContrasts below; prevent 
that, and a warning about the intercept, by modifying these column names 
here; see 90_limma-model-term-name-fiasco.r and .md for more


```r
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

voom transformation


```r
v <- voom(y, modMat, plot = TRUE) # take a couple moments
```

![plot of chunk voom-plot](figure/03-dea-with-limma-voom-voom-plot-1.png) 

Linear modelling and forming contrasts


```r
fit <- lmFit(v, modMat)

cont_matrix <-
  makeContrasts(
    # A - Within genotype comparisons
    # A1 - Wound vs. Control
    wound_Q903 = txWound,
    wound_H898 = txWound + gTypeH898res_txWound,
    # A2 - Gallery vs. Wound
    feed_Q903 = txGallery - txWound,
    feed_H898 = txGallery - txWound +
      gTypeH898res_txGallery - gTypeH898res_txWound,
    # A3 - Gallery vs. Control (combined effect of wounding and feeding)
    combined_effect_Q903 = txGallery,
    combined_effect_H898 = txGallery + gTypeH898res_txGallery,
    # B - Between genotype comparisons
    # B1 - control vs control
    ctrl_vs_ctrl = gTypeH898res,
    # B2 - wound vs wound
    wound_vs_wound = gTypeH898res + gTypeH898res_txWound,
    # B3 - gallery vs gallery
    gallery_vs_gallery = gTypeH898res + gTypeH898res_txGallery,
    # C - Comparison of comparison
    # C1 -  Mwounding vs. wounding
    wounding_diff = gTypeH898res_txWound,
    # C2 - feeding vs. feeding i.e. (H898G - H898W) - (Q903G - Q903W)
    feeding_diff = gTypeH898res_txGallery - gTypeH898res_txWound,
    # C3 - combined effect vs. combined effect i.e. (H898G - H898C) - (Q903G - Q903C)
    combined_diff = gTypeH898res_txGallery,
    # D - Special case 
    # D1 - Induced vs Constitutive (Q903G - Q903W) - (H898C - Q903C)
    induced_vs_const = txGallery - txWound - gTypeH898res,
    levels = modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

summary(decideTests(fit3, p.value = 0.01, lfc = 2))
```

```
##    wound_Q903 wound_H898 feed_Q903 feed_H898 combined_effect_Q903
## -1          1          0      4632         4                 3896
## 0       65599      65600     59054     65596                58117
## 1           0          0      1914         0                 3587
##    combined_effect_H898 ctrl_vs_ctrl wound_vs_wound gallery_vs_gallery
## -1                    2         3679           3769               4390
## 0                 65597        57989          57971              55130
## 1                     1         3932           3860               6080
##    wounding_diff feeding_diff combined_diff induced_vs_const
## -1             0           99           175             5804
## 0          65600        64872         65151            56611
## 1              0          629           274             3185
```

```r
#    wound_Q903 wound_H898 feed_Q903 feed_H898 combined_effect_Q903
# -1          1          0      4632         4                 3896
# 0       65599      65600     59054     65596                58117
# 1           0          0      1914         0                 3587
#    combined_effect_H898 ctrl_vs_ctrl wound_vs_wound gallery_vs_gallery
# -1                    2         3679           3769               4390
# 0                 65597        57989          57971              55130
# 1                     1         3932           3860               6080
#    wounding_diff feeding_diff combined_diff induced_vs_const
# -1             0           99           175             5804
# 0          65600        64872         65151            56611
# 1              0          629           274             3185

focus_terms <- colnames(cont_matrix)

names(focus_terms) <- focus_terms
```

to get individual t statistics, etc., must use topTable() on each coef
separately


```r
statInf_focus_terms <-
  adply(focus_terms, 1,
        function(x) name_rows(topTable(fit3, coef = x,
                                       number = Inf, sort.by = "none")))

statInf_focus_terms <-
  rename(statInf_focus_terms, c(X1 = "focus_term", .rownames = "contig"))

statInf_focus_terms <-
  statInf_focus_terms[ ,c("contig",
                          setdiff(names(statInf_focus_terms),"contig"))]
str(statInf_focus_terms)
```

```
## 'data.frame':	852800 obs. of  8 variables:
##  $ contig    : chr  "WPW_Inoculation_Trinity_C500_comp100847_c0_seq1" "WPW_Inoculation_Trinity_C500_comp100915_c0_seq1" "WPW_Inoculation_Trinity_C500_comp101147_c0_seq1" "WPW_Inoculation_Trinity_C500_comp101332_c0_seq1" ...
##  $ focus_term: Factor w/ 13 levels "wound_Q903","wound_H898",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ logFC     : num  0.2973 -0.0506 -0.1176 0.0413 0.6248 ...
##  $ AveExpr   : num  0.424 0.395 2.714 -0.229 1.671 ...
##  $ t         : num  0.961 -0.169 -0.722 0.132 1.932 ...
##  $ P.Value   : num  0.3475 0.8672 0.4783 0.8963 0.0671 ...
##  $ adj.P.Val : num  0.959 0.997 0.997 0.997 0.713 ...
##  $ B         : num  -5.44 -5.87 -6.16 -5.8 -4.4 ...
```

```r
test_that("stat inf on the focus terms has 852,800 rows",
          expect_equal(65600 * 13, nrow(statInf_focus_terms)))

(t_medians <- aggregate(t ~ focus_term, statInf_focus_terms, median))
```

```
##              focus_term            t
## 1            wound_Q903  0.093017448
## 2            wound_H898  0.012387609
## 3             feed_Q903 -0.164613806
## 4             feed_H898 -0.093503376
## 5  combined_effect_Q903  0.021539496
## 6  combined_effect_H898 -0.089969465
## 7          ctrl_vs_ctrl  0.162680935
## 8        wound_vs_wound  0.011875199
## 9    gallery_vs_gallery  0.118430987
## 10        wounding_diff -0.067642621
## 11         feeding_diff  0.052106222
## 12        combined_diff  0.007929199
## 13     induced_vs_const -0.321207705
```

```r
all.equal(t_medians$t, c(0.093017448, 0.012387609, -0.164613806, -0.093503376, 
            0.021539496, -0.089969465, 0.162680935, 0.011875199, 
            0.118430987, -0.067642621, 0.052106222, 0.007929199, -0.321207705))
```

```
## [1] TRUE
```

```r
write.table(statInf_focus_terms,
            "../results/limma-results-focus-terms.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)
```

