library(edgeR)
library(ggplot2)
library(plyr)
library(reshape2)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# Load counts from Sailfish
rawSailfishCounts <- read.delim("consolidated-Sailfish-results.txt")
str(rawSailfishCounts) # 'data.frame':  491928 obs. of  24 variables:
test_that("Sailfish data has 491928 rows upon import",
          expect_equal(491928, nrow(rawSailfishCounts)))
test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(rawSailfishCounts)))

# Load design matrix
desMat <- read.delim("White_Pine_Weevil_design_matrix.tsv", stringsAsFactors = FALSE)
desMat <-
  mutate(desMat,
         gTypeCode = factor(gTypeCode, levels = c("Q903", "H898")),
         gType = factor(gType, levels = c("susc", "res")),
         txCode = factor(txCode),
         tx = factor(tx))
desMat$grp <-
  with(desMat, factor(grp,
                      levels = paste(levels(gTypeCode),
                                     rep(levels(tx), each = 2), sep = ".")))
str(desMat) # 'data.frame':  24 obs. of  7 variables:
test_that("design matrix has 24 rows upon import", expect_equal(24, nrow(desMat)))

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = rawSailfishCounts, group = desMat$grp)
lib_size <- data.frame(raw = y$samples$lib.size)

# exploring the phenomenon of low expression
non_zero_freq <- as.data.frame(with(y, table(rowSums(cpm(y) > 1))))
names(non_zero_freq) <- c("num.nonzero", "freq")
non_zero_freq$num.nonzero <-
  with(non_zero_freq,
       as.numeric(levels(num.nonzero)[num.nonzero]))
p <- ggplot(subset(non_zero_freq, num.nonzero > 1),
            aes(x = as.factor(num.nonzero), y = freq)) +
  geom_bar(stat = "identity")
p + coord_flip() + xlab("frequency or number of samples")

# Filtering low expression genes
# We are setting an arbitary threshold and only keeping contigs with at least
# 1 count-per-million (cpm) in at least half of the biological replicates
# in 1 timepoint (i.e. 2 samples)
y <- y[(rowSums(cpm(y) > 1) >= 2), ]
test_that("After low expression filter, we have 65609 rows",
          expect_equal(65609, nrow(y)))
# 65609 (down from 491928) ~= we have about 13% of original rows

# Library depth is now changed with the filtering of the low count 
# contigs so we need to reset the libray depth.
(lib_size$filtered <- colSums(y$counts))

# TMM Normalization by Depth
y <- calcNormFactors(y)

# make model matrix
#modMat <- model.matrix(~ gType/tx - 1, desMat)
modMat <- model.matrix(~ tx/gType - 1, desMat)

# voom transformation
v <- voom(y, modMat, plot = TRUE)

# Linear modelling
fit <- lmFit(v, modMat)
fit <- eBayes(fit)
tt <- topTable(fit, coef = grep(":", colnames(modMat)))

data.wide <- cbind(desMat, t(y$counts[rownames(tt)[1:4], ]))
data.tall <- melt(data.wide,
                  id.vars = names(desMat),
                  variable.name = 'contig', value.name = 'Expression')
str(data.tall)

#x <- data.frame(desMat, gExp = y$counts[rownames(tt)[1], ])
p <- ggplot(data.tall, aes(x = gType, y = Expression))
p + geom_point() + facet_grid(tx ~ contig)

# Linear design matrix, simplest approach
design <- model.matrix(~ 0 + grp, desMat)
colnames(design) <- levels(desMat$grp)

# Create contrast matrix
cont.matrix <- makeContrasts(
  Control.H898_vs_Q903 = H898.Control - Q903.Control,
  Induced.H898_vs_Q903 = 
    (H898.Gallery - (H898.Wound - H898.Control)) - 
    (Q903.Gallery - (Q903.Wound - Q903.Control)),
  levels = design)

# Voom transformation
v <- voom(y, modMat, plot = TRUE)



# MDS analysis
plotMDS(v, top=Inf)

# PCA analysis
pca <- as.data.frame(prcomp(t(v$E))$x)
ggplot(pca, aes(PC1, PC2, color = Group)) + geom_point(size=3) + 
  scale_color_manual(name = "", 
                     values=c("H898.Control" = "#4d004b", 
                              "H898.Gallery" = "#88419d", 
                              "H898.Wound" = "#8c96c6", 
                              "Q903.Control" = "#7F0000", 
                              "Q903.Gallery" =  "#d7301f", 
                              "Q903.Wound" = "#fc8d59")) + 
  labs(title="Voom + Limma Principal Component Analysis") + theme_bw()

# Linear modelling
fit <- lmFit(v, design)

fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2)

linear.results <- decideTests(fit2, method="global", 
                              adjust.method="fdr", p.value=0.01)

summary(linear.results)

###

vennDiagram(linear.results, include="both")
vennDiagram(linear.results, include="up")
vennDiagram(linear.results, include="down")

###
