### Differential Expression Analysis on Sitka Spruce Weevil Experiment with limma + voom

# Load counts from RSEM
rawSailfishCounts <- read.table("consolidated-Sailfish-results.txt", header = TRUE)

# Load edgeR library
library(edgeR)

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts=rawSailfishCounts)

# Filtering for low expression genes
# We are setting an arbitary threshold and only keeping contigs with at least
# 1 count-per-million (cpm) in at least half of the biological replicates
# in 1 timepoint (i.e. 2 samples)
y <- y[(rowSums(cpm(y) > 1) >= 2), ]

# Library depth is now changed with the filtering of the low count 
# contigs so we need to reset the libray depth.
y$samples$lib.size <- colSums(y$counts)

# TMM Normalization by Depth
y <- calcNormFactors(y)

# Create design matrix
targets <- cbind(colnames(rawSailfishCounts))
targets <- cbind(targets, c(rep("H898", 12), rep("Q903", 12)))
targets <- cbind(targets, rep(c(rep("Control", 4), 
                                rep("Gallery", 4), 
                                rep("Wound", 4)), 2))
colnames(targets) <- c("Sample", "Genotype", "Treatment")
targets <- as.data.frame(targets)

Genotype <- factor(targets$Genotype, levels=c("H898", "Q903"))
Treatment <- factor(targets$Treatment, levels=c("Control", "Gallery", "Wound"))

Group <- factor(paste(targets$Genotype, targets$Treatment, sep="_"))
targets <- cbind(targets, Group=Group)

# Linear design matrix, simplest approach
linear.design <- model.matrix(~ 0 + Group)
colnames(linear.design) <- levels(targets$Group)

# Create contrast matrix
linear.cont.matrix <- makeContrasts(
  Gallery_vs_Control_H898 = H898_Gallery - H898_Control,
  Wound_vs_Control_H898 = H898_Wound - H898_Control,
  Gallery_vs_Wound_H898 = H898_Gallery - H898_Wound,
  Gallery_vs_Control_Q903 = Q903_Gallery - Q903_Control,
  Wound_vs_Control_Q903 = Q903_Wound - Q903_Control,
  Gallery_vs_Wound_Q903 = Q903_Gallery - Q903_Wound,
  Control_H898_vs_Q903 = H898_Control - Q903_Control,
  Gallery_H898_vs_Q903 = H898_Gallery - Q903_Gallery,
  Wound_H898_vs_Q903 = H898_Wound - Q903_Wound,
  levels = linear.design)

# Voom transformation
linear.v <- voom(y, linear.design, plot=TRUE)

# MDS analysis
plotMDS(linear.v, top=Inf)

# PCA analysis
# library(ggplot2)
# pca <- as.data.frame(prcomp(t(linear.v$E))$x)
# ggplot(pca, aes(PC1, PC2, color = Group)) + geom_point(size=3) + 
#   scale_color_manual(name = "", 
#                      values=c("H898_Control" = "#4d004b", 
#                               "H898_Gallery" = "#88419d", 
#                               "H898_Wound" = "#8c96c6", 
#                               "Q903_Control" = "#7F0000", 
#                               "Q903_Gallery" =  "#d7301f", 
#                               "Q903_Wound" = "#fc8d59")) + 
#   labs(title="Voom + Limma Principal Component Analysis") + theme_bw()

# Linear modelling
linear.fit <- lmFit(linear.v, linear.design)

linear.fit2 <- contrasts.fit(linear.fit, linear.cont.matrix)

linear.fit2 <- eBayes(linear.fit2)

linear.results <- decideTests(linear.fit2, method="separate", 
                              adjust.method="BH", p.value=0.01)

summary(linear.results)
