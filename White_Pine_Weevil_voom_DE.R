### Differential Expression Analysis on Sitka Spruce Weevil Experiment with limma + voom

# Load edgeR library
library(edgeR)

rawCounts <- read.table("RSEM_raw_counts_cleaned.txt", header=TRUE)

# Load counts into DGEList object from edgeR package.
y <- DGEList(counts=rawCounts)

# Keep only genes with at least 1 count-per-million reads (cpm) in 
# at least 4 samples.
z <- y[(rowSums(cpm(y) > 1) >= 4), ]

# Library depth is now changed with the filtering of the low count 
# contigs so we need to reset the libray depth.
z$samples$lib.size <- colSums(z$counts)

# TMM Normalization by Depth
z <- calcNormFactors(z)

# Create design matrix
targets <- cbind(colnames(rawCounts))
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

# The design generates the exact same design matrix as the linear approach
#treatmentPerGenotype.design <- model.matrix(~ Treatment %in% Genotype - 1, data=targets)
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


# Trying out a 2x3 factor design matrix approach
twoFactorModel.design1 <- model.matrix(~ Genotype/Treatment - 1, targets)
colnames(twoFactorModel.design1) <- c("H898", "Q903", "H898_Gallery", "Q903_Gallery",
                                      "H898_Wound", "Q903_Wound")

twoFactorModel.design2 <- model.matrix(~ Treatment/Genotype - 1, targets)
colnames(twoFactorModel.design2) <- c("Control", "Gallery", "Wound", 
                                      "Control_Q903", "Gallery_Q903", "Wound_Q903")

# Voom transformation
linear.v <- voom(z, linear.design, plot=TRUE)

twoFactorModel1.v <- voom(z, twoFactorModel.design1, plot=TRUE)
twoFactorModel2.v <- voom(z, twoFactorModel.design2, plot=TRUE)

# MDS analysis
#plotMDS(linear.v, top=Inf)

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

twoFactorModel1.fit <- lmFit(twoFactorModel1.v, twoFactorModel.design1)
twoFactorModel2.fit <- lmFit(twoFactorModel2.v, twoFactorModel.design2)

linear.fit2 <- contrasts.fit(linear.fit, linear.cont.matrix)

linear.fit2 <- eBayes(linear.fit2)
twoFactorModel1.fit2 <- eBayes(twoFactorModel1.fit)
twoFactorModel2.fit2 <- eBayes(twoFactorModel2.fit)

#linear.results <- decideTests(linear.fit2, method="global", adjust.method="BH", p.value=0.01)

linear.results <- decideTests(linear.fit2, 
                              adjust.method="BH", p.value=0.01)
twoFactorModel1.results <- decideTests(twoFactorModel1.fit2, 
                                       adjust.method="BH", p.value=0.01)
twoFactorModel2.results <- decideTests(twoFactorModel2.fit2, 
                                       adjust.method="BH", p.value=0.01)

summary(linear.results)
summary(twoFactorModel1.results)
summary(twoFactorModel2.results)
