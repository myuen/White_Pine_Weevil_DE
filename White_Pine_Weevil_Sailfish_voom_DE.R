### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

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

Group <- factor(interaction(targets$Genotype, targets$Treatment))
targets <- cbind(targets, Group=Group)

# Linear design matrix, simplest approach
design <- model.matrix(~ 0 + Group)
colnames(design) <- levels(targets$Group)

# Create contrast matrix
cont.matrix <- makeContrasts(
  Control.H898_vs_Q903 = H898.Control - Q903.Control,
  Induced.H898_vs_Q903 = 
    (H898.Gallery - (H898.Wound - H898.Control)) - 
    (Q903.Gallery - (Q903.Wound - Q903.Control)),
  levels = design)

# Voom transformation
v <- voom(y, design, plot=TRUE)

# MDS analysis
plotMDS(linear.v, top=Inf)

# PCA analysis
library(ggplot2)
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
