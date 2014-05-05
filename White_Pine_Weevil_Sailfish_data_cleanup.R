# Read the raw output file from Sailfish
H898C1 <- read.table("Sailfish-results/898C1_quant.sf", header = FALSE)
H898C2 <- read.table("Sailfish-results/898C2_quant.sf", header = FALSE)
H898C3 <- read.table("Sailfish-results/898C3_quant.sf", header = FALSE)
H898C4 <- read.table("Sailfish-results/898C4_quant.sf", header = FALSE)

H898G1 <- read.table("Sailfish-results/898G1_quant.sf", header = FALSE)
H898G2 <- read.table("Sailfish-results/898G2_quant.sf", header = FALSE)
H898G3 <- read.table("Sailfish-results/898G3_quant.sf", header = FALSE)
H898G4 <- read.table("Sailfish-results/898G4_quant.sf", header = FALSE)

H898W1 <- read.table("Sailfish-results/898W1_quant.sf", header = FALSE)
H898W2 <- read.table("Sailfish-results/898W2_quant.sf", header = FALSE)
H898W3 <- read.table("Sailfish-results/898W3_quant.sf", header = FALSE)
H898W4 <- read.table("Sailfish-results/898W4_quant.sf", header = FALSE)

Q903C1 <- read.table("Sailfish-results/903C1_quant.sf", header = FALSE)
Q903C2 <- read.table("Sailfish-results/903C2_quant.sf", header = FALSE)
Q903C3 <- read.table("Sailfish-results/903C3_quant.sf", header = FALSE)
Q903C4 <- read.table("Sailfish-results/903C4_quant.sf", header = FALSE)

Q903G1 <- read.table("Sailfish-results/903G1_quant.sf", header = FALSE)
Q903G2 <- read.table("Sailfish-results/903G2_quant.sf", header = FALSE)
Q903G3 <- read.table("Sailfish-results/903G3_quant.sf", header = FALSE)
Q903G4 <- read.table("Sailfish-results/903G4_quant.sf", header = FALSE)

Q903W1 <- read.table("Sailfish-results/903W1_quant.sf", header = FALSE)
Q903W2 <- read.table("Sailfish-results/903W2_quant.sf", header = FALSE)
Q903W3 <- read.table("Sailfish-results/903W3_quant.sf", header = FALSE)
Q903W4 <- read.table("Sailfish-results/903W4_quant.sf", header = FALSE)

# We are only taking column 7 (i.e. estimate numbers of reads) for 
# downstream DE analysis
rawSailfishCounts <- cbind(H898C1[,7], H898C2[,7], H898C3[,7], H898C4[,7],
                           H898G1[,7], H898G2[,7], H898G3[,7], H898G4[,7],
                           H898W1[,7], H898W2[,7], H898W3[,7], H898W4[,7],
                           Q903C1[,7], Q903C2[,7], Q903C3[,7], Q903C4[,7],
                           Q903G1[,7], Q903G2[,7], Q903G3[,7], Q903G4[,7],
                           Q903W1[,7], Q903W2[,7], Q903W3[,7], Q903W4[,7])

# All result files from Sailfish are sorted in the same order.  We are 
# porting the row name from one of the result file to our consolidated 
# results file
rownames(rawSailfishCounts) <- H898C1[,1]

# Simplify column names
# Demystify Column Names:
# First 4 alphanumeric character represents the susceptible/resistance genotype;
# H898 = Resistance Genotype; Q903 = Susceptible Genotype
#
# Followed by single alphabet for treatment;
# C = Control, G = Gallery, W = Wounding
#
# Last number in column name represents the biological replicates.
#
colnames(rawSailfishCounts) <- c("H898C1", "H898C2", "H898C3", "H898C4", 
                                 "H898G1", "H898G2", "H898G3", "H898G4", 
                                 "H898W1", "H898W2", "H898W3", "H898W4", 
                                 "Q903C1", "Q903C2", "Q903C3", "Q903C4", 
                                 "Q903G1", "Q903G2", "Q903G3", "Q903G4", 
                                 "Q903W1", "Q903W2", "Q903W3", "Q903W4")

# We discovered more ribosomal RNA post-assembly using BLAST.  The following 
# line removes putaive ribosomal RNA from the table.
rRNA <- read.table("putativeRibosomalRNA.id", header = FALSE)
rawSailfishCounts <- rawSailfishCounts[!(rownames(rawSailfishCounts) %in% rRNA[,1]),]

write.table(rawSailfishCounts, "consolidated-Sailfish-results.txt", 
            quote = FALSE, sep = "\t")