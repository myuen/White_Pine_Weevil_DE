### Differential Expression Analysis on Sitka Spruce Weevil Experiment with limma + voom

# Load counts from RSEM
rawCounts <- read.table("WPW_Inoculation_Trinity_C500.gene.counts.matrix", header = TRUE)

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
colnames(rawCounts) <- c("H898C1", "H898C2", "H898C3", "H898C4", 
                         "H898G1", "H898G2", "H898G3", "H898G4", 
                         "H898W1", "H898W2", "H898W3", "H898W4", 
                         "Q903C1", "Q903C2", "Q903C3", "Q903C4", 
                         "Q903G1", "Q903G2", "Q903G3", "Q903G4", 
                         "Q903W1", "Q903W2", "Q903W3", "Q903W4")

# We discovered more ribosomal RNA post-assembly using BLAST.  The following 
# line removes putaive ribosomal RNA from the table.
rRNA <- read.table("blast/putativeRibosomalRNA.id", header = TRUE)
rawCounts <- rawCounts[!(rownames(rawCounts) %in% rRNA$Query),]

write.table(rawCounts, "RSEM_raw_counts_cleaned.txt", quote=FALSE, sep="\t")
