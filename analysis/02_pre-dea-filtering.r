library(edgeR)
library(ggplot2)
library(plyr)
library(testthat) # facilitate tests that will catch changes on re-analysis

### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# Load counts from Sailfish
rawSailfishCounts <- read.delim("../data/consolidated-Sailfish-results.txt")
str(rawSailfishCounts) # 'data.frame':  491928 obs. of  24 variables:
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
