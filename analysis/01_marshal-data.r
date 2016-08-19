library(purrr)
library(stringr)
library(testthat) # facilitate tests that will catch changes on re-analysis


# Marshal output files from Sailfish
rawFiles <- list.files("data/Sailfish-results-8June2016/", full.names = TRUE)
test_that("Exactly 24 Sailfish output files are found",
          expect_equal(24, length(rawFiles)))

## Extract, e.g. 898C1 from Sailfish-results/898C1_quant.sf
tmp <- str_extract(basename(rawFiles), "\\d{3}[A-Z]{1}\\d")

## Prepend "H" to 898 and "Q" to 903
tmp <- str_replace(tmp, "898", "H898")
tmp <- str_replace(tmp, "903", "Q903")

## Use as names for good side effects later
names(rawFiles) <- tmp


## Read one file to learn how many rows we expect and to grab rownames
tmp <- read.table(rawFiles[1], row.names = 1, header = TRUE,
                  ## specifying colClasses speeds this up 2x
                  colClasses = rep(c("character", "numeric"), c(1, 4)))
(n <- nrow(tmp)) # 108152
rwNames <- rownames(tmp)
test_that("First Sailfish output file has 108152 rows",
          expect_equal(108152, n))


## Read in all Sailfish data
system.time(
  myBigList <- map(rawFiles, function(x) {
    content <- read.table(x, header = TRUE, colClasses = c("character", rep("numeric", 4)))

    # Only extract columns that will be use for analysis
    content <- select(content, Name, NumReads)
    # Extract library name from filename
    libName <- str_extract(x, "\\d{3}[A-Z]\\d")
    # Pad library name with proper prefix
    libName <- str_replace(libName, "898", "H898")
    libName <- str_replace(libName, "903", "Q903")

        # Rename colname to "CDS" and respective library name
    colnames(content) <- c("CDS", libName)
    return(content)
  })
)
#    user  system elapsed
#  21.878   0.062  21.938 
# purrr not much faster than original implementation with aaply
str(myBigList) # List of 24

# Flatten big list of 24 into single data frame
rawSailfishCounts <- flatten_df(myBigList)
str(rawSailfishCounts) # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	108152 obs. of  25 variables:


# We discovered more ribosomal RNA post-assembly using BLAST.  The following 
# line removes putative ribosomal RNA from the table.
rRNA <- scan("data/ribosomalRNAContaminated.txt", what = "character")
summary(rawSailfishCounts$CDS %in% rRNA)
#    Mode   FALSE    TRUE    NA's 
# logical  108118      34       0 

microbial <- scan("data/microbialContaminated.txt", what = "character")
summary(rawSailfishCounts$CDS %in% microbial)
# Mode   FALSE    TRUE    NA's 
# logical  108142      10       0 

fungal <- scan("data/fungalContaminated.txt", what = "character")
summary(rawSailfishCounts$CDS %in% fungal)
#    Mode   FALSE    TRUE    NA's 
# logical  104915    3237       0 

contaminants <- c(rRNA, microbial, fungal)
contaminants <- unique(sort(contaminants))
length(contaminants)
# [1] 3277
summary(rawSailfishCounts$CDS %in% contaminants)
# Mode   FALSE    TRUE    NA's 
# logical  104875    3277       0 


rawSailfishCounts <- rawSailfishCounts[!(rawSailfishCounts$CDS %in% contaminants), ]
str(rawSailfishCounts) # Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	104875 obs. of  25 variables
summary(rawSailfishCounts$CDS %in% contaminants)
#    Mode   FALSE    NA's 
# logical  104875       0 


write.table(rawSailfishCounts, "data/consolidated-Sailfish-results.15July.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


## write design matrix
libNames <- colnames(rawSailfishCounts)[-1] # exclude the first one which is colname for CDS
expDes <- data.frame(sample = libNames, 
                     gType = factor(str_extract(libNames, "[A-Z]\\d{3}"), levels = c("Q903", "H898")),
                     txCode = factor(str_sub(libNames, 5, 5), levels = c('C', 'W', 'G')),
                     # C = Control, W = Wounding, G = Gallery 
                     tx = factor(NA),
                     bioRep = as.numeric(str_sub(libNames, 6, 6)))

expDes$tx <- recode(expDes$txCode, 'C' = "Control", 'W' = "Wound", 'G' = "Gallery")

expDes$grp <- with(expDes, interaction(gType, tx))
expDes
str(expDes)

write.table(expDes, "data/White_Pine_Weevil_exp_design.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
