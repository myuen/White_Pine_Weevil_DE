# read here about rules with multiple outputs
# http://www.cmcrossroads.com/article/rules-multiple-outputs-gnu-make?page=0%2C0
# but none of it seemed perfect for me, so I just wrote two identical rules :(

#! /usr/local/bin/Rscript --vanilla --default-packages=utils

# aggregates Sailfish-processed data
../data/White_Pine_Weevil_exp_design.tsv: ../data/Sailfish-results/*_quant.sf 01_marshal-data.r
	Rscript -e 'library(methods); knitr::stitch_rmd("01_marshal-data.r")'
	rm 01_marshal-data.html

# writes experimental design to proper file (same rule as above)
../data/consolidated-Sailfish-results.txt: ../data/Sailfish-results/*_quant.sf ../data/putativeRibosomalRNA.id 01_marshal-data.r
	Rscript -e 'library(methods); source("01_marshal-data.r")'

# filters out some contigs with practically no reads
../data/consolidated-filtered-Sailfish-results.txt: ../data/consolidated-Sailfish-results.txt 02_pre-dea-filtering.r
	Rscript -e 'library(methods); source("02_pre-dea-filtering.r")'

# figure documenting the filtering (same rule as above)
figure/White_Pine_Weevil_preDE_filtering.png: ../data/consolidated-Sailfish-results.txt 02_pre-dea-filtering.r
	Rscript -e 'library(methods); source("02_pre-dea-filtering.r")'

# conduct differential expression analysis and output statistical inference for
# the model parameters
../results/limma-results-model-terms.tsv: ../data/consolidated-filtered-Sailfish-results.txt ../data/White_Pine_Weevil_exp_design.tsv helper01_load-counts.r helper02_load-exp-des.r 03_dea-with-limma-voom.r
	Rscript -e 'library(methods); source("03_dea-with-limma-voom.r")'