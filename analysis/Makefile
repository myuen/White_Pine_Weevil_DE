# I read here about rules with multiple outputs
# http://www.cmcrossroads.com/article/rules-multiple-outputs-gnu-make?page=0%2C0
# but none of it seemed perfect for me, so I often wrote two identical rules :(

# the testthat package appears to depend on methods but currently does not
# declare that properly; my weird Rscript calls below address that; watch this
# issue for developments: https://github.com/hadley/testthat/issues/122

# aggregates Sailfish-processed data
../data/White_Pine_Weevil_exp_design.tsv: ../data/Sailfish-results/*_quant.sf 01_marshal-data.r
	Rscript -e 'library(methods); knitr::stitch_rmd("01_marshal-data.r")'
	rm 01_marshal-data.html

# writes experimental design to proper file (same rule as above)
../data/consolidated-Sailfish-results.txt: ../data/Sailfish-results/*_quant.sf ../data/putativeRibosomalRNA.id 01_marshal-data.r
	Rscript -e 'library(methods); knitr::stitch_rmd("01_marshal-data.r")'
	rm 01_marshal-data.html

# filters out some contigs with practically no reads
../data/consolidated-filtered-Sailfish-results.txt: ../data/consolidated-Sailfish-results.txt 02_pre-dea-filtering.r
	Rscript -e 'library(methods); knitr::stitch_rmd("02_pre-dea-filtering.r")'
	rm 02_pre-dea-filtering.html

# figure documenting the filtering (same rule as above)
figure/White_Pine_Weevil_preDE_filtering.png: ../data/consolidated-Sailfish-results.txt 02_pre-dea-filtering.r
	Rscript -e 'library(methods); knitr::stitch_rmd("02_pre-dea-filtering.r")'
	rm 02_pre-dea-filtering.html

# conduct differential expression analysis and output statistical inference for
# the model parameters
../results/limma-results-model-terms.tsv: ../data/consolidated-filtered-Sailfish-results.txt ../data/White_Pine_Weevil_exp_design.tsv helper01_load-counts.r helper02_load-exp-des.r 03_dea-with-limma-voom.r
	Rscript -e 'library(methods); knitr::stitch_rmd("03_dea-with-limma-voom.r")'
	rm 03_dea-with-limma-voom.html 03_dea-with-limma-voom.Rmd

# documenting my work as I sorted out puzzling behavior re: column names of the
# model matrix and making contrasts
90_limma-model-term-name-fiasco.md: 90_limma-model-term-name-fiasco.r
	Rscript -e 'library(methods); knitr::spin("90_limma-model-term-name-fiasco.r", precious = TRUE)'
	rm 90_limma-model-term-name-fiasco.html 90_limma-model-term-name-fiasco.Rmd
