White_Pine_Weevil_DE
====================

White Pine Weevil Differential Expression Analysis Experiment

People involved in this particular analysis:

  * Mack Yuen <myuen@mail.ubc.ca>
  * Justin Whitehill <whiteh5@msl.ubc.ca>
  * Jenny Bryan <jenny@stat.ubc.ca>

From Mack's original email:

"We have 24 RNA-Seq libraries coming from 2 genotypes: Susceptible (S) and Resistance (R), and 3 treatment conditions: Control (C), Gallery (G) and Wound (W), each with biological replicates.  We are interested in comparing: 1) the differential expression for each condition between the genotypes and (e.g. susceptible-gallery vs. resistance-gallery) 2) the differential expression between the treatment condition within the genotypes (e.g. susceptible-gallery vs. susceptible-wound)."

Goals of the differential expression analysis, in rough order in interest:

  * In the Control condition, compare the susceptible Q903 strain to the resistant H898
    - Justin has greater interest in "things that are unique in H898", i.e. expression in H898 >> Q903
  * In each genotype, compare the Wound treatment to the Gallery
    - Justin has greater interest in "things that are unique to Gallery", i.e. expression under Gallery tx >> Wound tx
  * For each condition, compare the susceptible Q903 strain to the resistant H898
    - JB musing: make a scatterplot matrix of these effects?
  * For each genotype, resistant}, compare treatments
    - JB musing: again, scatterplot matrix?

This contrast defined in Mack's original code conveys what they call "induced", which seems to get at the phenomenon of expression variation in the Gallery condition that is above and beyond changes seen in the Wound condition (expression below presumes a cell means type of parametrization):

```
Induced.H898_vs_Q903 = 
    (H898.Gallery - (H898.Wound - H898.Control)) - 
    (Q903.Gallery - (Q903.Wound - Q903.Control))
```

The main reference for `limma + voom` for the analysis of RNA-seq data:

Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology, 15(2), R29. doi:10.1186/gb-2014-15-2-r29 [link](http://genomebiology.com/2014/15/2/R29)

The [`limma` User's Guide](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) has Chapter 15 on RNA-seq data with code.

Look at this: [Case study: using a Bioconductor R pipeline to analyze RNA-seq data](http://bioinf.wehi.edu.au/RNAseqCaseStudy/).
