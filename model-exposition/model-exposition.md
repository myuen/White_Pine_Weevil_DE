---
title: "Stating our model"
author: "Jenny Bryan"
date: "23 May, 2014"
output: html_document
---

We need to decide how to parametrize our model and, within that, to specify which effects we are interested in, i.e. which effects do we want to test for equality with zero? Here JB is using written and verbal material to work out how to translate the biological question(s) into statistical ones.

From Mack's original email:

"We have 24 RNA-Seq libraries coming from 2 genotypes: Susceptible (S) and Resistance (R), and 3 treatment conditions: Control (C), Gallery (G) and Wound (W), each with biological replicates.  We are interested in comparing: 1) the differential expression for each condition between the genotypes and (e.g. susceptible-gallery vs. resistance-gallery) 2) the differential expression between the treatment condition within the genotypes (e.g. susceptible-gallery vs. susceptible-wound)."

First let me just work out how to link to my own PNGs.

![](https://github.com/myuen/White_Pine_Weevil_DE/blob/master/model-exposition/model-exposition.001.png)

![](../master/model-exposition/model-exposition.001.png)

![](model-exposition.001.png)

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

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r}
summary(cars)
```

You can also embed plots, for example:

```{r, echo=FALSE}
plot(cars)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
