White_Pine_Weevil_DE
====================

White Pine Weevil Differential Expression Analysis Experiment

People involved in this particular analysis:

  * Mack Yuen <myuen@mail.ubc.ca>
  * Justin Whitehill <whiteh5@msl.ubc.ca>
  * Jenny Bryan <jenny@stat.ubc.ca>

From Mack's original email:

"We have 24 RNA-Seq libraries coming from 2 genotypes: Susceptible (S) and Resistance (R), and 3 treatment conditions: Control (C), Gallery (G) and Wound (W), each with biological replicates.  We are interested in comparing: 1) the differential expression for each condition between the genotypes and (e.g. susceptible-gallery vs. resistance-gallery) 2) the differential expression between the treatment condition within the genotypes (e.g. susceptible-gallery vs. susceptible-wound)."

Jenny listing the comparisons of interest:
  * for condition in {control, gallery, wound}, compare susceptible to resistant
  * for genotype in {susceptible, resistant}, compare treatments

[1] Create a new factor with 6 levels using interaction() with factors genotype and tx condition as inputs. Use y ~ x - 1 to force a cell means parametrization of the linear model. Make the contrasts you want "by hand". Literally, make a contrast matrix C. Do matrix multiplication of that times the estimated parameter vector to get estimates. Extract the estimated variance-covariance matrix of the parameter vector from the fitted model. Use the C^{T}VC formula to get estimated variance-covariance matrix for your vector of contrasts. Then take the square root of diagonal elements to get estimated standard errors. Make t statistics and then p-values.

[2] Use the "fit this one-way ANOVA model within each level of a factor" approach via a formula like y ~ a/b and y ~ b/a, where a and b are genotype and tx condition. Between the two fitted models -- which are, at a high level, equivalent! -- you will have estimated parameters that correspond exactly to the comparisons you care about. Therefore, you can "mine" the usual test statistics and p-values to get your inferential results.

Honestly, I would do it both ways to make sure I was getting the same results! At least for a couple of genes. Then use one approach or the other to scale up to all genes.

I think the cheat sheet I wrote for myself should be very helpful to you, now that I understand what you want. It's so similar!

I've paged through Applied Linear Models Fourth Edition now. Leads on sections especially relevant to you:

"ANOVA Model I -- Cell Means Model" is equation 16.2, page 673

"ANOVA Model I -- Factor Effects Model" is equation 16.62, page 693
It is amplified in section 16.11 on p. 696 and following.
This is the "effects sum to zero" approach.

*No where* that I can find do they discuss the parametrization that is the default in R. Namely intercept = reference level and other parameters are effects of the non-reference levels. Strange but true.

"Cell Means Model" in a two-factor study is described in eqn 19.15, p. 814 -->

The section 20.3 Analysis of Factor Effects when Interactions are Important looks relevant to you, especially if you choose option [1] and fit cell means.

(Chapters 21 Two-Factor Studies -- One Case Per Treatment and 22 Two-Factor Studies -- Unequal Sample Sizes and Unequal Treatment Importance … there but for the grace of God go you. Luckily I don't think you need to worry about those very real concerns.)

