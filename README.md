# ptlmapper
Useful genetic tools allowing to discriminate individuals by comparing the distribution of their phenotype and to map the genetic loci responsible for specific distributions.

## Description

In living organisms, cell-to-cell variability can have dramatic consequences. For example, in some individuals, the occurrence of rare cellular events is responsible for disease. This single-cell variability can be quantified. It can also be described statistically by considering quantitative traits as random variables that follow probability density functions, each single-cell having a probability to express a given value of a trait. It is now known that genotypes can shape the statistical properties of single-cell phenotypic traits. Finding genetic loci involved in the control of single-cell quantitative trait density functions is therefore a promising challenge. However, methods to map such loci have not yet been developed beyond the approaches based on QTL (Quantitative Trait Loci) mapping, which are focused on the main moments (mean and variance) of the distribution. We propose here a novel method that allows us : i) to discriminate individuals by comparing the distribution of their phenotype and ii) to map the genetic loci responsible for specific distributions. This method is based on a metric that compares individuals on the basis of their single-cell trait distributions. Individuals are then projected into a space vector, and multivariate analysis searches for discrimination based on the genotype. Our method was validated on a simulated dataset of an auto-regulated gene model and on a published yeast morphological traits dataset. We have observed that our method is of particular interest to detect loci that have a small effect on the shape of single-cell trait density function.
  
## Installation

To get the current development version from github:

```R
install.packages("devtools")
devtools::install_github("fchuffar/ptlmapper")
```

## Usage

To learn more about how to use this method, refers to the vignette of the package: 

```R
browseVignettes("ptlmapper")
```
