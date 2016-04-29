# ptlmapper
This package contains functions that are used to map single-cell Probabilistic Trait Loci (scPTL). The principle of scPTL mapping is to consider a distribution of single-cell traits as the phenotype of a (multicellular) individual, to acquire this phenotype in many individuals, and to scan the genome for DNA variants that change this phenotype.

## Description

Cellular quantitative traits can be described statistically by considering them as random variables that follow specific probability density functions. Each cell then has a given probability to express a given value of a trait. It is now known that genotypes can shape the statistical properties of single-cell phenotypic traits. Finding genetic loci involved in the control of single-cell quantitative trait density functions is therefore a promising area of investigation. We propose here a method that allows to i) discriminate individuals by comparing the distribution of their phenotype and ii) map the genetic loci responsible for differences in the distributions. This method is based on the Kantorovich metric that compares individuals on the basis of their single-cell trait distributions. Individuals are then projected into a vectorial space, and multivariate analysis is then applied to search for discrimination based on the genotype. Our method was validated on a simulated dataset of an auto-regulated gene model and on a published dataset of yeast morphological traits. See reference Chuffart et al. 2016 for details.
  
## Installation

To get the current development version from github:

```R
install.packages("devtools")
devtools::install_github("fchuffar/ptlmapper")
```

## Usage

To learn on how to use this method, please go to the vignette of the package: 

```R
browseVignettes("ptlmapper")
```

## Reference

Chuffart F, Richard M, Jost D, Burny C, Duplus-Bottin H, Ohya Y and Yvert G. Exploiting single-cell quantitative data to map genetic variants having probabilistic effects. PLoS Genetics 2016 (in revision)