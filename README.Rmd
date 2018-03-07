
[![Travis-CI Build Status](https://travis-ci.org/bradfordcondon/WGPAN.svg?branch=master)](https://travis-ci.org/bradfordcondon/WGPAN)
[![DOI](https://zenodo.org/badge/64423500.svg)](https://zenodo.org/badge/latestdoi/64423500)

# WGPAN Whole-Genome Pairwise Alignment Noising
written by Bradford Condon, under advisement of Mark Farman
University of Kentucky


Special thanks to Hadley Wickham for writing [a fantastsic guide to distributing R packages](http://r-pkgs.had.co.nz/intro.html).

Note: this module is currently under development, converting it from scripts to a simple standalone package.   This package uses`formatR::tidy_dir("R")` to format code. Apologies about the inconsistent camel_case and snakeCase usaage, I was but a fledgeling programmer.


## Introduction

The WGPAN (Whole-Genome Pairwise Alignment Noising) package was designed to deal with neighbor-joining trees built from whole-genome pairwise alignment.  Traditional bootstrapping methods are invalid, because there is not a single alignment of all genomes to sample with replacement from.  As such, a new method for assigning confidence to nodes in a tree was required.  WGPAN allows you to assign confidence values to nodes that are relevant for **within clade confidence only**.  You can therefore test hypotheses regarding user-defined clades, **not** the relationship among those clades.

In addition, I've included some general utilities to make life a little easier when manipulating trees with the excellent package `APE`.  The vignette and documentation is forthcoming for these methods. 

# Installation

#### Using github
```
library(devtools)
install_github("bradfordcondon/WGPAN")
```
This package is not yet accepted in CRAN or BioConductor.  Once it is, instructions on installing it this way will be made available.

#### Vignettes

For example usage of the bootstrapping methods, please see the [vignette](https://github.com/bradfordcondon/WGPAN/blob/master/vignettes/WGPAN-vignette.Rmd)
