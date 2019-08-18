# `ncbitax` - a simple NCBI taxonomy parser

A very light-weight straight-forward parser for NCBI 
taxonomy, for simple lineage retrieval, and common ancestor searches.

## Installation, Quick Start & Usage

### With [`devtools`](https://cran.r-project.org/package=devtools) from git

```
library(devtools)
install_github("raim/ncbitax")
```

### Source code from `github`

* Download:

```
git clone git@github.com:raim/ncbitax.git
```
* Install:

```
R CMD build ncbitax
sudo R CMD INSTALL ncbitax_0.0.1.tar.gz
```

## TODO

* improve speed by using more efficient indexed tree searches, eg.
via separate package or Rcpp

* add current NCBI taxonomy as data object or generate a subset,
eg. cyanobacteria
* provide download & unpack function
* check whether functionalities are present in taxize package
* parse other formats, eg. newick, phylo object
* adapt and copy cmpTaxa.R and tax2newick.R rscripts from gIScanner
to here
* add a high-level function: mapTaxa to find the shortest
paths from query to target taxon IDs
