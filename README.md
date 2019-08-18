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

* add current NCBI taxonomy as data object or generate a subset,
eg. cyanobacteria
* provide download & unpack function
* check whether functionalities are present in taxize package
* improve speed by using more efficient indexed tree searches, eg.
via separate package or Rcpp

