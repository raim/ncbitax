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

## Use of this Package

This package allows to browse the NCBI taxonomy and extract sub-trees
via sets of NCBI taxonomy IDs, eg. from a local blast database or
genome-scanning results. The results can then be interpreted
along the taxonomic tree or alignment-based phylogenetic trees,
mapped to NCBI taxonomy IDs.

See the bash script `scripts/setup.sh` how to download NCBI taxonomy,
database, the Silva 16S rRNA tree and a protein consensus tree from 
cyanobacteria and map all trees to NCBI taxonomy identifiers.


## TODO

* improve speed by using more efficient indexed tree searches, eg.
via separate package or Rcpp
* convert `scripts/setup.sh` to an Rmarkdown vignette for the
R package, and move shell, Perl and R scripts to R package location
for scripts

* generate a small sample taxonomy to allow examples to run


* getID: catch missing and resolve multiple hits
* add current NCBI taxonomy as data object or generate a subset,
eg. cyanobacteria
* provide download & unpack function
* check whether functionalities are present in taxize package
* parse other formats, eg. newick, phylo object
* adapt and copy cmpTaxa.R and tax2newick.R rscripts from gIScanner
to here
* add a high-level function: mapTaxa to find the shortest
paths from query to target taxon IDs
