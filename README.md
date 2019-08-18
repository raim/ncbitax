# `ncbitax` - a simple NCBI taxonomy parser

A very light-weight straight-forward parser for NCBI 
taxonomy, for simple lineage retrieval, and common ancestor searches.

## Installation, Quick Start & Usage

### With `devtools` from git

### Manual from git 

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
