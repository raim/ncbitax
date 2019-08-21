#!/usr/bin/Rscript
## map 16S rRNA to taxon (acc -> gi -> taxon)
library("ape")

## CHANGE THIS PATH TO RUN THIS CORRECTLY FROM taxonomy_2018.sh
tpath="/data/introns/prokaryotes_2018/taxomony"

version <- 111
version <- 132
spaceholder <- "XYZYX" # replaces with spaces in version 132

## input/output file names for version 132
treef  <- file.path(tpath, paste0("LTPs",version,"_SSU_tree.newick"))
gb2gif <- file.path(tpath, paste0("LTPs",version,"_SSU_gb2gi.csv"))
gi2txf <- file.path(tpath, paste0("LTPs",version,"_SSU_gi2tax.csv"))
otree  <- file.path(tpath, paste0("LTPs",version,"_SSU_tree_taxonIDs.newick"))
odat   <- file.path(tpath, paste0("LTPs",version,"_SSU_tree2tax.dat"))

## tree variants
quote <- "" # version 111
if ( version>=132 ) {
    quote <- "'"
    ## loading file version with white space holder
    treef <- sub("tree.newick", "tree_nospaces.newick", treef) 
}

tree <- read.tree(treef, quote=quote)
gb2gi <- read.csv(gb2gif,header=FALSE, stringsAsFactors=FALSE)
gi2tx<- read.csv(gi2txf,header=FALSE,sep="\t", stringsAsFactors=FALSE)

## TODO: include taxons in xf below
colnames(gb2gi) <- c("acc","v","gi")
colnames(gi2tx) <- c("gi","tax")

if ( version==111 ) {
    range <- regexpr("__.*__",tree$tip.label)
    acc <- regmatches(tree$tip.labe,range)
acc <- gsub("^_+","",acc)
    acc <- gsub("_+.*","",acc)
}
if ( version==132 ) {
    acc <- gsub("_.*","", tree$tip.label)
}
ac2tx <- rep("",length(acc))
names(ac2tx) <- acc
for ( i in 1:length(ac2tx) ) {
  gi <- gb2gi[gb2gi[,"acc"]==acc[i],"gi"]
  ac2tx[i] <- gi2tx[gi2tx[,"gi"]==gi,"tax"]
}

## replace white-space holder in tree labels = species names
specids <- gsub(spaceholder, " ", tree$tip.label)


## store names and replace by taxonomy IDs
t.names <- cbind(names(ac2tx), specids, ac2tx)
colnames(t.names) <- c("acc","id","taxon")
tree$tip.label <- ac2tx

write.tree(tree, otree)
write.table(t.names, odat, row.names=FALSE,quote=F,sep="\t",col.names=TRUE)

if ( !interactive() ) quit(save="no")
