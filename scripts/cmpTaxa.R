#!/usr/bin/Rscript

## calculates taxonomic relations
## USAGE eg.
#datpath="/home/raim/data/introns/prokaryotes_2018/taxonomy"
#outpath="/home/raim/work/CoilProject/halime/taxonomy/species2silva"
#infile="/data/topoisomerases/cyanoids.txt"
#
#scripts/cmpTaxa.R  -txd $datpath  -tn acc -tf ${datpath}/LTPs132_SSU_tree2tax.dat -sf $infile -out $outpath -xf ${outpath}/species2silva_tree.RData -dbug 

## get directory of utils
library(ncbitax)

## parse commandline into normal variables
# logical parameters, override with --<parameter>, e.g, --inc
dbug <- FALSE
sub <- FALSE

args  <- R.utils::commandArgs(excludeReserved=TRUE, asValues=TRUE)
for ( i in 2:length(args) ) {
  arg <- names(args)[i]
  #cat(paste(arg, as.character(args[arg]), "\n"))
  assign(arg, as.character(args[arg]))
}

#tf <- "taxonomy/phispy_taxons.dat"
#sf <- "genomes/bacteria/taxonomy.dat"
#xf <- "taxonomy/phispyXspecies.phy"

if ( !exists("tf",mode="character") ) {
  cat(paste("no target taxons provided, use --tf=<filename>\n"))
  if ( !interactive() ) quit(save="no")
}
if ( !exists("sf",mode="character") ) {
  cat(paste("no source taxons provided, use --sf=<filename>\n"))
  if ( !interactive() ) quit(save="no")
}
#if ( !exists("txd",mode="character") ) {
#  cat(paste("no taxononmy dir, use -txd=<dir>\n"),file=stderr())
#  if ( !interactive() ) quit(save="no")
#}
#if ( !exists("xf",mode="character") ) xf <- ""
if ( !exists("out",mode="character") ) out <- "./sp2ph"
if ( !exists("tn",mode="character") ) tn <- "" ## name, acc
dbug <- as.logical(dbug)
sub <- as.logical(sub)

## check if either taxonomy newick tree or NCBI taxonomy directory
## are present

if ( !exists("txd",mode="character") & !exists("xf",mode="character") )
    stop("either newick tree or directory with NCBI taxonomy files required!")


## find the closest relatives between two lists of ncbi taxon ids

## map species to 16S rRNA,
## where species is not present, load taxonomy tree to find closest relative
## that is present in 16S rRNA

### PARSE QUERY AND TARGET TAXON IDs

targt <- read.csv(tf, header=TRUE, sep="\t")
sourc <- read.csv(sf, header=TRUE, sep="\t")

if ( !"taxon"%in%colnames(targt) | !"taxon"%in%colnames(sourc) )
    stop("both query and target require a column with header `taxon'")

## only one per source taxa
sourc <- sourc[!duplicated(sourc[,"taxon"]),,drop=FALSE]

## taxon IDs, can be several!?
targtx <- lapply(strsplit(as.character(targt[,"taxon"]),";"),
                 function(x) gsub(" +","",x))
srcetx <- lapply(strsplit(as.character(sourc[,"taxon"]),";"),
                 function(x) gsub(" +","",x))

## source taxon IDs
allsrc <- unlist(srcetx)
## several taxon IDs, per target
alltar <- unlist(targtx) 

## combined list of taxon IDs
taxids <- unique(c(alltar,allsrc))


## READ OR GENERATE TAXONOMY

## read or generate common taxonomy tree
## get tree for all taxa from NCBI taxonomy, including singletons!
## (cutting the full taxonomy tree here saves lots of time below)


## TODO: write/read newick tree instead
##       * the written tree can not be parsed!?
##       read tree with singleton nodes?
##       * rank information missing from tree, only provided by
##       taxonomy parser

## check if an RData file with the tree exists already
xfile <- ""
if ( exists("xf",mode="character") )
    if ( file.exists(xf) )
        xfile <- xf

if ( xfile != "" ) {

    cat(paste("loading previously generated common taxonomy tree:",xfile,"\n"))

    ## load previously parsed taxonomy as RData
    ## TODO: load tree
    load(xfile)
    
    ##stop("TREE PARSING NOT IMPLEMENTED YET",
    ##     "please provide taxonomy directory")
    ##require("ape")
    ##tree <- read.tree(file=xf)

} else if ( exists("txd",mode="character") ) {

    cat(paste("generating common taxonomy tree, can take long:",txd,"\n"))

    ## generate common taxomony tree for all taxonomy IDs
    taxo <- tax2newick(taxids, txd, names=FALSE,
                       ranks=TRUE, full=TRUE,
                       verb=dbug, dbug=dbug)
    ## store as RData
    ## if xf was provided but didn't exist as an RData file yet
    ## TODO: store as newick tree such that parsing above works!?
    if ( exists("xf",mode="character") ) {
        ##require("ape");
        ##write.tree(tree, file=xf)
        save(taxo, file=xf)
    } 
}

## get tree and rank from tax2newick object
tree <- taxo$tree
rank <- taxo$rank
## get all nodes. the order is used as index in tree$edge
nodes <- c(tree$tip.label,tree$node.label)


## TODO: generate new NCBI tax object, and use this below instead
## of tree$edges
ntax <- list()
ntax$parents <- edg2idx(taxo$tree)
ntax$rank <- taxo$rank

#warnings(file=stderr()) ## reports missing taxon ids

## filter not found nodes
if ( sum(!alltar %in% nodes)>0 )
  cat(paste("\nMISSING TARGET SPECIES:",
            paste(alltar[!alltar %in% nodes],collapse=";"), "\n"))
alltar <- alltar[alltar %in% nodes]

cat(paste("taxonomy loaded; searching closest neighbors ..."),file=stderr())


## get closest target species for each source
## (and collect names, IDs and common ancestor)
src2tar <- src2nam <- src2tax <- rep(list(NA), length(allsrc))
src2anc <- s2t.dst <- reltype <- rep(NA,length(allsrc))

## node index of target taxa
idst <- unique(alltar)
tid <- rep(NA,length(idst))
for ( j in 1:length(idst) )
    if ( idst[j]%in%nodes ) {
        tid[j] <- which(nodes==idst[j])
    } else {warning("target", j, ", taxon ID", idst[j],
                    " not found in nodes; obsolete tax id?")}

for ( i in 1:length(allsrc) ) {
    
    if ( dbug )
        cat(paste0("searching targets for taxon #",
                   i, "/", round(i/length(allsrc),2),
                   ", ID ", allsrc[i],"\t"))
    
    ## first, check if query taxon is present at all in taxonomy
    ## and get numeric index of its node, as used in tree$edge
    ## (note: only required for option sub, which we may skip alltogether)
    nid <- which(nodes==allsrc[i])
    if ( length(nid)!=1 ) {
        warning("taxon ID ", allsrc[i], " not found in nodes; obsolete tax id?")
        next
    }
    
    
    ## 1) is any target equal to the source?
    if ( allsrc[i] %in% alltar ) {
        
        src2tax[[i]] <- src2anc[i] <- allsrc[i]
        s2t.dst[i] <- 0
        reltype[i] <- "ident"
        if ( dbug ) cat("found direct hit")
    } else {
        
        ## if not, scan taxonomy tree
        
        ## get target taxa and ancestor
        ## via the lowest branch between src <-> target
        pars <- get.parents(nid, ntax) 
        
        ## 2) is any target a direct ancestor? <- find parent species of strains
        if ( any(pars%in%tid) )  {
            
            ## take the first parent!
            src2tax[[i]] <- src2anc[i] <-
                nodes[pars[pars%in%tid]][1] ## ancestor=target
            s2t.dst[i] <- which(pars%in%tid)[1] - 1 ## taxonomic distance
            reltype[i] <- "parent"
            if ( dbug ) cat("found parent")
            
        } else {
            
            ## 3) scan children of ancestors <- find sibling strains
            ## TODO: record siblings, even if parent is present?
            
            ## to save time ??, extract a subtree
            ## NOTE: this doesn't save time apparently
            ## TODO: what is root here?
            if ( sub ) 
                st <- getAllParents(c(nid,tid),
                                    ntax, root=length(tree$tip.label)+1,
                                    dbug=FALSE)
            else st <- ntax
            
            ## note: first "parent" is query taxon
            for ( j in 2:length(pars) ) {
                skip <- pars[1:(j-1)] # skip prev. to save time!
                cdr <- get.children(pars[j], st, skip=skip)
                present <- cdr %in% tid
                if ( any(present) ) {
                    src2tax[[i]] <- nodes[cdr[present]] #nodes[tid[tid%in%cdr]] 
                    src2anc[i] <- nodes[pars[j]]
                    ## taxonomic distance to ancestor
                    ## TODO: add distance ancestor-sibling ?
                    s2t.dst[i] <- -j # -1 
                    reltype[i] <- "sibling"
                    if ( dbug ) cat(paste("found", sum(present),
                                          "sibling(s) at", -j))
                    break
                }
            }
        }
    }
    
    ## get names and rank of target taxa
    idx <- unlist(lapply(targtx, function(x) any(src2tax[[i]] %in% x)))
    if ( sum(idx)==0 )
        stop("taxon ", i, " ID ", alltar[i], ": NOTHING FOUND, SHOULDN'T",
             "HAPPEN - BUG??")
    
    src2tar[[i]] <- c(as.character(targt[idx,"id"]))
    if ( tn!= "" ) src2nam[[i]] <- c(as.character(targt[idx,tn]))
    
    ## if target has "no rank" get next ancestor with rank
    if ( rank[src2anc[i]] == "no rank" ) {
        id <- nodes[get.parents(id=which(nodes==src2anc[i]), tax)]
        idx <- rank[id] != "no rank"
        if ( any(idx) ) {
            ## take first with rank
            src2anc[i] <- id[which(idx)[1]]
            ##s2t.dst[i] <- s2t.dst[i] + which(idx)[1] - 1
        }
    }
    
    if ( dbug ) cat("; done\n")
}

trgs <- NULL
if ( tn != "" )
  trgs <- unlist(lapply(src2nam, function(x) paste(x,collapse=";")))
tids <- unlist(lapply(src2tar, function(x) paste(x,collapse=";")))
txids <- unlist(lapply(src2tax, function(x) paste(x,collapse=";")))

cat(paste("... done. writing results .... "), file=stderr())
results <- cbind.data.frame(sourc,
                            taxdist=s2t.dst,
                            relation=reltype,
                            target.taxon=txids,
                            ancestor=src2anc,
                            rank=rank[src2anc],
                            target=trgs,
                            tid=tids,
                            stringsAsFactors=FALSE)

## sort by distance to ancestor
#results <- results[order(as.numeric(results[,"taxdist"])),]

file.name <- paste(out,".dat",sep="")
if ( file.exists(file.name) ){
  cat(paste("WARNING: overwriting", file.name, "\n"),file=stderr())
  unlink(out)
}
write.table(results,file.name,row.names=FALSE,quote=FALSE,sep="\t")

cat(paste("... done.\n"),file=stderr())

if (!interactive() ) {
  file.name <- paste(out,".RData",sep="")
  save.image(file=file.name)
  quit(save="no")
}


