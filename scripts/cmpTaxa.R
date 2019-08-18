#!/usr/bin/Rscript

## calculates taxonomic relations

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
sourc <- sourc[!duplicated(sourc[,"taxon"]),]

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
    taxo <- tax2newick(taxids, txd, species=FALSE,
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
        pars <- get.parents(nid, tree$edge) 
        
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
            if ( sub ) 
                st <- getAllParents(c(nid,tid),
                                    tree$edge, rootid=length(tree$tip.label)+1,
                                    dbug=FALSE)
            else st <- tree$edge
            
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
        id <- nodes[get.parents(id=which(nodes==src2anc[i]), tree$edge)]
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




#### NOTE 20180916:
#### NOT DONE: kept code in case it contains useful improvement attempts
#### MAY REQUIRE ADAPTATION TO NEW PARSING ABOVE!!!


## get taxonomic distance between all nodes
## TODO: extract subtree for source & targets only!
tree <- taxo$tree
taxd <- dist.nodes(tree)
nm <- c(tree$tip.label, tree$node.label)
dimnames(taxd) <- list(nm,nm)


cat(paste("... done\nfinding closest neighbors ..."),file=stderr())

## FIND CLOSEST SPECIES: mapping src -> target
if ( sum(!alltar %in% colnames(taxd))>0 ) {
  na <- paste(alltar[!alltar %in% colnames(taxd)],collapse=";")
  cat(paste("target taxon ids", na, "not available\n"),file=stderr())
  alltar <- alltar[alltar %in% colnames(taxd)]
}
src2tar <- src2nam <- src2tax <- rep(list(NA), length(allsrc))
s2t.dst <- rep(NA, length(allsrc))
for ( i in 1:length(allsrc) ) {
  src <- allsrc[i]
  if ( src %in% rownames(taxd) ) {
    ds <- taxd[src,alltar]
    s2t.dst[i] <- min(ds) ## get minimal distance
    s2t <- alltar[ds==min(ds)] ## get all target taxa with min dist
    ## get ids of target taxa
    ids <- lapply(s2t, function(x) as.character(targt[targt[,"taxon"]==x,"id"]))
    src2tar[[i]] <- unlist(ids)
    ## get taxa of target taxa
    src2tax[[i]] <- names(which(ds==min(ds)))
     ## get names of target taxa
    if ( tn!= "" ) {
      nms <- lapply(s2t,function(x) as.character(targt[targt[,"taxon"]==x,tn]))
      src2nam[[i]] <- unlist(nms)
    }
  }else cat(paste("source", src,"not available\n"),file=stderr())
  
}


trgs <- NULL
if ( tn != "" )
  trgs <- unlist(lapply(src2nam, function(x) paste(x,collapse=";")))
tids <- unlist(lapply(src2tar, function(x) paste(x,collapse=";")))
# paste(src2tar[[i]],collapse=";")

cat(paste("... done. getting ancestor and rank .... "),file=stderr())

## FIND COMMON ANCESTOR
anc <- rep(NA,length(allsrc))
for ( i in 1:length(allsrc) ){
  ids <- unique(c(allsrc[i], src2tax[[i]]))
  nids <- rep(NA,length(ids)) ## get numeric indices of nodes
  for ( j in 1:length(ids) ) 
    nids[j] <- which(c(tree$tip.label,tree$node.label)==ids[j])
  ### get common ancestor of all src<->target mappings
  #nanc <- getAncestor(nids, tree$edge, rootid=length(tree$tip.label)+1)
  #
  ### map back numeric index to node labels
  #anc[i] <- c(tree$tip.label,tree$node.label)[as.numeric(nanc)]

  st <- getAllParents(nids,tree$edge,dbug=FALSE,
                      rootid=length(tree$tip.label)+1)
  
  ## get lowest branch between src <-> target
  pars <- get.parents(nids[1],st)
  tid <- nids[2:length(nids)]
  ## is any target a direct ancestor?
  if ( any(tid %in% pars) ) {    
    src2tax[[i]] <- c(tree$tip.label,tree$node.label)[tid[tid%in%pars]]
    anc[i] <- src2tax[[i]]
  } else {
    ## scan children of ancestors
    for ( j in 2:length(pars) ) {
      cdr <- get.children(pars[j], st)
      if ( any(tid %in% cdr) ) {
        anc[i] <- c(tree$tip.label,tree$node.label)[pars[j]]
        src2tax[[i]] <- c(tree$tip.label,tree$node.label)[tid[tid%in%cdr]]
        break
      }
    }
  }
  ## if ancestor has "no rank" get next ancestor with rank from full tree
  if ( rank[anc[i]] == "no rank" ) {
    id <- get.parents(id=which(c(tree$tip.label,tree$node.label)==anc[i]),
                       tree$edge)
    id <- c(tree$tip.label,tree$node.label)[id]
    if ( any(rank[id]!="no rank") )
      anc[i] <- id[which(rank[id]!="no rank")[1]]
  }
  ## TODO: get lowest branch!? and replace both src2tax/nam/tar[[i]] and anc[i]
  ## 
  #cat(paste(i, "ancestor is", anc[i], rank[anc[i]], s2t.dst[i], "\n"))
}

cat(paste("... done. writing results .... "),file=stderr())
results <- cbind(taxdist=s2t.dst,
                 name=as.character(sourc[,"species"]),
                 taxon=sourc[,"taxon"],
                 target=trgs,
                 tid=tids,ancestor=anc,rank=rank[anc])

## sort by distance
results <- results[order(as.numeric(results[,"taxdist"])),]

if ( file.exists(out) ){
  cat(paste("WARNING: overwriting", out, "\n"),file=stderr())
  unlink(out)
}
write.table(results,out,row.names=FALSE,quote=FALSE,sep="\t")

cat(paste("... done.\n"),file=stderr())

if (!interactive() ) {
  save.image(file="cmpTaxa.Rdata")
  quit(save="no")
}



## get 16S rRNA distance: too big Error: cannot allocate vector of size 2.8 Gb
#treed <- cophenetic.phylo(tree)

## tax2spec 
tx2sp <- rep("",nrow(taxd))
names(tx2sp) <- rownames(taxd)
for ( i in 1:length(tx2sp) ) {
  idx <- which(spec[,"tax"]==names(tx2sp)[i])
  tx2sp[i] <- as.character(spec[idx[1],"species"])
}

## map 16S rRNA to taxon (acc -> gi -> taxon)
range <- regexpr("__.*__",tree$tip.label)
acc <- regmatches(tree$tip.labe,range)
acc <- gsub("^_+","",acc)
acc <- gsub("_+.*","",acc)
ac2tx <- rep("",length(acc))
names(ac2tx) <- acc
for ( i in 1:length(ac2tx) ) {
  gi <- as.character(gb2gi[as.character(gb2gi[,"acc"])==acc[i],"gi"])
  ac2tx[i] <- as.character(gi2tx[as.character(gi2tx[,"gi"])==gi,"tax"])
}

## map species taxa to 16S rRNA
sp2tr <- matrix(NA, length(tx2sp),ncol=2)
colnames(sp2tr) <- c("taxon","distance")
maxd <- 7#6#5#4
for ( i in 1:length(tx2sp) ) {
  tx <- names(tx2sp)[i]
  idx <- which(ac2tx==tx)
  if ( length(idx)==0 ) {
    cat(paste(i, tx2sp[i], "not found; scanning along lineage\n"))
    for ( j in 1:maxd ) {
      rel <- names(which(taxd[i,]==j))
      if ( length(rel)>0 ) {
        #cat(paste("found relatives at dist.", j, "\n"))
        if ( sum(ac2tx %in% rel)>0 ) {
          idx2 <- which(ac2tx %in% rel)
          cat(paste(sum(ac2tx %in% rel), "relatives at dist.", j, " in 16S\t"))
          cat(paste(paste(tree$tip.label[idx2],collapse=";"),"\n"))
          sp2tr[i,] <- c(paste(idx2,collapse=";"),j)
          break
        }
      }
    }
  } else if ( length(idx)>1 ) {
    cat(paste(i, tx2sp[i], "multiple found, taking first\n"))
    tx <- unique(ac2tx[idx])
    if ( length(tx)>1 ) 
      stop(paste(i, tx2sp[i], "conflicting multiple hits!\n"))
    sp2tr[i,] <- c(idx[1],0)
  }
  else sp2tr[i,] <- c(idx,0)
}
cat(paste(sum(is.na(sp2tr[,1])), "species could not be assigend\n"))


## for each species, find closest relative among a list of target
## organisms, each with one or more taxon ids

taxalist <- strsplit(sp2tr[,"taxon"],";")
targlist <- strsplit(as.character(targt[,"taxon"]),";")
#map <- matrix(NA, ncol=nrow(sp2tr), nrow=nrow(targt))
for ( i in 1:nrow(sp2tr) ) {
  for ( j in 1:nrow(targt) ) {
    taxa <- unlist(strsplit(sp2tr[i,"taxon"],";"))
    targ <- unlist(strsplit(as.character(targt[j,"taxon"]),";"))
    if ( sum(taxa %in% targ)>0 ) {
      tmp <- paste(taxa[taxa %in% targ],collapse=";")
      cat(paste(i,j,"found targets for", tmp, "\n"))
    }
  }
}




## map taxonomy IDs in bacteria to 16S rRNA tree
leaves <- rep(NA, nrow(spec))
for ( i in 1:nrow(spec) ) {
  gi <- as.numeric(spec[i,"gi"])
  taxid <- as.numeric(spec[i,"tax"])
  if ( !taxid %in% as.numeric(gi2tx[,2]) ) 
    cat(paste(i, "NOT found",spec[i,"seq"], "\n"))
  else {
    cat(paste(i, "found",spec[i,"seq"], "\n"))
    leaves[i] <- which(as.numeric(gi2tx[,2])==taxid)
  }
}

bacteria <- extract.clade(tree,"Bacteria")

saureus <- grep("Staphylococcus",bacteria$tip.label)
zoom(bacteria,saureus)

strep <-  extract.clade(tree,"Staphylococcaceae")
saureus <- grep("Staphylococcus_aureus",strep$tip.label)


#tax <- read.tree("taxonomy/ncbi_complete_with_taxIDs.newick")

## TODO:
## map 16s rRNA to tax id: AccNum in LTPs108_SSU.csv, AccNum2GI in GbAccList.0127.2013 and GI2taxid in taxonomy/gi_taxid_nucl.dmp
## OOR add Genbank ID to bacteria/taxonomy.dat 

