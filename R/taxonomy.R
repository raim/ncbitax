
## PHYLOGENY/TAXONOMY: TREE PARSING AND GENERATION
## similar to functions provided by ape and phylobase, but
## the functions here work on a named edge list
## ("tax") only

## TODO: remove singular root chain
trimRoot <- function(tax, ids) {}

## TODO: cut leaves, where one of the ids
trimLeaves <- function(tree, ids) {}

## TODO: add a node to an existing tree object, get edges via taxonomy
addNode <- function(tree, id, tax) {
  ## add lineage edges via getParents and
  ## check whether children (collapsed singletons!)
  ## exist
}

#' get lowest common ancestor of nodes 'ids'
#'
#' takes an NCBI taxonomy ID and the NCBI taxonomy tree, and reports
#' the lowest common ancestor of all taxa
#' TODO: either use tax (named edges) or rename and write function for tax
#' @param ids node IDs for which ancestors should be reported
#' @param tax taxonomy tree
#' @param rootid ID of the taxonomic class to be considered as root
#' @param dbug print debug messages
#' @export
getAncestor <- function(ids, tax, rootid="1", dbug=FALSE) {

  anc <- unique(ids)

  ## return if only one node in 'ids'
  if ( length(anc)==1 ) return(anc)

  ## get subtree for only "ids" nodes/leaves
  st <- getAllParents(ids, tax, rootid=rootid, dbug=dbug)
  root <- unique(st[!st[,1]%in%st[,2],1])
  if ( length(root)!=1 ) stop("there can only be one root")
  if ( sum(!st[,2]%in% st[,1])==1 ) # all ids in single lineage: get highest
    anc <- ids[ids %in% st[,1]][1]
  else { # several lineages
    ## order 'cladewise'
    st <- getChildren(root, st)
    ## get first branch in 'cladewise' ordered tree (i.e. skip 
    ## singelton chain to root)
    anc <- names(which(table(st[,1])>1)[1])
  }  
  anc
}

## loop over getParents for severals taxon IDs; avoids multiple
## lineage retrieval by passing existing edges to the recursive function
## call
## TODO: what is the order of edges returned!?
#' @export
getAllParents <- function(ids, tax, rootid="1", dbug=FALSE) {
  edges <- NULL
  for ( id in ids ) {
    if ( dbug )
      cat(paste("parents for id", id, ";",
                round(which(ids==id)/length(ids),2),"done\n"), file=stderr())
    edges <- rbind(edges,getParents(id,tax,skip=edges,rootid=rootid,dbug=dbug))
  }
  edges
}

## simple lineage retrieval, returns a single lineage vector
#' @export
get.parents <- function(id, tax, skip=NULL) {
  if ( id%in%tax[,2] & !id%in%skip ) id<-c(id,get.parents(tax[tax[,2]==id,1],tax,skip=skip))
  id
}

## lineage retrieval with 'existence check', returns edges!
## extract parents of one node from complete taxonomy tree
## 'skip' edges can be passed, and if the node 'id' is already present,
## the recursive function is not called
#' @export
getParents <- function(id, tax, skip=NULL, rootid="1", verb=TRUE, dbug=FALSE) {

    ## full taxonomy or just the edge table?
    ## TODO: change function to use full taxonomy list by parseNCBITaxonomy
    all <- NULL
    if ( typeof(tax)=="list" ) {
        all <- tax
        tax <- tax$taxon
    }
    
  childCol <- "child"
  if ( !childCol %in% colnames(tax) ) childCol <- 2
  parntCol <- "parent"
  if ( !parntCol %in% colnames(tax) ) parntCol <- 1
   
    if ( !id %in% tax[,childCol] ) 
        ## check "merged" in full taxonomy list
        if ( "merged"%in%names(all) )
            id <- updateIDs(id=id, tax=tax, verb=verb)
    
    ## does this ever happen for valid tax ids?
    ## (function should never be called for root node)
    if ( !id %in% tax[,childCol] ) {
        if ( id != rootid )
            warning(paste(id, "not found and skipped.\n"))
        else if ( dbug ) cat(paste("found root\n"))
        return(NULL)
    }
  
  ## is id already present?
  if ( id %in% c(skip) ) {
    if ( dbug ) cat(paste(id,"id present.\n"),file=stderr())
    return(NULL)
  }
  
  ## new edges, start from id as child
  nedges <- tax[tax[,childCol]==id,]
  parent <- nedges[parntCol]
  
  ## check if lineage is already present in 'skip' edges?
  linedone <- parent %in% skip[,childCol]
  if ( dbug & linedone ) cat(paste(parent,"lineage present.\n"),file=stderr())

  ## recursive parent retrieval until root is reached
  if ( parent != rootid & !linedone ) 
    nedges <- rbind(getParents(parent,tax,skip=skip, rootid=rootid,dbug=dbug),
                    nedges)
  nedges
}

## simple children retrieval: returns a list of all children
#' @export
get.children <- function(id, tax, skip=NULL) {
  if ( id %in% skip ) return(NULL)
  srtedg <- id
  for ( idx in which(tax[,1]==id) ) ## get indices of children
    srtedg <- c(srtedg, get.children(tax[idx,2],tax, skip=skip))
  srtedg
}

## get all chilrden of a node
## does a depth-first search thus returns edges ordered
## "cladewise" (in package ape)
#' @export
getChildren <- function(id, tax) {

  childCol <- "child"
  if ( !childCol %in% colnames(tax) ) childCol <- 2
  parntCol <- "parent"
  if ( !parntCol %in% colnames(tax) ) parntCol <- 1

  idcs <- which(tax[,parntCol]==id)
  srtedg <- NULL ## leave, if no idcs
  ## recursive call for all children
  for ( idx in idcs ) 
    srtedg <- rbind(srtedg,tax[idx,],getChildren(tax[idx,childCol],tax))
  srtedg
}

## load tbi taxonomy files
## TODO: children only have 1 parent -> convert edge list to hash
#' @export
parseNCBITaxonomy <- function(taxd, species=FALSE, ranks=FALSE, phylo=FALSE,
                              order=FALSE, merged=FALSE, deleted=FALSE,
                              verb=TRUE) {
  
    ## read NCBI taxonomy
    ## TODO: get IDs via scientific names instead!
    if ( verb )
        cat(paste("parsing NCBI taxonomy ",
                  " from '",taxd,"'\n",sep=""),file=stderr())
    
    taxf <- file.path(taxd, "nodes.dmp")
    tax <- read.table(taxf,header=FALSE, sep="|")
    edges <- cbind(parent=as.character(tax[,2]),child=as.character(tax[,1]))

    ## TODO - 201908: switch edges to hash
    if ( FALSE ) {
       parent <- as.character(tax[,2])
       names(parent) <-  tax[,1] 
    }
    
    ## get ranks
    rank <- NULL
    if ( ranks ) {
        rank <- gsub("\t","",as.character(tax[,3]))
        names(rank) <- edges[,"child"]
    }
    
    ## rm self-references, (in nodes.dmp root is parent of itself)
    edges <- edges[edges[,1]!=edges[,2],]
    
    ## read species names
    nam <- NULL
    if ( species ) {
        if ( verb )
            cat(paste("parsing scientific names, takes long\n"),file=stderr())
        namf <- gsub("nodes","names",taxf)
        ##fill=T because "line 1979/2105 did not have 5 elements" ??
        ## TODO: does this cause errors??
        nam <- read.table(namf,header=FALSE, sep="|", fill=TRUE,quote ="")
        nam <- nam[grep("scientific name",nam[,4]),]
        nam <- cbind(as.character(nam[,1]),
                     gsub("\t","",as.character(nam[,2])))

        ## replace name table by named vector
        nms <- nam[,2]
        names(nms) <- nam[,1]
        nam <- nms
    }

    ## parse merged/updated nodes
    mrg <- NULL
    if ( merged ) {
        if ( verb )
            cat(paste("parsing merged/obsolete IDs\n"),file=stderr())
        mrgf <- gsub("nodes","merged",taxf)
        mrg <- read.table(mrgf,header=FALSE, sep="|", fill=FALSE,quote ="")
        nms <- mrg[,2]
        names(nms) <- mrg[,1]
        mrg <- nms
    }

    ## TODO: parse deleted nodes
    
    ## TODO: convert to tree
    if ( phylo ) {
        if ( verb ) cat(paste("generating phylo object\n"),file=stderr())
        ## TODO: reorder rank if order=TRUE!
        return(list(phylo=tax2phylo(edges=edges, root="1", nodes=nam,
                                    order=order),
                    rank=rank, merged=mrg))
    } else return(list(taxon=edges, names=nam, rank=rank, merged=mrg))
    
}

## takes a list of edges with named nodes (e.g. NCBI taxonomy),
## and converts this to a ape::phylo (TODO: phylobase::phylo4) tree object,
## optionally renames edges if a node vector is given
tax2phylo <- function(edges, root="1", nodes=NULL, order=FALSE) {

  colnames(edges) <- c("parent","child")
  ## order edges 'cladewise' (ape, eq. to 'preorder' in phylobase)
  ## by depth-first search
  if ( order ) edges <- getChildren(root,tax=edges) 
  
  ## map edges to node and leave vectors
  inodes <- edges[which(edges[,"child"]  %in% edges[,"parent"]),"child"]
  leaves <- edges[which(!edges[,"child"]%in% edges[,"parent"]),"child"]
  allnodes <- 1:(length(inodes)+length(leaves)+1)
  ## order: leaves, root, internal nodes
  names(allnodes) <- c(leaves,root,inodes)
  edges <- cbind(parent=allnodes[edges[,"parent"]],
                 child=allnodes[edges[,"child"]])
  rownames(edges) <- names(inodes) <- names(leaves) <- NULL
  ## set node names
  lnames <- leaves
  names(lnames) <- leaves
  nnames <- inodes
  names(nnames) <- inodes
  rname <- root
  ## replace ID by node names
  if ( !is.null(nodes) ) {
    for ( inode in inodes )
      nnames[inode] <- nodes[nodes[,1]==inode,2]
    for ( leave in leaves )
      lnames[leave] <- nodes[nodes[,1]==leave,2]
    rname <- nodes[nodes[,1]==root,2]
  }


  ## generate phylo object (package ape)
  ## 1) full tree with internal node
  ## order: leaves; root, internal nodes
  tree <- list(edge=edges,
               Nnode=length(c(root,inodes)), tip.label=lnames,
               edge.length=rep(1,nrow(edges)), node.label=c(rname,nnames))
  class(tree) <- "phylo"
  if ( order ) 
    attributes(tree)$order <- "cladewise"
  tree
}

## get NCBI taxonomy from  taxon IDs as tree in newick format
#' @export
tax2newick <- function(ids, txd, species=FALSE, full=FALSE,ranks=FALSE,
                       verb=FALSE, dbug=FALSE) {

  if ( missing(ids) ) stop("No IDs provided!\n")
  if ( missing(txd) ) stop("No directory for NCBI taxonomy files provided!\n")

  ## read NCBI taxonomy
  txn <- parseNCBITaxonomy(txd, species=species, verb=verb,ranks=ranks)
  tax <- txn$taxon
  nam <- txn$names
  rank <- txn$rank
  
  if ( verb ) cat(paste("collecting lineages, may take long!",
                        length(ids), "taxa ..."),file=stderr())

  ## collect parent lineages for all IDs
  edges <- getAllParents(ids, tax, rootid="1",dbug=dbug)

  ## are there still duplicated?
  if ( sum(duplicated(edges[,"child"]))>0 )
    stop("duplicated edges found; should not happen!?")
  #edges <- edges[!duplicated(edges[,"child"]),]

  if ( verb )
    cat(paste(" ... done. Sorting",nrow(edges),"edges and",
              "creating phylo object ...\n"),file=stderr())

  ## find root - there can only be one, it should be "1" for NCBI!
  ## TODO: stop on or allow multiple?)
  root <- unique(edges[which(!edges[,"parent"] %in% edges[,"child"]),"parent"])
  ## convert named edges to phylo object!
  tree <- tax2phylo(edges,root=root,nodes=nam, order=TRUE) 

  ## collapse 'singleton nodes' with only one child
  if ( !full ) {
    if ( verb ) cat(paste("done. collapsing singletons!\n"),file=stderr())
    #require("ape")
    tree <- ape::collapse.singles(tree)
  }
  ## TODO: is this redundant by getChildren reordering?
  ## problem: reorder.phylo can't allocate too big vector!! 16.0 Gb WHY?
  #tree <- reorder.phylo(tree,"cladewise") 

  if ( verb ) cat(paste("done!\n"),file=stderr())
  
  list(tree=tree, rank=rank)
}


## NOTE: 201908 - using ranks and merged in full NCBI taxonomy list
## return requested taxonomy ranks, starting from taxid
getRank <- function(txid, tax, ranks=c("phylum","species"), names=FALSE) {

    ## make sure to handle these as characters
    txid <- as.character(txid)

    ## replace obsolete IDs
    txid[txid%in%names(tax$merged)] <-
        tax$merged[txid[txid%in%names(tax$merged)]]

    ## remove non-present IDs
    na <-!txid%in%tax$taxon[,2]
    if ( sum(na)>0 )
        warning(sum(na), " IDs not found in tree: ",
                paste(txid[na], collapse=";"))
    
    ## TODO: more efficient while loop that ends at highest rank!
    ## TODO: call getAllParents for multiple txids to reduce tree,
    ## then use this edge list below
    all <- matrix(NA, nrow=length(txid), ncol=length(ranks))
    all[!na,] <- t(sapply(txid[!na], function(id) {

        ## get complete lineage
        ## TODO: make this more efficient for similar taxa!
        lineage <- getParents(id, tax=tax)
        
        ## get ranks of lineage
        rnks <- tax$rank[c(lineage[,1],id)]
        
        ## look for requested ranks and return
        ares <- rep(NA, length(ranks))
        names(ares) <- ranks
        res <- unlist(sapply(ranks, function(x) names(rnks[rnks==x])))
        ridx <- names(res)
        if ( names )
            res <- tax[["names"]][res]
        ares[ridx] <- res
        ares
    }))
    #all <- t(all)
    colnames(all) <- ranks
    all
}

#' update a list of taxonomy IDs
#'@export
updateIDs <- function(ids, tax, verb=TRUE) {

    ## update all available
    idx <- which(ids%in%names(tax$merged))
    if ( verb & length(idx)>0 )
        cat(paste("updating", length(idx), "IDs\n"))
    if ( length(idx)>0 )
        ids[idx] <- tax$merged[ids[idx]]
    ids    
}
