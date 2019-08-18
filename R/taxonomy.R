#' ncbitax: a simple NCBI taxonomy parser
#'
#' Similar to functions provided by ape and phylobase,
#' but specific for NCBI Taxonomy. 
#'@author Rainer Machne \email{machne@hhu.de}
#'@docType package
#'@name ncbitax
#'@importFrom utils read.table installed.packages
NULL # this just ends the global package documentation

## PHYLOGENY/TAXONOMY: TREE PARSING AND GENERATION
## similar to functions provided by ape and phylobase,
## but specific for NCBI Taxonomy


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

edg2idx <- function(edges, parentCol=1, childCol=2) {
    
    if ( sum(duplicated(edges[,childCol]))>0 )
        stop("duplicated children (DAG) not allowed")

    idx <- edges[,parentCol]
    names(idx) <- edges[,childCol]
    idx
}
idx2edg <- function(idx) {
    edges <- matrix(NA, ncol=2, nrow=length(idx))
    cbind(parent=idx, child=names(idx))
}

#' simple lineage retrieval
#'
#' returns a simple lineage vector
#' @param id a taxon ID
#' @param tax NCBI taxonomy object
#' @param skip edges table to skip: if the node \code{id} is already present,
#' the recursive function is not called
#' @examples
#' get.parents("1148", tax)
#' get.parents("1140", tax)
#' @export
get.parents <- function(id, tax, skip=NULL) {
    if ( id%in%names(tax$parents) & !id%in%skip )
        id<-c(id,get.parents(tax$parents[names(tax$parents)==id],
                             tax,skip=skip))
  id
}

#' simple children retrieval
#'
#' returns a vector of all children
#' @param id a single taxon ID
#' @param tax NCBI taxonomy object
#' @param skip vectors of ids to "skip", if the node 'id' is already
#' present in \code{skip}, the recursive function is not called
#' @examples
#' get.children("1148", tax)
#' @export
get.children <- function(id, tax, skip=NULL) {
  if ( id %in% skip ) return(NULL)
  srtedg <- id
  for ( idx in which(tax$parents==id) ) ## get indices of children
    srtedg <- c(srtedg, get.children(names(tax$parents)[idx],tax, skip=skip))
  srtedg
}

#' get lowest common ancestor 
#'
#' retrieves the lowest common ancestor
#' @param ids vector of taxon IDs 
#' @param tax NCBI taxonomy object
#' @param root taxon ID of the taxonomic class to be considered as root
#' @param names return taxon names instead of ids
#' @param dbug print debug messages
#' @examples
#' getAncestor(ids=c("1148","1140"), tax)
#' @export
getAncestor <- function(ids, tax, root="1", names=FALSE, dbug=FALSE) {

    anc <- unique(ids)

    ## return if only one node in 'ids'
    if ( length(anc)==1 ) return(anc)

    ## get subtree for only "ids" nodes/leaves
    st <- getAllParents(ids, tax, root=root, dbug=dbug)
    
    root <- unique(st[!st[,1]%in%st[,2],1])
    if ( length(root)!=1 ) stop("there can only be one root")
    
    if ( sum(!st[,2]%in% st[,1])==1 ) # all ids in single lineage: get highest
        anc <- ids[ids %in% st[,1]][1]
    else { # several lineages
        ## order 'cladewise'
        tmp <- list()
        tmp$parents <- st[,1]
        names(tmp$parents) <- st[,2] 
        st <- getChildren(root, tmp)
        ## get first branch in 'cladewise' ordered tree (i.e. skip 
        ## singelton chain to root)
        anc <- names(which(table(st[,1])>1)[1])
    }
    if ( names ) getName(anc, tax)
    else anc
}

#' get common parents for a list of taxon IDs
#' 
#' Loops over getParents for severals taxon IDs and avoids multiple
#' lineage retrieval operations by passing already found edges to the
#' \code{skip} option of the recursive function.
## TODO: what is the order of edges returned!?
## TODO: work-horse recursion with storage of prior hits - they only
## speed increase in the whole package happens here -> BETTERER?
#' @param ids vector of taxon IDs 
#' @param tax NCBI taxonomy object
#' @param root taxon ID of the taxonomic class to be considered as root
#' @param names return taxon names instead of ids
#' @param dbug print debug messages
#' @examples
#' getAllParents(c("1148","1140"), tax)
#' @export
getAllParents <- function(ids, tax, root="1", names=FALSE, dbug=FALSE) {

    edges <- NULL
    for ( id in ids ) {
        if ( dbug )
            cat(paste("parents for id", id, ";",
                      round(which(ids==id)/length(ids),2),"done\n"),
                file=stderr())
        edges <- rbind(edges,
                       getParents(id,tax,skip=edges,root=root,names=FALSE,
                                  dbug=dbug))
    }
    if ( names ) getName(edges, tax)
    else edges
}


#' lineage retrieval 
#'
#' Returns an edge table of the parent lineage until \code{root}.
#' @param id a single taxon ID
#' @param tax NCBI taxonomy object
#' @param root taxon ID of the taxonomic class to be considered as root
#' @param skip edges table to skip, if the node \code{id} is already present,
#' the recursive function is not called (used internally
#' from \code{\link{getAllParents}})
#' @param names return taxon names instead of ids
#' @param verb print progress messages
#' @param dbug print debug messages
#' @examples
#' getParents("1140", tax)
#' getParents("1148", tax, root="1783272")
#' getParents("1048", tax, root="1783272") # root not in lineage!
#' @export
getParents <- function(id, tax, skip=NULL, root="1", names=FALSE,
                       verb=TRUE, dbug=FALSE) {

    ## this occurs e.g. when specifying a "root" that
    ## does not appear in the lineage of "id"
    ## TODO: avoid, or provide better warning?
    if ( !id %in% names(tax$parents) ) {
        if ( id != root )
            warning(paste(id, "not found and skipped.\n"))
        else if ( dbug ) cat(paste("found root\n"))
        return(NULL)
    }
    
    ## is id already present?
    if ( id %in% c(skip) ) { #c(skip,names(skip)) ) {
        if ( dbug ) cat(paste(id,"id present.\n"),file=stderr())
        return(NULL)
    }
  
    ## new edges, start from id as child
    parent <- tax$parents[names(tax$parents)==id] # tax[tax[,childCol]==id,]
    
    ## edge list
    nedges <- cbind(parent, names(parent))
    
    ## check if lineage is already present in 'skip' edges?
    linedone <- parent %in% skip[,2]
    if ( dbug & linedone ) cat(paste(parent,"lineage present.\n"),file=stderr())
    
    ## recursive parent retrieval until root is reached
    if ( parent != root & !linedone ) 
        nedges <- rbind(getParents(parent, tax,
                                   skip=skip, root=root,dbug=dbug),
                        nedges)
    colnames(nedges) <- c("parent","child")

    if ( names ) getName(nedges, tax)
    else nedges
}



#' get all chilrden of a node
#'
#' Returns edge table of all descendents. Performs
#' a depth-first search thus returns ordered edges
#' (order "cladewise" in package ape).
#' @param id a single taxon ID
#' @param tax NCBI taxonomy object
#' @param names return taxon names instead of ids
#' @export
getChildren <- function(id, tax, names=FALSE) {

    idcs <- which(tax$parent==id)
    srtedg <- NULL ## leave, if no idcs
    ## recursive call for all children
    for ( idx in idcs ) 
        srtedg <- rbind(srtedg,
                        c(tax$parent[idx],names(tax$parent[idx])),
                        getChildren(names(tax$parent[idx]),tax))
    if ( names ) getName(srtedg, tax)
    else srtedg
}

#' find taxon ID for species names
#'
#' Simply uses base R's \code{\link{grep}} to grep for a
#' name pattern in the official NCBI taxon names and returns
#' all matching taxon IDs as a table with columns ID and name.
#' @param pattern regular expression, see argument \code{pattern}
#' in \code{\link{grep}}
#' @param tax NCBI taxonomy object
#' @param ... arguments to \code{\link{grep}}
#' @examples
#' grepName("Synechocystis", tax)
#'@export
grepName <- function(pattern, tax, ...) {
    idx <- grep(pattern, tax$names, ...)
    cbind(ID=names(tax$names)[idx], NAME=tax$names[idx])
}

#' gets taxon names from taxon IDs
#' @param ids vector or edge table of taxon IDs 
#' @param tax NCBI taxonomy object
#' @examples
#' getName(c("1148","1140"), tax)
#' @export
getName <- function(ids, tax) {
    if ( class(ids)=="matrix" )
        apply(ids, 2, function(x) tax[["names"]][x])
    else
        tax[["names"]][ids]
}

#' gets taxon ranks from taxon IDs
#' @param ids vector or edge table of taxon IDs 
#' @param tax NCBI taxonomy object
#' @examples
#' getRank(c("1148","1140"), tax)
#' @export
getRank <- function(ids, tax) {
    if ( class(ids)=="matrix" )
        apply(ids, 2, function(x) tax[["rank"]][x])
    else
        tax[["rank"]][ids]
}

#' load tbi taxonomy files
#' @param taxd path to unpacked NCBI taxonomy (taxdump.tar.gz)
#' @param names load species names, takes long to load!
#' @param ranks load rank names
#' @param merged load list of merged (outdated) taxon IDs and their replacement
#' @param verb print progress messages
#' @examples
#' taxdb <- "/data/taxonomy/db"
#' tax <- parseNCBITaxonomy(taxdb)
#' @export
parseNCBITaxonomy <- function(taxd, names=TRUE, ranks=TRUE, merged=TRUE,
                              verb=TRUE) {
  
    ## read NCBI taxonomy
    if ( verb )
        cat(paste("parsing NCBI taxonomy ",
                  " from '",taxd,"'\n",sep=""),file=stderr())
    
    ## check if read.table is installed, warn otherwise
    fread.available <- is.element("data.table", installed.packages()[,1])
    if ( !fread.available )
        warning("parsing NCBI taxonomy is much faster if package",
                " read.table is installed (fread)")
    
    taxf <- file.path(taxd, "nodes.dmp")
    if ( fread.available )
        tax <- data.table::fread(taxf,header=FALSE, sep="|",data.table=FALSE)
    else
        tax <- read.table(taxf,header=FALSE, sep="|", comment.char="")
    edges <- cbind(parent=as.character(tax[,2]),child=as.character(tax[,1]))

    ## MAIN EDGE LIST as hash `parent[child]`
    ## TODO: use edg2idx
    parent <- as.character(tax[,2])
    names(parent) <-  tax[,1] 
    ## rm self-references, (in nodes.dmp root is parent of itself)
    parent <- parent[parent!=names(parent)]
    
    ## get ranks
    rank <- NULL
    if ( ranks ) {
        rank <- gsub("\t","",as.character(tax[,3]))
        names(rank) <- tax[,1] 
    }
    

    ## parse merged/updated nodes
    mrg <- NULL
    if ( merged ) {
        if ( verb )
            cat(paste("parsing merged/obsolete IDs\n"),file=stderr())
        mrgf <- gsub("nodes","merged",taxf)
        if ( fread.available )
            mrg <- data.table::fread(mrgf,header=FALSE, sep="|",
                                     fill=FALSE,quote ="", data.table=FALSE)
        else mrg <- read.table(mrgf,header=FALSE, sep="|",
                               fill=FALSE,quote ="",comment.char="")
        nms <- mrg[,2]
        names(nms) <- mrg[,1]
        mrg <- nms
    }

    ## TODO: parse deleted nodes

    ## parse species names
    ## this takes long, and is only required to access
    ## names for taxon IDs
    nam <- NULL
    if ( names ) {
        if ( verb )
            cat(paste("parsing scientific names, takes long\n"),file=stderr())
        namf <- gsub("nodes","names",taxf)
        ##fill=T because "line 1979/2105 did not have 5 elements" ??
        ## TODO: does this cause errors??
        if ( fread.available )
            nam <- data.table::fread(namf,header=FALSE, sep="|", fill=TRUE,
                                     quote ="", data.table=FALSE)
        else
            nam <- read.table(namf,header=FALSE, sep="|", fill=TRUE,quote ="",
                              comment.char="")
        nam <- nam[grep("scientific name",nam[,4]),]
        nam <- cbind(as.character(nam[,1]),
                     gsub("\t","",as.character(nam[,2])))

        ## replace name table by named vector
        nms <- nam[,2]
        names(nms) <- nam[,1]
        nam <- nms
    }
    
    ## TODO: order for more efficient searches?
    ## TODO: check class NCBItaxonomy in other functions, or
    ## use as class method for parent/child/ancestor?
    
    tax <- list(parents=parent, names=nam, rank=rank, merged=mrg)
    class(tax) <- "NCBItaxonomy"
    tax
    
}

## TODO: phylobase::phylo4 tree
#' converts an NCBI taxonomy object to an \code{ape::phylo} tree
#'
#' Only a sub-tree until \code{root] is returned. Pptionally renames
#' edges to NCBI taxon names
#' @param tax NCBI taxonomy object
#' @param root taxon ID of the taxonomic class to be considered as root
#' @param names replace taxon IDs by names
#' @param order order edges, required for attribute "cladewise" of
#' the returned \code{ape::phylo} tree
#' @examples
#' phy <- tax2phylo(tax, root="1890428", order=TRUE, names=TRUE)
#' ape::plot.phylo(phy)
#' @export
tax2phylo <- function(tax, root="1", names=FALSE, order=FALSE) {


    ## get edges table
    ## order edges 'cladewise' (ape, eq. to 'preorder' in phylobase)
    ## by depth-first search
    if ( order ) edges <- getChildren(root,tax=tax)
    else edges <- cbind(tax$parents, names(tax$parents))
    colnames(edges) <- c("parent","child")
  
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
    if ( names ) {
        nnames[inodes] <- tax$names[as.character(inodes)]
        lnames[leaves] <- tax$names[as.character(leaves)]
        rname <-  tax$names[as.character(root)]
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

#' convert NCBI taxonomy sub trees to newick format
#' @param ids vector of taxon IDs 
#' @param tax NCBI taxonomy object
#' @param names replace taxon IDs with names
#' @param full if  \code{TRUE} \code{singleton} nodes with only
#' one child will NOT be collapsed
#' @param ranks add rank information for nodes
#' @param verb print progress messages
#' @param dbug print debug messages
#' @examples
#' new <- tax2newick(ids=c("1148","1140"), tax=tax, names=TRUE)
#' @export
tax2newick <- function(ids, tax, names=FALSE, full=FALSE,ranks=FALSE,
                       verb=FALSE, dbug=FALSE) {

    if ( missing(ids) ) stop("No IDs provided!\n")
    if ( missing(tax) ) stop("No NCBI taxonomy provided!\n")

    ## read NCBI taxonomy
    if ( typeof(tax)=="character" )
        tax <- parseNCBITaxonomy(tax, names=names, verb=verb,ranks=ranks)

    
    if ( verb ) cat(paste("collecting lineages, may take long!",
                          length(ids), "taxa ..."),file=stderr())
    
    ## collect parent lineages for all IDs
    edges <- getAllParents(ids, tax, root="1",dbug=dbug)
    
    ## are there still duplicated?
    if ( sum(duplicated(edges[,"child"]))>0 )
        stop("duplicated edges found; should not happen!?")
    ##edges <- edges[!duplicated(edges[,"child"]),]
    
    if ( verb )
        cat(paste(" ... done. Sorting",nrow(edges),"edges and",
                  "creating phylo object ...\n"),file=stderr())
    
    ## find root - there can only be one, it should be "1" for NCBI!
    ## TODO: stop or allow multiple?
    root <- unique(edges[which(!edges[,"parent"] %in%
                               edges[,"child"]),"parent"])

    ## convert named edges to phylo object!
    pars <- edges[,1]
    names(pars) <- edges[,2]
    tmp <- list()
    tmp$parents <- pars
    tmp$names <- tax$names
    tree <- tax2phylo(tax=tmp, root=root, names=names, order=TRUE) 
    
    ## collapse 'singleton nodes' with only one child
    if ( !full ) {
        if ( verb ) cat(paste("done. collapsing singletons!\n"),file=stderr())
        tree <- ape::collapse.singles(tree)
    }
    ## TODO: is this redundant by getChildren reordering?
    ## problem: reorder.phylo can't allocate too big vector!! 16.0 Gb WHY?
    ##tree <- reorder.phylo(tree,"cladewise") 
    
    if ( verb ) cat(paste("done!\n"),file=stderr())
    
    list(tree=tree, rank=tax$rank)
}


#' get NCBI taxonomy rank lineages
#'
#' Returns a table with taxonomic lineages at the requested
#' rank levels. A wrapper around \code{\link{getAllParents}}.
#' @param ids vector of taxon IDs 
#' @param tax NCBI taxonomy object
#' @param ranks vector of ranks to retrieve
#' @param names return taxon names instead of ids
#' @param reduce call \code{\link{reduceTaxonomy}} before search; untested
#' but this makes subsequent search faster by using \code{\link{getAllParents}}
#' first to reduce taxonomy to sub-tree for \code{ids}
#' @examples
#' getLineage(c("1148","1140"), tax, names=TRUE,
#'          ranks=c("superkingdom","phylum","species"))
#' @export
getLineage <- function(ids, tax, ranks=c("phylum","species"),
                       names=FALSE, reduce=FALSE) {
    
    ## make sure to handle these as characters
    ids <- as.character(ids)

    ## remove non-present IDs
    na <-!ids%in%names(tax$parents)
    if ( sum(na)>0 )
        warning(sum(na), " IDs not found in tree: ",
                paste(ids[na], collapse=";"))

    if ( reduce )
        tax <- reduceTaxonomy(ids, tax)
    
    ## TODO: more efficient while loop that ends at highest rank!
    ## TODO: call getAllParents for multiple txids to reduce tree,
    ## then use this edge list below
    all <- matrix(NA, nrow=length(ids), ncol=length(ranks))
    all[!na,] <- t(sapply(ids[!na], function(id) {
        
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
    colnames(all) <- ranks
    all
}


#' reduced taxonomy tree for faster searches
#'
#' This reduces the full NCBI taxonomy to the subset
#' containing all taxons in \code{ids}. This makes all
#' subsequent searches much faster.
#' @param ids vector of taxon IDs 
#' @param tax NCBI taxonomy object
#' @examples
#' rtax <- reduceTaxonomy(c("1148","1140"), tax)
#' lapply(tax, length)
#' lapply(rtax, length)
#' @export
reduceTaxonomy <- function(ids, tax) {

    ## exception: call updateIDs here, since it would
    ## hamper all subsequent searches to loose outdated/merged IDs here
    ids <- updateIDs(ids, tax)

    ## get all parents up to root
    pars <- getAllParents(ids, tax, root="1")

    ## generate reduced taxonomy object
    rtax <- tax

    ## edge list: named vector parent[child]
    rtax$parents <- pars[,1]
    names(rtax$parents) <- pars[,2]
    ## scientific names
    if ( "names" %in% names(tax) )
        rtax$names <- rtax$names[c(pars)] #[c(pars)%in%names(rtax$names)]
    if ( "rank" %in% names(tax) )
        rtax$rank <- rtax$rank[c(pars)] #[c(pars)%in%names(rtax$names)]

    ## but keep full merged list!
    
    class(rtax) <- "NCBItaxonomy"
    rtax
}

#' update a list of taxonomy IDs
#'
#' Looks in merged node list from file \code{merged.dmp}, whether
#' the passed \code{ids} have been replaced (merged) and returns
#' a list of updated (where in merged) ids.
#' It is recommended to run this on all IDs before using other
#' functions of this package!
#' @param ids vector of taxon IDs
#' @param tax NCBI taxonomy object
#' @param verb print progress messages
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
