
## TESTING speed improvement by reducing taxonomy before multiple searches

#source("~/programs/ncbitax/R/taxonomy.R")
library(ncbitax)

## example IDs for testing
## TODO: use full list for sampling
PATH <- "/data/topoisomerases"
examples <- file.path(PATH,"all_taxonomy_IDs.dat")


## location of downloaded NCBI database
taxdb <- "~/work/CoilProject/halime/taxonomy/db"
## LOAD NCBI taxonomy
tax <- parseNCBITaxonomy(taxdb, names=TRUE, ranks=TRUE, merged=TRUE)

## test 10 and perm 100 takes 10-15 min
test <- 50 #FALSE #   ## number of random taxon IDs
perm <- 10 # number of permutations
random.species <- FALSE # TRUE # 

## load taxonomy IDs from ALL species with full genomes in blast DB
if ( random.species ) {
    ## ~RANDOM IDs: about 30% faster
    examplids <-  as.character(read.csv(examples)[,1])
} else {
    ## random children of a specific taxon : about 15% faster
    examplids <- getChildren(getAncestor(ids=c("1148","1140"), tax), tax)[,2]
}


tot <- matrix(NA, ncol=2, nrow=perm)
for ( r in 1:perm ) {
    
    origids <-examplids
    
    if ( test>0 & test<Inf )
        origids <- sample(origids, test)
    
    ## update IDs
    taxids <- updateIDs(origids, tax)
    

    ## reduce taxonomy to sub-tree containing all taxids
    ## THIS CAN TAKE LONG, but getRank may then be faster - TODO: test speed!
    start.red <- Sys.time()
    rtax <- reduceTaxonomy(taxids, tax)
    end.red <- Sys.time()
    phy <- getRank(ids=taxids, tax=rtax,
                   ranks=c("superkingdom","kingdom",
                           "phylum","family","species"),names=TRUE,
                   reduce=FALSE)
    end.rank <- Sys.time()

    ## search on full taxonomy
    phy2 <- getRank(ids=taxids, tax=tax,
                    ranks=c("superkingdom","kingdom",
                            "phylum","family","species"),names=TRUE,
                   reduce=FALSE)
    end.full <- Sys.time()
    tot[r,] <- c(reduced=difftime(end.rank, start.red, units="secs"),
                 full=difftime(end.full, end.rank, units="secs"))
                 
    cat(paste("DONE, timing:\n",
              "reducing taxonomy\t",
              round(difftime(end.red,start.red,units="secs"),3), "sec\n",
              "collecting ranks \t",
              round(difftime(end.rank,end.red,units="secs"),3), "sec\n",
              "search on full   \t",
              round(difftime(end.full,end.rank,units="secs"),3), "sec\n"))
    
    unlist(lapply(tax, length))
    unlist(lapply(rtax, length))
}
    
plot(tot)
dff <- apply(tot,1, diff)
hist(tot[,1]/tot[,2],breaks=50)
