
library(ape)
library(ncbitax)
library(segmenTools) # just for plotdev

tax.path <- "/data/taxonomy"

## read meta info from shih et al. 2013 genomes
shih <- read.csv(file.path(tax.path, "cyanobacteria","shih13_metadata.csv"),
                 sep=";", stringsAsFactors=FALSE,skip=2)
shih <- shih[-1,] # additional column info
shih <- shih[shih[,1]!="",]

shih[,"Organism"] <- sub(" +\\(T\\)","",shih[,"Organism"])
shih[,"Organism"] <- trimws(shih[,"Organism"])
gid <- shih[,"Gold.card.ID"]

## read GOLD database info
gold <- read.csv(file.path(tax.path,"gold","gold.csv"),sep=";",
                 stringsAsFactors=FALSE)

idx1 <- match(shih[,"Gold.card.ID" ],
              gold[,"LEGACY.GOLD.ID"])

idx2 <- match(paste0("PRJNA",shih[,"NCBI.proj.ID"]),
              gold[,"NCBI.BIOPROJECT.ACCESSION"])

idx3 <- match(sub(" \\(T\\)", "",shih[,"Organism"]),
              gold[,"PROJECT.NAME"])

## check and merge indices
idx <- cbind(gold=idx1,ncbi=idx2,name=idx3,merged=idx1)

## 3 conflicts, manual inspection shows that GOLD ID fits
## and NCBI Projekt number is wrong
shih[which(idx1!=idx2),]

## 
cat(paste("found", sum(!is.na(idx[,"merged"])), "via GOLD ID\n"))

## take ncbi match, where no gold ID match is available
cat(paste("searching", sum(is.na(idx[,"merged"])), "via NCBI Project number\n"))
idx[is.na(idx[,"merged"]),"merged"] <- idx[is.na(idx[,"merged"]),"ncbi"] 
## take name match, where neither gold ID nor ncbi project is available
## -> only one, manually checked, Cyanobacterium aponinum PCC 10605
cat(paste("searching", sum(is.na(idx[,"merged"])), "via organism name\n"))
idx[is.na(idx[,"merged"]),"merged"] <- idx[is.na(idx[,"merged"]),"name"] 

## final mapping
idx <- idx[,"merged"]
## adding NCBI taxonomy ID 
shih <- cbind(shih, taxon=gold[idx,"NCBI.TAXON.ID"])

## map taxon IDs and names in tree
tree <- read.tree(file.path(tax.path, "cyanobacteria",
                            "cyanoGEBAspeciesTree.txt"))

## auto-search
shiNames <- gsub("PCC", "PCC ", tree$tip.label)
shi2spec <- rep(NA, nrow(shih))

for ( i in 1:length(shiNames) ) {
  idx <- grep(shiNames[i],shih[,"Organism"])
  if ( length(idx)>0) {
    #cat(paste(shiNames[i],"\t",shih[idx,"Organism"],"\n"))
    shi2spec[idx] <- i
    next
  }
  ## split first 3 and last 4
  str <- unlist(strsplit(shiNames[i],""))
  idx3 <- grep(paste(c("^",str[1:3]),collapse=""), shih[,"Organism"])
  idx4 <- grep(paste(c(str[4:7],"$"),collapse=""), shih[,"Organism"])
  ## match by last four numbers
  if ( length(idx4) ) {
      cat(paste(shiNames[i],"\t",shih[idx4,"Organism"],"\n"))
      shi2spec[idx4] <- i
    }
}

map <- cbind(shih[,"Organism"], tree$tip.label[shi2spec], rep(NA,nrow(shih)))

## manually add missing
miss <- c(AmaMBIC="Acaryochloris marina MBIC11017",
          Cwat003="Crocosphaera watsonii WH 0003",          
          Maer843="Microcystis aeruginosa NIES-843",        
          ProNAT1="Prochlorococcus marinus NATL1A",                  
          ProNAT2="Prochlorococcus marinus NATL2A",                  
          ProMED4="Prochlorococcus marinus subsp. pastoris CCMP1986",
          PdideP1="Prochloron didemni P1",                           
          SynB107="Synechococcus sp. BL107",                         
          SynJ23B="Synechococcus sp. JA-2-3B",                       
          SynJ33A="Synechococcus sp. JA-3-3Ab",                      
          SynR307="Synechococcus sp. RCC307",                        
          TeloBP1="Thermosynechococcus elongatus BP-1",              
          UNFixCy="Unicellular cyanobacterium UCYN-A",               
          Amax328="Arthrospira maxima CS-328",                       
          ArthPla="Arthrospira platensis NIES-39",                   
          Artpara="Arthrospira platensis Paraca",                    
          Mchthon="Coleofasciculus chthonoplastes PCC 7420",  #Microcoleus chthonoplastes PCC 7420      
          Lyn8106="Lyngbya sp. CCY 9116", # via gold table info
          MvagFGP="Microcoleus vaginatus FGP-2",                     
          "Moorea producta 3L", # "Moorea producens 3L" in gold
          Tery101="Trichodesmium erythraeum IMS101",                 
          Crac505="Cylindrospermopsis raciborskii CS-505",           
          RaphiBr="Raphidiopsis brookii D9",
          PCC9605="Fischerella sp. PCC 9605",
          PCC9431="Fischerella sp. PCC 9431",
          PCC9339="Fischerella sp. PCC 9339")

## add missing manually
idx <- match(map[,1], miss)


map[!is.na(idx),2] <- names(miss)[idx[!is.na(idx)]]

map[map[,2]=="",]
tree$tip.label[!tree$tip.label%in%map[,2]]

rownames(map) <- map[,2]

## replace names
shi <- tree
shi$tip.label[tree$tip.label%in%rownames(map)] <-
    map[tree$tip.label[tree$tip.label%in%rownames(map)],1]

## drop those not found
shi <- drop.tip(shi,tree$tip.label[!tree$tip.label%in%rownames(map)])

## set root at Gloeobacter violaceus
shi <- root(shi,grep("7421",shi$tip.label),resolve.root=T)
#shi <- root(shi,grep("Trichodesmium",shi$tip.label),resolve.root=T)

plotdev(file.path(tax.path,"cyanobacteria","shih13tree"),
        type="pdf", width=5, height=7)
par(mai=c(.1,.1,.1,.1))
plot(shi, show.node.label=TRUE,edge.width=2,no.margin=TRUE, cex=.7)
     #tip.color=1:length(shi$tip.label))
add.scale.bar()
dev.off()

## TODO: switch location of Oscillatoria to lie next to chamaesiphon

## update taxonomy IDs
library(ncbitax)
tax <- parseNCBITaxonomy(file.path(tax.path,"ncbi"))

## NOTE: currently all IDs available, but perhaps changes with later
## taxonomies
shih$taxon <- updateIDs(shih$taxon, tax)

## bind short names
shih <- cbind(shih, tree.name=map[,2])
write.table(shih,file.path(tax.path, "cyanobacteria",
                           "shih13_metadata_annotated.tsv"),
            sep="\t", quote=FALSE, na="")
## write tree with 
write.tree(shi,file=file.path(tax.path, "cyanobacteria",
                              "cyanoGEBAspeciesTree_names.txt"))
