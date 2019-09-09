
source("~/programs/ncbitax/R/taxonomy.R")
tax <- parseNCBITaxonomy("/data/taxonomy/ncbi")

getID(c("Actinobacteria","Cyanobacteria"), tax, all=FALSE)
getID(c("Actinobacteria","Cyanobacteria"), tax, all=TRUE)
getID(c("Actinobacteria","Cyanobactera"), tax, all=TRUE)
getID(c("Actinobacteria","Cyanobactera"), tax, all=FALSE)
