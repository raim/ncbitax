

### TEST data.table::fread instead of read.table

## location of downloaded NCBI database
taxdb <- "~/work/CoilProject/halime/taxonomy/db"


namf <- file.path(taxdb, "names.dmp")
nam <- read.table(namf,header=FALSE, sep="|", fill=TRUE,
                  quote ="",comment.char="",
                  stringsAsFactors=TRUE)
nam <- nam[grep("scientific name",nam[,4]),]
nam <- cbind(as.character(nam[,1]),
                     gsub("\t","",as.character(nam[,2])))
nms <- nam[,2]
names(nms) <- nam[,1]
  

nam2 <- data.table::fread(namf,header=FALSE, sep="|", fill=TRUE,
                          quote ="",
                          data.table=FALSE, stringsAsFactors=TRUE)
nam2 <- nam2[grep("scientific name",nam2[,4]),]
nam2 <- cbind(as.character(nam2[,1]),gsub("\t","",as.character(nam2[,2])))
nms2 <- nam2[,2]
names(nms2) <- nam2[,1]

## check for same result
sum(nms!=nms2)
