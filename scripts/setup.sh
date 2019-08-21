
## DOWNLOAD AND SETUP TAXONOMY AND PHYLOGENETIC INFO

## TODO: convert this to Rmarkdown vignette, move
## required sripts to R package compatible locations

SRCPATH=~/programs/ncbitax/scripts
TAXPATH=/data/taxonomy/

### NCBI TAXONOMY
## download taxonomy DB and collect taxonomy information for all
## blasted genomes
mkdir -p $TAXPATH/ncbi
cd $TAXPATH/ncbi
date > $TAXPATH/ncbi/downloaddate.txt #comes w/o version, record date

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxvf taxdump.tar.gz
## NCBI: taxid refseq ID mapping
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
gunzip gi_taxid_nucl.dmp.gz
## NBCI: genbank to nucleotide accession (gi) mapping 
wget ftp://ftp.ncbi.nih.gov/genbank/livelists/GbAccList.0407.2019.gz
#gunzip GbAccList.0407.2019.gz # <- HUGE FILE

### PROKARYOTES - 16S rRNA tree

## download 16S rRNA tree from silva and extract IDs
mkdir -p $TAXPATH/silva
cd $TAXPATH/silva
## SILVA: species AccNum mapping
wget https://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_132/LTPs132_SSU.csv
## SILVA: 16S rRNA tree
wget https://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_132/LTPs132_SSU_tree.newick

## map accession numbers from SILVA 16S Tree to taxon ids
## 1) find accession number in GbAccList from NCBI and get GI identifier
gunzip  $TAXPATH/ncbi/GbAccList.0407.2019.gz
cut -f 1 $TAXPATH/silva/LTPs132_SSU.csv  | perl $SRCPATH/retrieveGIDs.pl -f $TAXPATH/ncbi/GbAccList.0407.2019 -d "," -c 0 > $TAXPATH/silva/LTPs132_SSU_gb2gi.csv 2> $TAXPATH/silva/LTPs132_SSU_gb2gi.log
rm  -f $TAXPATH/ncbi/GbAccList.0407.2019
## 2) find GI in gi_taxid_nucl.dmp
gunzip $TAXPATH/ncbi/gi_taxid_nucl.dmp
cut -d , -f 3 $TAXPATH/silva/LTPs132_SSU_gb2gi.csv | perl $SRCPATH/retrieveGIDs.pl -f $TAXPATH/ncbi/gi_taxid_nucl.dmp -d "\t" -c 0  > $TAXPATH/silva/LTPs132_SSU_gi2tax.csv 2> $TAXPATH/silva/LTPs132_SSU_gi2tax.log
rm -f $TAXPATH/ncbi/gi_taxid_nucl.dmp

## NOTE: see LTPs132_SSU_gb2gi.log
## for species with multiple 16S rRNA sequences !

## generate tree with taxon IDs: the script generates two files;
## taxonomy/LTPs132_SSU_tree2tax.dat  maps the old species names to taxon IDs
## taxonomy/LTPs132_SSU.tree_taxonIDs.newick is the 16S rRNA tree with 
##    taxon IDs and species names (but maintaining clade names)

## 1) replace spaces in tree file
sed s'/ \+/XYZYX/g' $TAXPATH/silva/LTPs132_SSU_tree.newick > $TAXPATH/silva/LTPs132_SSU_tree_nospaces.newick ## NOTE: XYZYX is replaced with spaces in the convert16Stree.R
## replace path in convert16Stree.R
cp  -a $SRCPATH/convert16Stree.R  $TAXPATH/silva/tmp.R
tmppath=$TAXPATH/silva
sed -i "s|tpath\=.*|tpath\='$tmppath'|" $TAXPATH/silva/tmp.R
$TAXPATH/silva/tmp.R # TODO: use command-line option for version

### NOTE: end general taxonomy update script here
###       start taxon ID input here


### CYANOBACTERIA protein consensus tree

## download phylogenetic tree for cyanobacteria from Shih et al. 2013
mkdir -p $TAXPATH/cyanobacteria
cd $TAXPATH/cyanobacteria
## Supp. table S1 with species meta information 
wget https://www.pnas.org/highwire/filestream/611359/field_highwire_adjunct_files/1/sd01.xls
## convert to csv
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator=; eol=unix sheet='cyano_metdata'" $TAXPATH/cyanobacteria/sd01.xls $TAXPATH/cyanobacteria/shih13_metadata.csv

## Shih et al. 2013 - protein consensus tree
## sent by Patrick Shih <pmshih@berkeley.edu> on 20140401, contains
## 4 outgroups
cp -a /data/lehmann13/data/genomes/cyanoGEBAspeciesTree.txt $TAXPATH/cyanobacteria/

## GOLD ID vs. NCBI tax mapping
mkdir $TAXPATH/gold
## manually from https://gold.jgi.doe.gov/download?mode=site_excel - save as
##  $TAXPATH/gold/goldData.xlsx
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "locale=C format=automatic separator=; eol=unix sheet='site data'" $TAXPATH/gold/goldData.xlsx $TAXPATH/gold/gold.csv

## for each gold ID in shih13_metadata.csv data get the ncbi taxonomy
## info via gold.csv, update taxonomy ID and add the name of the species
## in the tree file
R --vanilla < $SRCPATH/mapShih13.R 

