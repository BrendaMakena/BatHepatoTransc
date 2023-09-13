library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")


## where did this table acually come from MAKE A NOTE, so that we are
## able to cite this!!!
proteins <- read.csv("inputdata/rousettus_proteins.csv")

getBM(mart=mart,
      attributes = c("entrezgene_accession", "entrezgene_id",
                     "go_id"),
      filters =  "entrezgene_accession",                       
      values = "PPP6R3") ### value in "Locus" testing one

### This works only for the short abreviated "gene names", not for the
### accessions starting with LOCxxxxxxx (x are numbers)


### But let's see whether this is a problem and how big...

### Before you use this you might want to subset your "proteins" to
### only those really tested in the DE analysis!

allLocusGO <- getBM(mart=mart,
                    attributes = c("entrezgene_accession",
                                   "entrezgene_id",
                                   "entrezgene_description",
                                   "go_id"),
                    filters =  "entrezgene_accession",
                    values = proteins$Locus) ###all values in "Locus"

table(unique(proteins$Locus)%in%allLocusGO$entrezgene_accession)
## FALSE  TRUE 
##  3540 15529 

## we miss (only?) 3540 genes in the overall set
## how is this for the expression-analysed genes?

table(found=unique(proteins$Locus)%in%allLocusGO$entrezgene_accession,
      LOC=grepl("LOC", unique(proteins$Locus)))

## even better: almost all the ones we're are really those with "LOC"

## how is this for the expression-analysed genes?
