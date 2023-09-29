## annotation of transcripts:

## Hepatocystis transcripts are queried here and annotations are
## printed to the console

## Bat GO term transcript annotations are obtained from the ENSEMBLE
## biomart and stored in an intermediated file "intermediateData/GOtermAnnot.RDS"


library(rtracklayer)
library(biomaRt)

## rerun script 1 for generating counts table or read from
## intermediate data
readcount <- FALSE

## loading the dataframe (file containing the read counts) from the features count step
if(readcount){
  source("scripts/3_featurecounts.R")
}else{
  tagseqRNAfeatureCounts <- readRDS("intermediateData/countTable.RDS")
}

## loading the proteins annotations table from ncbi
rousettus_proteins <- read.csv("inputdata/rousettus_proteins.csv")

## the merged gtf file contains all the information for hepatocystis
merged_gtf <- readGFF("/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/RousettusHepatocystis_merged.gtf")

## for hepatocystis annotation this is enought, all information is in the gff file
merged_gtf <- merged_gtf[merged_gtf$gene_id %in% rownames(tagseqRNAfeatureCounts),]


# adding species column to the gtf file
merged_gtf$species <- ifelse(grepl("HEP_", merged_gtf$gene_id),
                             "hepatocystis",
                             "rousettus")

# mitochondrial genome not mapped to by our transcripts  
rowSums(tagseqRNAfeatureCounts %>% 
  filter(grepl("HEP_MIT",
               rownames(tagseqRNAfeatureCounts),
               ignore.case = TRUE)))

## numeric(0)
sum(grepl("HEP_MIT", rownames(tagseqRNAfeatureCounts), ignore.case = TRUE))
## [1] 0


# what the highest expressed hepatocystis transcript is
merged_gtf[merged_gtf$gene_id %in% "HEP_00048200",]   #its a 28S rRNA

# what our hepatocystis transcripts are
merged_gtf %>%
        subset(gene_id %in%
(rownames(tagseqRNAfeatureCounts %>% 
  filter(grepl("HEP_",
               rownames(tagseqRNAfeatureCounts),
               ignore.case = TRUE),
         rowSums(tagseqRNAfeatureCounts)>5))),
select = c(gene_id,gene_biotype,product)) %>%
  as.data.frame()
## they are all rRNA


# what our rousettus transcripts are
merged_gtf %>%
  subset(gene_id %in%
  (rownames(tagseqRNAfeatureCounts %>% 
  filter(!grepl("HEP_",
  rownames(tagseqRNAfeatureCounts),
  ignore.case = TRUE),
  rowSums(tagseqRNAfeatureCounts)>5))),
  select = c(gene_id,gene_biotype,product)) %>%
  as.data.frame()
## proper protein coding genes

# selecting columns of interest in the merged gtf file
merged_gtf <- merged_gtf[c("gene_id")]
                             

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")

### query a human biomart for the gene symbols available in the
### "locus" column of the protein table

### Those are usually human (universal) gene symbols. Only in a few
### cases they are real non-transferrable locus idetnifers. We don't
### query those

rousettus_proteins <-  rousettus_proteins[!grepl("^LOC", rousettus_proteins$Locus), ]

allLocusGO <- getBM(mart=mart,
                    attributes = c("entrezgene_accession",
                                   "entrezgene_id",
                                   "entrezgene_description",
                                   "go_id"),
                    filters =  "entrezgene_accession",
                    values = rousettus_proteins$Locus) ###all values in "Locus"

saveRDS(allLocusGO, "intermediateData/GOtermAnnot.RDS")
