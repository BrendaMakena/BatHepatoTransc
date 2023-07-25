# annotation of hepatocystis transcripts
library(rtracklayer)

# rerun script 1 for generating counts table or read from intermediate data
readcount <- FALSE

# loading the dataframe (file containing the read counts) from the features count step
if(readcount){
  source("R/1_featurecounts.R")
}else{
  tagseqRNAfeatureCounts <- readRDS("intermediateData/countTable.RDS")
}

# loading metadata file
#metadata <- read.csv("inputdata/tagseqRNA_metadata2023.csv")

# loading from the metadata script
source("R/2_metadata.R")

# Check if the 'metadata' variable exists
if (exists("metadata")) {
  print(metadata)
} else {
  print("Metadata variable not found")
}

# loading the proteins annotations table from ncbi
hepatocystis_proteins <- read.csv("inputdata/hepatocystis_proteins.csv")
rousettus_proteins <- read.csv("inputdata/rousettus_proteins.csv")

merged_gtf <- readGFF("/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/RousettusHepatocystis_merged.gtf")

# our hepatocystis transcripts with >5 counts
tagseqRNAfeatureCounts %>% 
  filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
         rowSums(tagseqRNAfeatureCounts)>5)

rowSums(tagseqRNAfeatureCounts %>% 
          filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
                 rowSums(tagseqRNAfeatureCounts)>5))



tagseqRNAfeatureCounts %>% 
  filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
         rowSums(tagseqRNAfeatureCounts)>5) == hepatocystis_proteins$Locus.tag


table(rownames(tagseqRNAfeatureCounts) %in% merged_gtf$gene_id)
table(merged_gtf$gene_id %in% rownames(tagseqRNAfeatureCounts))

merged_gtf <- merged_gtf[merged_gtf$gene_id %in% rownames(tagseqRNAfeatureCounts),]
colnames(merged_gtf)

table(merged_gtf$type)
table(merged_gtf$gene_biotype)

merged_gtf$species <- ifelse(grepl("HEP_", merged_gtf$gene_id),
                             "hepatocystis",
                             "rousettus")

table(merged_gtf$type, merged_gtf$species)
table(merged_gtf$gene_biotype, merged_gtf$species)
table(merged_gtf$locus_tag, merged_gtf$species)


# mitochondrial genome  


merged_gtf[merged_gtf$gene_id %in% "HEP_00048200",]

# only what in our dataset

                             
                             
                             

