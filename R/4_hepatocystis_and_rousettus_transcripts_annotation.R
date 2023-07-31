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

# the merged gtf file
merged_gtf <- readGFF("/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/RousettusHepatocystis_merged.gtf")

# our hepatocystis transcripts with >5 counts
tagseqRNAfeatureCounts %>% 
  filter(!grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
         rowSums(tagseqRNAfeatureCounts)>5)

rowSums(tagseqRNAfeatureCounts %>% 
          filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
                 rowSums(tagseqRNAfeatureCounts)>5))


# are our 29 hepatocystis transcripts IDs present in the hepatocystis proteins table from NCBI
all(tagseqRNAfeatureCounts %>% 
  filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
         rowSums(tagseqRNAfeatureCounts)>5) == hepatocystis_proteins$Locus.tag)
            # they are not

# check if our feature counts gene IDs are in the merged gtf file
table(rownames(tagseqRNAfeatureCounts) %in% merged_gtf$gene_id)

# check the genes present in the merged gtf that our transcripts were mapped to
table(merged_gtf$gene_id %in% rownames(tagseqRNAfeatureCounts))

merged_gtf <- merged_gtf[merged_gtf$gene_id %in% rownames(tagseqRNAfeatureCounts),]
colnames(merged_gtf)

# exploring the data 

# number of transcripts according to type and biotype
table(merged_gtf$type)
table(merged_gtf$gene_biotype)

# adding species column to the gtf file
merged_gtf$species <- ifelse(grepl("HEP_", merged_gtf$gene_id),
                             "hepatocystis",
                             "rousettus")

# type, biotype and locus tag of transcripts according to species
table(merged_gtf$type, merged_gtf$species)
table(merged_gtf$gene_biotype, merged_gtf$species)
table(merged_gtf$locus_tag, merged_gtf$species)  #locus tag is only in hepatocystis transcripts


# mitochondrial genome not mapped to by our transcripts  
rowSums(tagseqRNAfeatureCounts %>% 
  filter(grepl("HEP_MIT",
               rownames(tagseqRNAfeatureCounts),
               ignore.case = TRUE)))

sum(grepl("HEP_MIT", rownames(tagseqRNAfeatureCounts), ignore.case = TRUE))

# motochondrial genome in the assembly (annotation file)
rowSums(merged_gtf %>% 
          filter(grepl("HEP_MIT",
                       merged_gtf$gene_id,
                       ignore.case = TRUE)))

sum(grepl("HEP_MIT", merged_gtf$gene_id, ignore.case = TRUE))


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


# selecting columns of interest in the merged gtf file
merged_gtf <- merged_gtf[c("gene_id")]
                             

# what is only expressed by our transcripts                             
# what proteins our hepatocystis transcripts translate to

rownames(tagseqRNAfeatureCounts %>% 
      filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
             rowSums(tagseqRNAfeatureCounts)>5)) %in% 
        merged_gtf$gene_id %in% 
  hepatocystis_proteins$Locus.tag
      # the transcripts ids are not present in the proteins table

# alternatively
hepatocystis_proteins %>%
filter(Locus.tag %in%
  merged_gtf$gene_id[merged_gtf$gene_id %in%
  (rownames(tagseqRNAfeatureCounts %>% 
           filter(grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
                  rowSums(tagseqRNAfeatureCounts)>0)))])
       
# annotations for rousettus transcripts       
rousettus_proteins %>%
  filter(Locus %in%
  merged_gtf$gene_id[merged_gtf$gene_id %in%
 (rownames(tagseqRNAfeatureCounts %>% 
 filter(!grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
 rowSums(tagseqRNAfeatureCounts)>0)))])

# selecting columns from the proteins annotations table
as.data.frame(rousettus_proteins %>%
  filter(Locus %in%
  merged_gtf$gene_id[merged_gtf$gene_id %in%
  (rownames(tagseqRNAfeatureCounts %>% 
  filter(!grepl("HEP_",rownames(tagseqRNAfeatureCounts),ignore.case = TRUE),
  rowSums(tagseqRNAfeatureCounts)>0)))]) %>%
  dplyr::select(Accession, Locus, Protein.Name))

                             

