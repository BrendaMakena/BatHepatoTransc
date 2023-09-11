#counting reads using features count after mapping with STAR

#https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html 

#first install and load the package Rsubread which has featuresCount as a built in package
#BiocManager::install("Rsubread")
library("Rsubread")
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)

#Setting the path to the directory containing the BAM files from STAR mapping
bamdirectory <- "/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/BAMfiles"
dir()

#Creating a list of BAM file names in the directory
bamFiles <- list.files(bamdirectory, pattern = ".bam$", full.names = TRUE)


#Perform feature counts to generate feature counts for each BAM file, specifying merged GTF

# Creating an empty data frame to store the feature counts
tagseqRNAfeatureCounts <- data.frame()

# List of files
bamFiles

# Loop over the files
for (file in bamFiles) {
  # Extract file name without path
  file_name <- tools::file_path_sans_ext(basename(file))
        
  #the basename function is used to extract the file name without the entire path. 
  #This file name is then used as the column name when adding the counts to the feature_counts data frame.
  
  # Running featureCounts for each file
 tagseqRNAcounts <- featureCounts(files = file, nthreads = 8, 
                               annot.ext = "/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/RousettusHepatocystis_merged.gtf", 
                isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id")

# Extracting the gene IDs and counts from the result
tagseqRNAgene_ids <- tagseqRNAcounts$annotation$GeneID
tagseqRNAfile_counts <- tagseqRNAcounts$counts

#Creating a new data frame with gene IDs and counts
tagseqRNAfile_data <- data.frame(GeneID = tagseqRNAgene_ids, tagseqRNAfile_counts)

#Merging the new data frame with the feature_counts data frame
if (nrow(tagseqRNAfeatureCounts) == 0) {
  tagseqRNAfeatureCounts <- tagseqRNAfile_data
} else {
  tagseqRNAfeatureCounts <- merge(tagseqRNAfeatureCounts, tagseqRNAfile_data, by = "GeneID", all = TRUE)
}
}


#Setting the gene IDs column as the row names of the feature counts table
rownames(tagseqRNAfeatureCounts) <- tagseqRNAfeatureCounts[,"GeneID"]

#Removing the gene ID column from the feature counts table
tagseqRNAfeatureCounts$GeneID <- NULL


#renaming the column names (sample IDs) to shorten them

#removing all unwanted characters from the end
names(tagseqRNAfeatureCounts) = gsub(pattern = "_S.*", 
                                     replacement = "", 
                                     x = names(tagseqRNAfeatureCounts))

#removing all unwanted characters from the start
names(tagseqRNAfeatureCounts) = gsub(pattern = "^.*D", 
                                     replacement = "D", 
                                     x = names(tagseqRNAfeatureCounts))

#removing transcripts with 0 counts
tagseqRNAfeatureCounts <- tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts)>0,]

# removing replicates by summing them up

# converting the matrix to a data frame and transposing the data
tagseqRNAfeatureCounts <- as.data.frame(t(tagseqRNAfeatureCounts))

# making row names unique by adding a suffix to duplicate sample IDs
rownames(tagseqRNAfeatureCounts) <- make.unique(rownames(tagseqRNAfeatureCounts))

# extracting common sample IDs to merge technical replicates
merged_counts_df <- tagseqRNAfeatureCounts %>%
  rownames_to_column(var = "gene_id") %>%
  mutate(sample_id = sub("\\.\\d+$", "", gene_id)) %>%
  dplyr::select(-gene_id) %>%
  group_by(sample_id) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

# transposing the data frame to get back the original format
tagseqRNAfeatureCounts <- as.data.frame(t(merged_counts_df[,-1]))
colnames(tagseqRNAfeatureCounts) <- merged_counts_df$sample_id


# Saving the feature counts to a file 
saveRDS(tagseqRNAfeatureCounts, "intermediateData/countTable.RDS")


