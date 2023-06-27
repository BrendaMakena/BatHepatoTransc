#counting reads using features count after mapping with STAR

#https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/05_counting_reads.html 

#first install and load the package Rsubread which has featuresCount as a built in package
#BiocManager::install("Rsubread")
library("Rsubread")
featureCounts #to see the usage

#Setting the path to the directory containing the BAM files from STAR mapping
bamdirectory <- "/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/BAMfiles"
dir()

#Creating a list of BAM file names in the directory
bamFiles <- list.files(bamdirectory, pattern = ".bam$", full.names = TRUE)

#Setting the path to the merged GTF annotation files
#Hepatocystis gtf file download: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/459/845/GCA_902459845.2_HEP1/ 
#rouseatus gtf file download: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/176/215/GCF_014176215.1_mRouAeg1.p/GCF_014176215.1_mRouAeg1.p_genomic.gtf.gz 

#merged_gtf <- "/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/RousettusHepatocystis_merged.gtf"
#merged_gtf


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
 tagseqRNAcounts <- featureCounts(files = file, nthreads = 40, 
                               annot.ext = "/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/RousettusHepatocystis_merged.gtf", 
                isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id")

# Extracting the gene IDs and counts from the result
tagseqRNAgene_ids <- tagseqRNAcounts$annotation$GeneID
tagseqRNAfile_counts <- tagseqRNAcounts$counts

#Creating a new data frame with gene IDs and counts
tagseqRNAfile_data <- data.frame(GeneID = tagseqRNAgene_ids, tagseqRNAfile_counts)

#Merging the new data frame with the feature_counts data frame
if (nrow(tagseqRNAfeature_counts) == 0) {
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


# Saving the feature counts to a file 
saveRDS(tagseqRNAfeatureCounts, "intermediateData/countTable.RDS")


