#Plotting with predictor variables eg; age and sex
#pairwise correlation between the two variables using eg ggpairs or pairs
#to determine best mapping of Hepatocystis transcripts:

#1st plotting heatmap of Hepatocystis transcripts in the 169 samples
#using pheatmap
#annotating the rows and columns - which are spleen, liver, sex, age,
#infection +/-, parasitemia intensity

#loading the libraries
library(pheatmap)
library(gplots) # for the heatmap.2 function
library(ggplot2)
library(GGally)
library(dplyr)

#rerun script 1 for generating counts table or read from intermediate data
readcount <- FALSE

#loading the dataframe (file containing the read counts) from the features count step
if(readcount){
  source("R/1_featurecounts.R")
}else{
  tagseqRNAfeatureCounts <- readRDS("intermediateData/countTable.RDS")
}

#loading metadata file
metadata <- read.csv("inputdata/tagseqRNA_metadata2023.csv")

metadata$ID <- paste(metadata$SampleID, 
                     substring(metadata$Organ, 1,1),
      sep = ".")

# we have technical replicates for some samples
# colnames in the count data have been made unique by R 
# by appending ".1" to them

metadata$ID <- make.unique(metadata$ID)

metadata <- metadata[order(metadata$ID),]

tagseqRNAfeatureCounts <- tagseqRNAfeatureCounts[,order(colnames(tagseqRNAfeatureCounts))]

if(!all(colnames(tagseqRNAfeatureCounts) == metadata$ID)){
  stop("metadata and count data are not alligned")
}


#adding a new column for the hepatocystis transcriptome parasitemia to the metadata file
metadata$hepatocystis_transcriptome_parasitemia <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),])

#adding columns for thresholds >2 for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_2 <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
  rowSums(tagseqRNAfeatureCounts) >2,])

#adding columns for thresholds >5 for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_5 <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
    rowSums(tagseqRNAfeatureCounts) >5,])


#adding columns for thresholds below upper quartile for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_5 <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
    rowSums(tagseqRNAfeatureCounts) ,])

uQ <- rowSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),]) %>%
  quantile(0.75)


#adding columns for thresholds below upper quartile for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_uQ <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
    rowSums(tagseqRNAfeatureCounts)<uQ ,])


