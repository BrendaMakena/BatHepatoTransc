# loading the metadata file and modifying it to remove replicates and 
# add new columns for hepatocystis transcriptome intensities

# loading the libraries
library(pheatmap)
library(gplots) # for the heatmap.2 function
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)
# rerun script 1 for generating counts table or read from intermediate data
readcount <- FALSE

# loading the dataframe (file containing the read counts) from the features count step
if(readcount){
  source("R/1_featurecounts.R")
}else{
  tagseqRNAfeatureCounts <- readRDS("intermediateData/countTable.RDS")
}

# loading metadata file
metadata <- read.csv("inputdata/tagseqRNA_metadata2023.csv")

# adding a new column for the sample IDs with their organ identifiers
metadata$ID <- paste(metadata$SampleID, 
                     substring(metadata$Organ, 1,1),
      sep = ".")

# we have technical replicates for some samples
# colnames in the feature counts data have been merged by summing
# to remove metadata replicate samples

# making row names unique by adding a suffix to the duplicate sample IDs
metadata$ID <- make.unique(metadata$ID)

# filtering out rows containing technical replicates
metadata <- metadata %>%
  filter(!grepl("\\.\\d+$", ID))

# ordering the metadata according to counts file
metadata <- metadata[order(metadata$ID),]

tagseqRNAfeatureCounts <- tagseqRNAfeatureCounts[,order(colnames(tagseqRNAfeatureCounts))]

if(!all(colnames(tagseqRNAfeatureCounts) == metadata$ID)){
  stop("metadata and count data are not alligned")
}

# adding a new column for the hepatocystis transcriptome parasitemia to the metadata file
metadata$hepatocystis_transcriptome_parasitemia <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),])

# adding columns for thresholds >2 for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_2 <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
  rowSums(tagseqRNAfeatureCounts) >2,])

# adding columns for thresholds >5 for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_5 <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
    rowSums(tagseqRNAfeatureCounts) >5,])

# getting the upper quartile threshold for hepatocystis transcriptome
uQ <- rowSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),]) %>%
  quantile(0.75)

# adding columns for thresholds below upper quartile for hepatocystis transcriptome
metadata$hepatocystis_transcriptome_parasitemia_uQ <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
    rowSums(tagseqRNAfeatureCounts)<uQ ,])

# adding column for sequencing depth for each sample
metadata$Sequencing_depth <- colSums(tagseqRNAfeatureCounts)

# metadata correction factor
metadata$mean_correction_factor <- mean(metadata$Sequencing_depth) / 
                                    metadata$Sequencing_depth 

# Renaming the column Infectious_status_blood
colnames(metadata)[colnames(metadata) == "Infectious_status_blood"] <- "Infection_status_blood"

# Renaming values in the "Infectious_status_blood" column
metadata <- metadata %>%
  mutate(
    Infection_status_blood = case_when(
      Infection_status_blood == 'uninfected' ~ 'undetected*',
      Infection_status_blood == 'infected' ~ 'infected*',
      TRUE ~ Infection_status_blood
    )
  )

# Calculating RPM values for each feature
rpm <- tagseqRNAfeatureCounts[!grepl("^HEP_", 
           rownames(tagseqRNAfeatureCounts)), ] %>%
  mutate(across(everything(),
                ~ . / colSums(tagseqRNAfeatureCounts[!grepl("^HEP_", 
         rownames(tagseqRNAfeatureCounts)),]) * 1e6))

# getting the rpmh values for the features 
rpmh <- tagseqRNAfeatureCounts[!grepl("^HEP_", 
        rownames(tagseqRNAfeatureCounts)), ] %>%
  mutate(across(everything(), #or mutate(across(starts_with("DMR"),
      ~ . / colSums(tagseqRNAfeatureCounts[!grepl("^HEP_", 
                rownames(tagseqRNAfeatureCounts)), ]) * 1e6))

# adding a new column for the rpmh to the metadata file
metadata$rpmh <- colSums(rpmh)



