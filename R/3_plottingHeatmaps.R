# plotting heatmap of transcripts in the 169 samples - using pheatmap
# annotating the rows and columns - which are spleen, liver, sex, age,
# infection +/-, parasitemia intensity

# loading the libraries
library(pheatmap)
library(gplots) # for the heatmap.2 function
library(ggplot2)
library(GGally)
library(dplyr)
library(RColorBrewer)
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

colnames(metadata)
metadata$SampleID

colnames(tagseqRNAfeatureCounts)

metadata$ID <- paste(metadata$SampleID, 
                     substring(metadata$Organ, 1,1),
                     sep = ".")


# all(metadata$ID %in% colnames(tagseqRNAfeatureCounts))
# all(colnames(tagseqRNAfeatureCounts) %in% metadata$ID)

cbind(colnames(tagseqRNAfeatureCounts) == metadata$ID,
      colnames(tagseqRNAfeatureCounts), metadata$ID)

all(colnames(tagseqRNAfeatureCounts) == metadata$ID)

# we have technical replicates for some samples
# colnames in the count data have been made unique by R 
# by appending ".1" to them

metadata$ID <- make.unique(metadata$ID)

metadata <- metadata[order(metadata$ID),]

tagseqRNAfeatureCounts <- tagseqRNAfeatureCounts[,order(colnames(tagseqRNAfeatureCounts))]

table(grepl("\\.1",metadata$ID))

cbind(colnames(tagseqRNAfeatureCounts) == metadata$ID,
      colnames(tagseqRNAfeatureCounts), metadata$ID)

all(colnames(tagseqRNAfeatureCounts) == metadata$ID)


# getting the hepatocystis transcripts
# tagseqRNAfeatureCounts %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))

# getting the Rousettus transcripts
# tagseqRNAfeatureCounts %>% filter(!grepl("HEP_",GeneID,ignore.case = TRUE))


# getting the rowsums for the genes counts table
table(rowSums(tagseqRNAfeatureCounts))

# transcripts with counts >5 for both Hepatocystis and Rousettus
table(rowSums(tagseqRNAfeatureCounts)>5)

# transcripts with counts >5 for Hepatocystis
table(rowSums(tagseqRNAfeatureCounts)>5,
      grepl("HEP_",rownames(tagseqRNAfeatureCounts)))

# transcripts with counts >5 for Rousettus
table(rowSums(tagseqRNAfeatureCounts)>5,
      !grepl("HEP_",rownames(tagseqRNAfeatureCounts)))


# heatmap of 29 Hepatocystis transcripts >5
pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts)>5 &
                                        grepl("HEP_",rownames(tagseqRNAfeatureCounts)),]+1),
         show_colnames = T, show_rownames = T)


# getting the column sums for the hepatocystis transcripts
colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),])

# adding a new column for the hepatocystis transcriptome parasitemia to the metadata file
metadata$hepatocystis_transcriptome_parasitemia <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),])


# viewing column names for the feature counts table and metadata file
colnames(tagseqRNAfeatureCounts) 
colnames(metadata)


# heatmap of all transcripts
pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts), ]+1),
         show_colnames = T, show_rownames = T)

pheatmap(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 5, ],
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = brewer.pal(10, "RdYlBu"))


# Convert list to matrix
tagseqRNAfeatureCounts_matrix <- as.matrix(tagseqRNAfeatureCounts)

# Replace NA, Inf, and NaN values with a constant
processed_data <- replace(tagseqRNAfeatureCounts_matrix, is.na(tagseqRNAfeatureCounts_matrix) | is.infinite(tagseqRNAfeatureCounts_matrix) | is.nan(tagseqRNAfeatureCounts_matrix), 0.1)

# Perform log transformation on the processed data and generate the heatmap
pheatmap(log10(processed_data[rowSums(processed_data) > 5, ] + 1),
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "heatmap of all transcripts",
         fontsize = 10,
         fontsize_row = 8,
         fontsize_col = 8,
         width = 20,  # Adjusts the width of the plot
         height = 10) # Adjusts the height of the plot
         
pheatmap(log10(processed_data[rowSums(processed_data) > 10, ] + 1),
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "heatmap of all genes",
         fontsize = 10,
         fontsize_row = 4,
         fontsize_col = 6,
         width = 40,  # Adjusts the width of the plot
         height = 20) # Adjusts the height of the plot


### checking if the technical replicates have similar gene counts
    
    # 1st identifying the technical replicate columns and their corresponding copies

replicate_cols <- grep("\\.1$", colnames(tagseqRNAfeatureCounts), 
                       value = TRUE)

copy_cols <- sub("\\.1$", "", replicate_cols)

      # 2nd subseting the dataframe to include only the technical replicate columns and their corresponding copies

replicates_subset_df <- tagseqRNAfeatureCounts[, c(replicate_cols, copy_cols)]

    # Performing comparison and analysis on the replicates_subset_df

    # You can calculate summary statistics, perform t-tests, correlations, or any other relevant analysis to determine similarity

    # Example: Calculating mean values for each technical replicate and copy column

means <- apply(replicates_subset_df, 2, mean)

median <- apply(replicates_subset_df, 2, median)
    
# Printing the replicates_subset_df, means and medians for comparison

print(replicates_subset_df)

print(means)

print(median)




