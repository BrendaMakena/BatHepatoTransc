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
library(stats)

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

# viewing the column names
colnames(metadata)
metadata$SampleID

colnames(tagseqRNAfeatureCounts)

# adding ID column in metadata file 
metadata$ID <- paste(metadata$SampleID, 
                     substring(metadata$Organ, 1,1),
                     sep = ".")

# confirming ID column in metadata matches column names in feature counts file
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

pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 5, ]+1),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = brewer.pal(10, "RdYlBu"))

#correlation between transcripts

#calculating the Euclidean distance between transcripts
euclidean_distance <- dist(t(tagseqRNAfeatureCounts), method = "euclidean")
euclidean_distance

#for all transcripts - log transformed
all_transcripts_correlation<- dist(log10(t(tagseqRNAfeatureCounts)+1))
all_transcripts_correlation
dim(all_transcripts_correlation)
str(all_transcripts_correlation)

pheatmap(all_transcripts_correlation)

#for transcripts with rowsums >500 log transformed
transcripts_correlation<- dist(log10(t(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 500, ])+1))
transcripts_correlation

dim(transcripts_correlation)
str(transcripts_correlation)

pheatmap(transcripts_correlation)

transcripts_correlation_matrix <- as.matrix(1 - as.matrix(transcripts_correlation))
transcripts_correlation_matrix
dim(transcripts_correlation_matrix)  #gives the number of rows and columns in the matrix, representing the pairwise distances between transcripts
str(transcripts_correlation_matrix)

pheatmap(transcripts_correlation_matrix)


#making the heatmap with metadata added for annotation and clustered by organ

# Specifying the column order for clustering based on the 'Organ' column
metadata_ordered <- metadata[order(metadata$Organ), , drop = FALSE]

# Removing rows with NA values
metadata_ordered <- metadata_ordered[complete.cases(metadata_ordered), ]
#str(metadata_ordered)

# Extracting unique organ names from metadata file 
unique_organs <- unique(metadata_ordered$Organ)
#print(unique_organs)

# Generating a sample color palette
#n_colors <- 5
#color_palette <- viridis(n_colors)
# Printing the color palette
#print(color_palette)

# Defining the number of unique organs
n_colors <- length(unique_organs)

# Defining the color palette
color_palette <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
#color_palette <- c("darkgreen", "darkblue")  # Custom colors for "Liver" and "Spleen"
#color_palette <- c("spleen" = "darkgreen", "liver" = "darkblue")

#print(color_palette)

# Creating the heatmap
pheatmap(transcripts_correlation_matrix, 
         annotation_col = metadata_ordered,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         color = list(heatmap = color_palette))

  
  # Creating the heatmap using the heatmap.2 function

# Clustering the columns based on the "Organ" column in metadata_ordered
organ_clusters <- hclust(dist(transcripts_correlation_matrix))
column_order <- order.dendrogram(as.dendrogram(organ_clusters))

# Create the heatmap using the reordered columns
heatmap.2(as.matrix(transcripts_correlation_matrix),
          Rowv = TRUE,
          Colv = column_order,
          col = color_palette,
          main = "Transcripts Correlation",
          xlab = "Samples",
          ylab = "Transcripts",
          margins = c(7, 10),
          key = TRUE,
          keysize = 1,
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          labRow = row_names,
          labCol = col_names,
          srtCol = 45)

dev.off()


# calculating correlation coefficient between transcripts using 
# Pearson's correlation coefficient
pearsons_correlation_matrix <- cor(tagseqRNAfeatureCounts, 
                                   method = "pearson" )

pearsons_correlation_matrix
pheatmap(pearsons_correlation_matrix)

# calculating correlation coefficient between transcripts using 
# spearman's correlation coefficient
spearmans_correlation_matrix <- cor(tagseqRNAfeatureCounts, 
                                    method = "spearman")

spearmans_correlation_matrix
pheatmap(spearmans_correlation_matrix)


