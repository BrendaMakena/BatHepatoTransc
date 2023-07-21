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
source("R/2_metadata.R")


# heatmap of 29 Hepatocystis transcripts >5
pdf("plots/heatmap_of_29_Hepatocystis_transcripts_with_greater_than_5_counts.pdf")
pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts)>5 &
                                        grepl("HEP_",rownames(tagseqRNAfeatureCounts)),]+1),
         show_colnames = T, show_rownames = T)
dev.off()

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


pdf("plots/heatmap_all_greater_than_10_genes.pdf", 
    width = 40,
    height = 20)
         
pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 10, ] + 1),
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "heatmap of all genes",
         fontsize = 10,
         fontsize_row = 4,
         fontsize_col = 6)
         
dev.off()

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


Dist <- dist(log10(t(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 500, ])+1))

dim(Dist)
str(Dist)

pheatmap(Dist)

pheatmap(as.matrix(Dist)[technical_replicate_cols,])

rowSums(as.matrix(Dist)[technical_replicate_cols,])

as.matrix(Dist)[technical_replicate_cols,"DMR970.S"]

as.matrix(Dist)[,"DMR970.S"]

mean(as.matrix(Dist)[,"DMR970.S"])

min(as.matrix(Dist)[upper.tri(as.matrix(Dist))])

as.matrix(Dist)["DMR970.S.1","DMR970.S"]

as.matrix(Dist)["DMR992.L.1", "DMR992.L"]

as.matrix(Dist)["DMR993.L.1", "DMR993.L"]

fivenum(as.matrix(Dist)[upper.tri(as.matrix(Dist))])




technical_replicate_cols
plot(Dist)

rowSums(tagseqRNAfeatureCounts)
median(rowSums(tagseqRNAfeatureCounts))
mean(rowSums(tagseqRNAfeatureCounts))


ggplot(as.data.frame(rowSums(tagseqRNAfeatureCounts)+1), 
       aes(x = rowSums(tagseqRNAfeatureCounts)))+
  geom_histogram()+
  scale_y_log10()+
  scale_x_log10()


dim(tagseqRNAfeatureCounts)

max(rowSums(tagseqRNAfeatureCounts))
sum(tagseqRNAfeatureCounts)

tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts)==max(rowSums(tagseqRNAfeatureCounts)),]


# how deeply the samples are sequenced

ggplot(as.data.frame(colSums(tagseqRNAfeatureCounts)), 
       aes(x = colSums(tagseqRNAfeatureCounts)))+
  geom_histogram()


table(colSums(tagseqRNAfeatureCounts) <150000)

colnames(tagseqRNAfeatureCounts)[colSums(tagseqRNAfeatureCounts) <150000]
    #the two/3 samples to be excluded




