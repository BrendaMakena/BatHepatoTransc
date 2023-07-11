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
library(MASS)
library(tidyverse)

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

#creating a heatmap using the pheatmap function in R and annotating 
#it using the metadata file while coloring by the "Organ" column

# Extracting the relevant columns from the metadata file
metadata_subset <- metadata[, c("ID", "Organ")]

# Filtering the features in the feature counts file 
filtered_counts <- tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 500, ]

# Calculating the correlation matrix
correlation_matrix <- cor(log2(filtered_counts + 1))

# Defining the color palette for the "Organ" column
#color_palette <- c("spleen" = "darkgreen", "liver" = "darkblue")
#color_palette <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")


color_palette <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
pheatmap(correlation_matrix,
         annotation_col = metadata_subset,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         breaks = seq(min(correlation_matrix), max(correlation_matrix), length.out = length(color_palette) + 1),
         col = color_palette)


# Generating a smooth color palette
#smooth_color_palette <- colorRampPalette(color_palette)(100)


# Closing any existing plotting devices
dev.off()

# Reset the graphics system
graphics.off()

# Create a new plotting device
dev.new(width = 10, height = 10)

# Create a color palette for the organs
organ_palette <- c("Spleen" = "blue", "Liver" = "green")

# Generate the heatmap with color by organ
# Close any existing plotting devices
dev.off()

# Reset the graphics system
graphics.off()

# Create a new plotting device
dev.new(width = 10, height = 10)

# saving the plot to a file
pdf("plots/heatmap_with_metadata_annotations_coloured_by_organ.pdf")
png("plots/heatmap_with_metadata_annotations_coloured_by_organ.png")

# Generate the heatmap
heatmap.2(correlation_matrix,
          ColSideColors = color_palette[as.factor(metadata_subset$Organ)],
          main = "Transcripts Correlation",
          xlab = "Samples",
          ylab = "Transcripts",
          col = smooth_color_palette,
          key = TRUE,
          keysize = 1,
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cexRow = 0.8,
          cexCol = 0.8)

dev.off()


      #### ggplot for correlation between the transcripts #### 


# Selecting the columns for x-axis, y-axis, and color
x_column <- "Parasitemia_in_percent"
y_column <- "hepatocystis_transcriptome_parasitemia"
color_column <- "Organ" 

# Creating the ggplot as a scatter plot
pdf("plots/scatter plot of transcripts correlations coloured by organ.pdf")
#png("plots/scatter plot of transcripts correlations coloured by organ.png")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
                  y = hepatocystis_transcriptome_parasitemia, 
                  color = Organ)) +
  geom_point() +
  scale_y_log10()+
  geom_smooth(method = "lm") +
  labs(title = "Scatter Plot of transcripts correlations coloured by organ", 
       x = "blood parasitemia (% infected erythrocytes)", 
       y = "transcriptome infection estimate (#reads)")

dev.off()

# Creating the ggplot as a box plot
pdf("plots/Boxplot of transcripts correlations coloured by organ.pdf")
#png("plots/Boxplot of transcripts correlations coloured by organ.png")

ggplot(metadata, 
       aes_string(x = x_column, 
                  y = y_column, 
                  fill = color_column)) +
  geom_boxplot() +
  labs(title = "Boxplot of transcripts correlations coloured by organ", 
       x = x_column, 
       y = y_column)

dev.off()

model_parasitemia <- glm.nb(hepatocystis_transcriptome_parasitemia ~ Parasitemia_in_percent * Organ,
       metadata)


summary(model_parasitemia)


# Creating the ggplot as a scatter plot
pdf("plots/boxplot_of_parasitemia_infected_erythrocytes_coloured_by_organ.pdf")
#png("plots/scatter plot of transcripts correlations coloured by organ.png")

ggplot(metadata, 
       aes(x = Infectious_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "blood parasitemia detected", 
       y = "transcriptome infection estimate (#reads)")

dev.off()

# Creating the ggplot as a scatter plot
pdf("plots/scatter_plot_of_organ_paired_parasitemia.pdf")
#png("plots/scatter plot of transcripts correlations coloured by organ.png")

ggplot(metadata, 
       aes(x = Organ, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Infectious_status_blood)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "organ", 
       y = "transcriptome infection estimate (#reads)")

dev.off()

# Creating the ggplot as a scatter plot
pdf("plots/scatter_plot_of_blood_parasitemia_count.pdf")
#png("plots/scatter plot of transcripts correlations coloured by organ.png")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "blood parasitemia (% infected erythrocytes)", 
       y = "transcriptome infection estimate (#reads)")

dev.off()


# Creating the ggplot as a scatter plot
pdf("plots/boxplot_of_transcriptome_parasitemia_in_liver_samples.pdf")
#png("plots/scatter plot of transcripts correlations coloured by organ.png")

ggplot(metadata, 
       aes(x = Infectious_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "blood parasitemia detected", 
       y = "transcriptome infection estimate (#reads)")

dev.off()

pdf("plots/scatter_plot_for.pdf")

as_tibble(metadata) %>% 
  dplyr::select(hepatocystis_transcriptome_parasitemia,
                Organ, ID, Infectious_status_blood) %>%
  pivot_wider(values_from = c(hepatocystis_transcriptome_parasitemia),
              names_from = Organ,
              values_fill = NA) %>% 
  ggplot(aes(Liver, Spleen, color = Infectious_status_blood)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

dev.off()  


colnames(metadata)


