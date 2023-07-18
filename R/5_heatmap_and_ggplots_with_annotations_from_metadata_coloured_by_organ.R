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



# heatmap of all transcripts
pdf("plots/heatmap_of_all_transcripts.pdf")
pheatmap(log10(tagseqRNAfeatureCounts +1),
         show_colnames = TRUE,
         show_rownames = TRUE)

dev.off()

# heatmap of all transcripts with greater than 5 counts
pdf("plots/heatmap_of_all_transcripts_with_greater_than_5_counts.pdf")
pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 5, ]+1),
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = brewer.pal(10, "RdYlBu"))

dev.off()
  

 
       #### correlation between transcripts ####

#creating a heatmap using the pheatmap function in R and annotating 
#it using the metadata file while coloring by the "Organ" column

pdf("plots/heatmap_of_all_transcripts_with_greater_than_500_counts_coloured_by_organ.pdf")

pheatmap(
  log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts) > 500, ] + 1),
  show_colnames = TRUE,
  show_rownames = TRUE,
  annotation_col = data.frame(Organ = metadata$Organ),
  annotation_colors = list(Organ = c("spleen" = "red", "liver" = "blue")),
  scale = "none",
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 8)

dev.off()



      #### ggplot for correlation between the transcripts #### 

# Creating the ggplot as a scatter plot
pdf("plots/scatter_plot_of_transcripts_correlations_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
                  y = hepatocystis_transcriptome_parasitemia, 
                  color = Organ)) +
  geom_point() +
  scale_y_log10()+
  geom_smooth(method = "lm") +
  labs(title = "Scatter Plot of transcripts correlations coloured by organ", 
       x = "blood parasitemia (% infected erythrocytes)", 
       y = "transcriptome infection estimate (#reads)") +
  theme_bw()

dev.off()

# Creating the ggplot as a box plot
pdf("plots/Boxplot_of_transcripts_correlations_coloured_by_organ.pdf")

ggplot(metadata, 
       aes_string(x = Parasitemia_in_percent, 
                  y = hepatocystis_transcriptome_parasitemia, 
                  fill = Organ)) +
  geom_boxplot() +
  labs(title = "Boxplot of transcripts correlations coloured by organ", 
       x = Parasitemia_in_percent, 
       y = hepatocystis_transcriptome_parasitemia) +
  theme_bw()

dev.off()


model_parasitemia <- glm.nb(hepatocystis_transcriptome_parasitemia ~ Parasitemia_in_percent * Organ,
       metadata)


summary(model_parasitemia)


# Creating the ggplot as a scatter plot
pdf("plots/boxplot_of_parasitemia_infected_erythrocytes_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Infectious_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "blood parasitemia detected", 
       y = "transcriptome infection estimate (#reads)",
       title = "Transcripts correlations coloured by organ") +
  theme_bw()

dev.off()


# Creating the ggplot as a scatter plot
pdf("plots/scatter_plot_of_organ_paired_parasitemia.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Infectious_status_blood)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "organ", 
       y = "transcriptome infection estimate (#reads)",
       title = "Transcripts correlations coloured by organ")+
  theme_bw()

dev.off()

# Creating the ggplot as a scatter plot
pdf("plots/scatter_plot_of_blood_parasitemia_count.pdf")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "blood parasitemia (% infected erythrocytes)", 
       y = "transcriptome infection estimate (#reads)",
       title = "Blood parasitemia count") +
  theme_bw()

dev.off()


# Creating the ggplot as a scatter plot
pdf("plots/boxplot_of_transcriptome_parasitemia_in_liver_vs_spleen_samples.pdf")

ggplot(metadata, 
       aes(x = Infectious_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "blood parasitemia detected", 
       y = "transcriptome infection estimate (#reads)",
       title = "Transcriptome parasitemia in liver vs spleen samples") +
  theme_bw()

dev.off()



pdf("plots/scatter_plot_for_spleen_and_liver_infectious_status.pdf")
as_tibble(metadata) %>% 
  filter(!SampleID %in% c("DMR970", "DMR992", "DMR993")) %>% 
  dplyr::select(hepatocystis_transcriptome_parasitemia,
                Organ, SampleID, Infectious_status_blood) %>%
  pivot_wider(values_from = c(hepatocystis_transcriptome_parasitemia),
              names_from = Organ,
              values_fill = NA) %>% 
  filter(!is.na(Liver) & !is.na(Spleen)) %>%
  ggplot(aes(Liver, Spleen, color = Infectious_status_blood)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10()+ 
  labs(title = "Liver vs Spleen Parasitemia") +
  theme_bw()


dev.off()







