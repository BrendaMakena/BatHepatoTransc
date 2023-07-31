# plotting correlation plots for hepatocystis transcriptome intensities between spleen and liver

# loading the libraries
library(ggplot2)
library(GGally)
library(dplyr)
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
#metadata <- read.csv("inputdata/tagseqRNA_metadata2023.csv")

# loading from the metadata script
source("R/2_metadata.R")

# Check if the 'metadata' variable exists
if (exists("metadata")) {
  print(metadata)
} else {
  print("Metadata variable not found")
}



      #### ggplots for correlation between the transcripts #### 

# Scatter Plot of transcriptome infection estimates coloured by organ
pdf("plots/scatter_plot_of_transcriptome_infection_estimates_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
                  y = hepatocystis_transcriptome_parasitemia, 
                  color = Organ)) +
  geom_point() +
  scale_y_log10()+
  geom_smooth(method = "lm") +
  labs(title = "Scatter Plot of transcriptome infection estimates coloured by organ", 
       x = "blood parasitemia (% infected erythrocytes)", 
       y = "transcriptome infection estimate (#reads)") +
  theme_bw()

dev.off()

# Boxplot of transcriptome infection estimates coloured by organ
pdf("plots/Boxplot_of_transcriptome_infection_estimates_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Parasitemia_in_percent, 
           y = hepatocystis_transcriptome_parasitemia, 
           fill = Organ)) +
  geom_boxplot() +
  labs(title = "Boxplot of transcriptome infection estimates coloured by organ", 
       x = "blood parasitemia (% infected erythrocytes)", 
       y = "transcriptome infection estimate (#reads)") +
  theme_bw()

dev.off()

#parasitemia model
model_parasitemia <- glm.nb(hepatocystis_transcriptome_parasitemia ~ Parasitemia_in_percent * Organ,
       metadata)


summary(model_parasitemia)


# boxplot of transcripts correlations by infection status coloured by organ
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


# scatter plot of transcripts correlation by organ coloured by infectious status
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
       title = "Transcripts correlations by organ coloured by infection status")+
  theme_bw()

dev.off()

# scatter plot of blood parasitemia count coloured by organ
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
       title = "Blood parasitemia count coloured by organ") +
  theme_bw()

dev.off()


# box plot of infected vs uninfected coloured by organ
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
       title = "Transcriptome parasitemia in liver vs spleen samples based on infection status") +
  theme_bw()

dev.off()


# scatter plot of liver vs spleen parasitemia coloured by infected and uninfected
pdf("plots/scatter_plot_for_spleen_and_liver_infectious_status.pdf")
as_tibble(metadata) %>% 
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







