### plotting correlation plots for hepatocystis transcriptome
### intensities between spleen and liver

### loading the libraries
library(ggplot2)
library(GGally)
library(dplyr)
library(MASS)
library(tidyverse)
library(ggeffects)


# rerun script 1 for generating counts table or read from intermediate data
readcount <- FALSE

# loading the dataframe (file containing the read counts) from the features count step
if(readcount){
  source("scripts/3_featurecounts.R")
}else{
  tagseqRNAfeatureCounts <- readRDS("intermediateData/countTable.RDS")
}

# loading metadata file
metadata <- read.csv("intermediateData/metadata_expanded.csv")


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
       x = "blood parasitemia in % 
       (percent of infected erythrocytes)", 
       y = "# reads of Hepatocystis transcripts)") +
  theme_bw()

dev.off()

# Boxplot of transcriptome infection estimates coloured by organ
pdf("plots/Boxplot_of_transcriptome_infection_estimates_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = hepatocystis_transcriptome_parasitemia, 
           fill = Organ)) +
  geom_boxplot() +
  labs(title = "Boxplot of transcriptome infection estimates coloured by organ", 
       y = "transcriptome infection estimate (#reads)") +
  theme_bw()

dev.off()

#parasitemia model
model_parasitemia <- glm.nb(hepatocystis_transcriptome_parasitemia ~ 
                              Parasitemia_in_percent * Organ + 
                            offset(log(Sequencing_depth)),
       metadata)

summary(model_parasitemia)


pdf("plots/model_parasitemia.pdf")
ggpredict(model_parasitemia, 
              terms = c("Parasitemia_in_percent", "Organ"),
          condition = c(Sequencing_depth = mean(metadata$Sequencing_depth))) %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group),
              alpha = .1) +
  #ylim(c(0, 2e+07)) +
  geom_point(data = metadata, aes(x = Parasitemia_in_percent,
             y = hepatocystis_transcriptome_parasitemia*
             mean_correction_factor,
             color = Organ)) +
  scale_y_log10() +
  labs(x = "blood parasitemia in %
       (percentage of infected erythrocytes)",
       y = "# of reads of Hepatocystis transcripts") +
  theme_bw()

dev.off()
 

# boxplot of transcripts correlations by infection status coloured by organ
pdf("plots/boxplot_of_parasitemia_infected_erythrocytes_coloured_by_organ.pdf")

ggplot(metadata, 
       aes(x = Infection_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "infection status
       (*by microscopy of blood smears and PCR)", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcripts correlations coloured by organ") +
  theme_bw()

dev.off()

# scatter plot of transcripts correlation by organ coloured by infection status
pdf("plots/scatter_plot_of_organ_paired_parasitemia.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Infection_status_blood)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "organ", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcripts correlations by organ coloured by infection status") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold"))
  
dev.off()

# scatter plot of transcripts correlation by organ coloured by infection status and using the rpmh values
pdf("plots/scatter_plot_of_organ_paired_parasitemia_using_rpmh_values.pdf")

ggplot(metadata, 
       aes(x = Organ, 
           y = rpmh, 
           color = Infection_status_blood)) +
  geom_point()+
  geom_line(aes(group = SampleID))+
  scale_y_log10()+
  labs(x = "organ", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcripts correlations by organ coloured 
       by infection status (using rpmh values)") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold"))

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
  labs(x = "blood parasitemia in % 
       (percent of infected erythrocytes)", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Blood parasitemia count coloured by organ") +
  theme_bw()

dev.off()


# box plot of infected vs undetected coloured by organ
pdf("plots/boxplot_of_transcriptome_intensity_in_liver_vs_spleen_samples.pdf")

ggplot(metadata, 
       aes(x = Infection_status_blood, 
           y = hepatocystis_transcriptome_parasitemia, 
           color = Organ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25))+
  scale_y_log10()+
  labs(x = "infection status
       (* by microscopy of blood smears and PCR)", 
       y = "# of reads of Hepatocystis transcripts",
       title = "Transcriptome intensity in liver vs spleen samples based on infection status") +
  theme_bw()

dev.off()


# scatter plot of liver vs spleen parasitemia coloured by infected and uninfected 
#for samples with both spleen and liver tissues
pdf("plots/scatter_plot_for_spleen_and_liver_infection_status_in_samples_with_both_tissues.pdf")
as_tibble(metadata) %>% 
  dplyr::select(hepatocystis_transcriptome_parasitemia,
                Organ, SampleID, Infection_status_blood) %>%
  pivot_wider(values_from = c(hepatocystis_transcriptome_parasitemia),
              names_from = Organ,
              values_fill = NA) %>% 
  filter(!is.na(Liver) & !is.na(Spleen)) %>%
  ggplot(aes(Liver, Spleen, color = Infection_status_blood)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "gray") +  # Reference line
  scale_y_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  scale_x_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  labs(title = "Liver vs spleen parasitemia in samples with both tissues") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold")) #+
# coord_fixed(ratio = 1) +
# xlim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE))) +
# ylim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE)))

dev.off()


# scatter plot of liver vs spleen parasitemia coloured by infected and uninfected using rpmh values 
#for samples with both spleen and liver tissues
pdf("plots/scatter_plot_for_spleen_and_liver_infection_status_in_samples_with_both_tissues_using_rpmh_values.pdf")
as_tibble(metadata) %>% 
  dplyr::select(rpmh,
                Organ, SampleID, Infection_status_blood) %>%
  pivot_wider(values_from = c(rpmh),
              names_from = Organ,
              values_fill = NA) %>% 
  filter(!is.na(Liver) & !is.na(Spleen)) %>%
  ggplot(aes(Liver, Spleen, color = Infection_status_blood)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "gray") +  # Reference line
  scale_y_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  scale_x_log10(expand = c(0, 0)) +  # Set expand argument to avoid extra space 
  labs(title = "Liver vs spleen parasitemia in samples with both tissues using rpmh values") +
  theme_bw() +
  scale_color_manual(values = c("infected*" = "green", "undetected*" = "gold")) #+
# coord_fixed(ratio = 1) +
# xlim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE))) +
# ylim(c(1, max(metadata$hepatocystis_transcriptome_parasitemia, na.rm = TRUE)))

dev.off()



