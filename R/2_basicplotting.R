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


#getting the hepatocystis transcripts
#tagseqRNAfeatureCounts %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))

#getting the Rousettus transcripts
#tagseqRNAfeatureCounts %>% filter(!grepl("HEP_",GeneID,ignore.case = TRUE))


#getting the rowsums for the genes counts table
table(rowSums(tagseqRNAfeatureCounts))

#transcripts with counts >5 for both Hepatocystis and Rousettus
table(rowSums(tagseqRNAfeatureCounts)>5)

#transcripts with counts >5 for Hepatocystis
table(rowSums(tagseqRNAfeatureCounts)>5,
      grepl("HEP_",rownames(tagseqRNAfeatureCounts)))

#transcripts with counts >5 for Rousettus
table(rowSums(tagseqRNAfeatureCounts)>5,
      !grepl("HEP_",rownames(tagseqRNAfeatureCounts)))


#heatmap of 29 Hepatocystis transcripts >5
pheatmap(log10(tagseqRNAfeatureCounts[rowSums(tagseqRNAfeatureCounts)>5 &
         grepl("HEP_",rownames(tagseqRNAfeatureCounts)),]+1),
         show_colnames = T, show_rownames = T)

#getting the column sums for the hepatocystis transcripts
colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),])

#adding a new column for the hepatocystis transcriptome parasitemia to the metadata file
metadata$hepatocystis_transcriptome_parasitemia <- colSums(tagseqRNAfeatureCounts[
  grepl("HEP_",rownames(tagseqRNAfeatureCounts)),])


#viewing column names for the feature counts table and metadata file
colnames(tagseqRNAfeatureCounts) 
colnames(metadata)

#viewing the sample IDs in the metadata file
(metadata$ID)

#xy plot
ggplot(metadata,aes(x = hepatocystis_transcriptome_parasitemia,y = Parasitemia_in_percent))


#getting the median 
tapply(metadata$hepatocystis_transcriptome_parasitemia,
       metadata$Parasitemia_in_percent,median)

#getting the mean
tapply(metadata$hepatocystis_transcriptome_parasitemia,
       metadata$Parasitemia_in_percent,mean)

#plotting box plot
ggplot(metadata, aes(y = hepatocystis_transcriptome_parasitemia, 
                     x = factor(Parasitemia_in_percent))) + 
  geom_boxplot() +
  scale_y_log10()


#making a ggpairs plot
ggpairs(metadata) #gives error because the 'Sample_ID_sequencing' column, 
              #has more levels (unique values) than the default cardinality threshold allows

#removing the column
#ggpairs(metadata[-which(names(metadata) == 'Sample_ID_sequencing')])

#alternatively increasing the cardinality threshold
#ggpairs(metadata, cardinality_threshold = 169) #this takes a lot of time plus outputs all columns even the unwanted

#to make ggpairs plot with only desired columns
selected_columns <- c("SampleID","Organ","Infectious_status_blood",
                      "Parasitemia_in_percent","Age_2category","Reproductive_status",
                      "Sex","Season","hepatocystis_transcriptome_parasitemia")

ggpairs(metadata[selected_columns],cardinality_threshold = 114)


