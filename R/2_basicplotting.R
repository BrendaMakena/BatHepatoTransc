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

#loading the dataframe (file containing the read counts) from the features count step
url <- "https://github.com/BrendaMakena/BatHepatoTransc/raw/main/intermediateData/countTable.RDS"
download.file(url, destfile = "countTable.RDS")
tagseqRNAfeatureCounts <- readRDS("countTable.RDS")


#Printing the loaded table
print(tagseqRNAfeatureCounts)

#renaming the column names (sample IDs) to shorten them
colnames(tagseqRNAfeatureCounts)
   
       #removing all unwanted characters from the end
names(tagseqRNAfeatureCounts) = gsub(pattern = "_S.*", 
                 replacement = "", 
                 x = names(tagseqRNAfeatureCounts))
    
       #removing all unwanted characters from the start
names(tagseqRNAfeatureCounts) = gsub(pattern = "^.*D", 
                                              replacement = "D", 
                                              x = names(tagseqRNAfeatureCounts))
#renaming column 1 for gene ID
colnames(tagseqRNAfeatureCounts)[1] <- "GeneID"

#viewing the table column names
colnames(tagseqRNAfeatureCounts)


#Setting the gene IDs column as the row names of the feature counts table
rownames(tagseqRNAfeatureCounts) <- tagseqRNAfeatureCounts[,1]

#Removing the gene ID column from the feature counts table
#tagseqRNAfeatureCounts <- tagseqRNAfeatureCounts[,-1]


#getting the hepatocystis transcripts
#tagseqRNAfeatureCounts %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))

#getting the Rousettus transcripts
#tagseqRNAfeatureCounts %>% filter(!grepl("HEP_",GeneID,ignore.case = TRUE))


#loading metadata file
metadata <- read.csv("https://raw.githubusercontent.com/BrendaMakena/BatHepatoTransc/main/inputdata/tagseqRNA_metadata2023.csv")
#alternatively importing from the server directory
#metadata <- read.csv("inputdata/tagseqRNA_metadata2023.csv")


#getting the rowsums for the genes counts table
table(rowSums(tagseqRNAfeatureCounts[-1]))

#transcripts with counts >5 for both Hepatocystis and Rousettus
table(rowSums(tagseqRNAfeatureCounts[-1])>5)

#transcripts with counts >5 for Hepatocystis
table(rowSums(tagseqRNAfeatureCounts[-1])>5,
      grepl("HEP_",tagseqRNAfeatureCounts$GeneID))

#transcripts with counts >5 for Rousettus
table(rowSums(tagseqRNAfeatureCounts[-1])>5,
      !grepl("HEP_",tagseqRNAfeatureCounts$GeneID))


#heatmap of 29 Hepatocystis transcripts >5
pheatmap(log10(tagseqRNAfeatureCounts[-1][rowSums(tagseqRNAfeatureCounts[-1])>5 &
         grepl("HEP_",tagseqRNAfeatureCounts$GeneID),]+1),
         show_colnames = F, show_rownames = T)

#getting the column sums for the hepatocystis transcripts
colSums(tagseqRNAfeatureCounts[-1][
  grepl("HEP_",tagseqRNAfeatureCounts$GeneID),])

#adding a new column for the hepatocystis transcriptome parasitemia to the metadata file
metadata$hepatocystis_transcriptome_parasitemia <- colSums(tagseqRNAfeatureCounts[-1][
  grepl("HEP_",tagseqRNAfeatureCounts$GeneID),])


#viewing column names for the feature counts table and metadata file
colnames(tagseqRNAfeatureCounts) 
colnames(metadata)

#viewing the sample IDs in the metadata file
(metadata$SampleID)

#xy plot
ggplot(metadata,aes(x = hepatocystis_transcriptome_parasitemia,y = Parasitemia_blood_percentage))


#getting the median 
tapply(metadata$hepatocystis_transcriptome_parasitemia,
       metadata$Parasitemia_blood_percentage,median)

#getting the mean
tapply(metadata$hepatocystis_transcriptome_parasitemia,
       metadata$Parasitemia_blood_percentage,mean)

ggplot(metadata,aes(y = hepatocystis_transcriptome_parasitemia,
                    x = Parasitemia_blood_percentage))+ 
  geom_boxplot() + 
  scale_y_log10()


ggpairs(metadata[-1]) #gives error


#Merging the hepatocystis transcripts feature counts dataframe with the metadata file using sample IDs
#merged_df <- merge(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)), 
#metadata, by.x = "SampleID", by.y = "SampleID")

#checking the column names for the files first
colnames(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)))
colnames(metadata)

#Transposing the transcript data to have samples as rows and transcripts as columns
#tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)), <- t(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)))

#Setting the gene IDs column as the row names of the hepatocystis transcript data
rownames(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))) <- tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))[,1]

#Removing the gene ID column from the hepatocystis transcript data
tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)) <- tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))[,-1]

#Extracting the predictor variables from the metadata that correspond to the columns of the hepatocystis transcript data
#first match column names between metadata and hepatocystis transcript data
#tagseq_matching_cols <- intersect(colnames(metadata), colnames(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))))

#Extracting the predictor variables from metadata
#tagseq_predictor_variables <- metadata[, matching_cols]

#extracting the predictor variables
sample_ids <- metadata$SampleID  # Extract the sample IDs from the "Sample ID" column
matching_cols <- sample_ids %in% colnames(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)))
predictor_variables <- tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))[, matching_cols]


# Convert the data frame to a matrix
mat <- as.matrix(tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE)))

# Extract the gene IDs from the row names of the matrix
gene_ids <- rownames(mat)

# Create a matrix for row annotations
row_annotations <- matrix(NA, nrow = nrow(mat), ncol = length(predictor_variables))
colnames(row_annotations) <- predictor_variables

# Set the gene IDs as row annotations
row_annotations[, "SampleID"] <- gene_ids

# Plot the heatmap with annotations
pheatmap(mat,
         annotation_row = row_annotations,
         annotation_col = predictor_variables)






