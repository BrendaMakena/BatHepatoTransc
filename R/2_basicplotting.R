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
tagseq_RNA_feature_counts_table <- read.table("/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/tagseq_RNA_feature_counts.txt", 
                                              header = T, sep = "\t", row.names = 1)
#Printing the loaded table
print(tagseq_RNA_feature_counts_table)

#renaming the column names (sample IDs) to shorten them
colnames(tagseq_RNA_feature_counts_table)
   
       #removing all unwanted characters from the end
names(tagseq_RNA_feature_counts_table) = gsub(pattern = "_S.*", 
                 replacement = "", 
                 x = names(tagseq_RNA_feature_counts_table))
    
       #removing all unwanted characters from the start
names(tagseq_RNA_feature_counts_table) = gsub(pattern = "^.*D", 
                                              replacement = "D", 
                                              x = names(tagseq_RNA_feature_counts_table))
#renaming column 1 for gene ID
colnames(tagseq_RNA_feature_counts_table)[1] <- "GeneID"

#viewing the table
colnames(tagseq_RNA_feature_counts_table)


#Setting the gene IDs column as the row names of the feature counts table
rownames(tagseq_RNA_feature_counts_table) <- tagseq_RNA_feature_counts_table[,1]

#Removing the gene ID column from the feature counts table
#tagseq_RNA_feature_counts_table <- tagseq_RNA_feature_counts_table[,-1]


#getting the hepatocystis transcripts
tagseq_RNA_feature_counts_table %>% filter(grepl("HEP_",GeneID,ignore.case = TRUE))

#loading metadata file
metadata <- read.csv("/SAN/RNASeqHepatoCy/HostTranscriptome/host_merged_Raegyptiacus_and_HepatocystisAunin/STAR/featurescount/tagseqRNA_metadata2023.csv",
                     header = T,sep = ",")

#renaming sample id column
colnames(metadata)[1] = "SampleID"

#getting the rowsums for the genes counts table
table(rowSums(tagseq_RNA_feature_counts_table[-1]))

#transcripts with counts >5 for both Hepatocystis and Rousettus
table(rowSums(tagseq_RNA_feature_counts_table[-1])>5)

#transcripts with counts >5 for Hepatocystis
table(rowSums(tagseq_RNA_feature_counts_table[-1])>5,
      grepl("HEP_",tagseq_RNA_feature_counts_table$GeneID))

#transcripts with counts >5 for Rousettus
table(rowSums(tagseq_RNA_feature_counts_table[-1])>5,
      !grepl("HEP_",tagseq_RNA_feature_counts_table$GeneID))


#heatmap of 29 Hepatocystis transcripts >5
pheatmap(log10(tagseq_RNA_feature_counts_table[-1][rowSums(tagseq_RNA_feature_counts_table[-1])>5 &
         grepl("HEP_",tagseq_RNA_feature_counts_table$GeneID),]+1),
         show_colnames = F, show_rownames = T)

#getting the column sums for the hepatocystis transcripts
colSums(tagseq_RNA_feature_counts_table[-1][
  grepl("HEP_",tagseq_RNA_feature_counts_table$GeneID),])

#adding a new column for the hepatocystis transcriptome parasitemia to the metadata file
metadata$hepatocystis_transcriptome_parasitemia <- colSums(tagseq_RNA_feature_counts_table[-1][
  grepl("HEP_",tagseq_RNA_feature_counts_table$GeneID),])


#viewing column names for the feature counts table
colnames(tagseq_RNA_feature_counts_table) 

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






