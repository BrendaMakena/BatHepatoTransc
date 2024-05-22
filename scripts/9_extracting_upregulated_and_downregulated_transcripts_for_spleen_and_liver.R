# Script for extracting significant DETs, upregulated and downregulated DETs from the spleen and liver
# 3'Tagseqs mapped to de novo assembled transcriptome of E. labiatus infected with Hepatocystis


#### Extracting up and down regulated transcripts in each category in the liver and spleen samples ####

# getting the results table for liver DETs
list_of_results_liver  <- lapply(resultsNames(dds_liver), function(n){
  results(dds_liver, name = n)
})
names(list_of_results_liver) <- paste0("liver:", resultsNames(dds_liver))

# getting the results table for spleen DETs
list_of_results_spleen  <- lapply(resultsNames(dds_spleen), function(n){
  results(dds_spleen, name = n)
})
names(list_of_results_spleen) <- paste0("spleen:", resultsNames(dds_spleen))

## combined liver and spleen list of DETs results
list_of_results <- c(list_of_results_liver, list_of_results_spleen)




#### Creating tables showing numbers of upregulated and downregulated transcripts per category #### 

# For condition A vs condition B

### For liver DETs ### 

# Initializing variables to store counts
upregulated_liver_counts <- numeric(length(list_of_results_liver) - 1)  # excluding the first result which is the baseMean
downregulated_liver_counts <- numeric(length(list_of_results_liver) - 1)

# Iterating through each category except the first one
for (i in 2:length(list_of_results_liver)) {
  result_liver <- list_of_results_liver[[i]]  # Gets results for the current category
  
  # Extracting significant rows based on adjusted p-value threshold
  significant_liver_rows <- result_liver[!is.na(result_liver$padj) & result_liver$padj < 0.01, ]
  
  # Counting upregulated and downregulated transcripts
  upregulated_liver_counts[i - 1] <- sum(significant_liver_rows$log2FoldChange > 0)
  downregulated_liver_counts[i - 1] <- sum(significant_liver_rows$log2FoldChange < 0)
}

# Printing the counts
print(upregulated_liver_counts)
print(downregulated_liver_counts)

# Combining the counts vectors into a data frame
liver_counts_df <- data.frame(Upregulated = upregulated_liver_counts, Downregulated = downregulated_liver_counts)

# Setting the row names for the counts data frame
rownames(liver_counts_df) <- resultsNames(dds_liver)[-1]  # Exclude the first category

# Transposing the data frame to have categories as headers
counts_table_liver <- t(liver_counts_df)

# Printing the table
print(counts_table_liver)

# Saving the table to a CSV file
write.table(counts_table_liver, file = "/home/brenda/BatHepatoTransc/intermediateData/LiverDETS_table_by_category.csv", 
            sep = ",", col.names = NA)


### For spleen DETs ### 

# Initializing variables to store counts
upregulated_spleen_counts <- numeric(length(list_of_results_spleen) - 1)  # excluding the first result which is the baseMean
downregulated_spleen_counts <- numeric(length(list_of_results_spleen) - 1)

# Iterating through each category except the first one
for (i in 2:length(list_of_results_spleen)) {
  result_spleen <- list_of_results_spleen[[i]]  # Gets results for the current category
  
  # Extracting significant rows based on adjusted p-value threshold
  significant_spleen_rows <- result_spleen[!is.na(result_spleen$padj) & result_spleen$padj < 0.01, ]
  
  # Counting upregulated and downregulated transcripts
  upregulated_spleen_counts[i - 1] <- sum(significant_spleen_rows$log2FoldChange > 0)
  downregulated_spleen_counts[i - 1] <- sum(significant_spleen_rows$log2FoldChange < 0)
}

# Printing the counts
print(upregulated_spleen_counts)
print(downregulated_spleen_counts)

# Combining the counts vectors into a data frame
spleen_counts_df <- data.frame(Upregulated = upregulated_spleen_counts, Downregulated = downregulated_spleen_counts)

# Setting the row names for the counts data frame
rownames(spleen_counts_df) <- resultsNames(dds_spleen)[-1]  # Exclude the first category

# Transposing the data frame to have categories as headers
counts_table_spleen <- t(spleen_counts_df)

# Printing the table
print(counts_table_spleen)


# Saving the table to a CSV file
write.table(counts_table_spleen, file = "/home/brenda/BatHepatoTransc/intermediateData/SpleenDETS_table_by_category.csv", 
            sep = ",", col.names = NA)



#### To extract up and downregulated genes in each category in spleen ####

# For condition A vs condition B

# Initializing lists to store upregulated and downregulated genes for each category in spleen
upregulated_spleen_genes_list <- list()
downregulated_spleen_genes_list <- list()

# Iterating through each category
for (i in 2:length(list_of_results_spleen)) {
  result_spleen <- list_of_results_spleen[[i]]  # Gets results for the current category
  
  # Extracting significant rows based on adjusted p-value threshold
  significant_spleen_rows <- result_spleen[!is.na(result_spleen$padj) & result_spleen$padj < 0.01, ]
  
  # Extracting upregulated and downregulated genes
  upregulated_spleen_genes <- significant_spleen_rows[significant_spleen_rows$log2FoldChange > 0, ]
  downregulated_spleen_genes <- significant_spleen_rows[significant_spleen_rows$log2FoldChange < 0, ]
  
  # Storing the genes in the respective lists
  upregulated_spleen_genes_list[[i - 1]] <- upregulated_spleen_genes
  downregulated_spleen_genes_list[[i - 1]] <- downregulated_spleen_genes
}

# Printing or further processing the lists of genes
print("Upregulated spleen Genes List:")
print(upregulated_spleen_genes_list)

print("Downregulated spleen Genes List:")
print(downregulated_spleen_genes_list)


# Writing the lists to a CSV file
lapply(upregulated_spleen_genes_list, function(x) write.table(x, file = "/home/brenda/BatHepatoTransc/intermediateData/upregulated_spleen_genes_list.csv", 
                                                              sep = ",", col.names = NA, append = TRUE))


lapply(downregulated_spleen_genes_list, function(x) write.table(x, file = "/home/brenda/BatHepatoTransc/intermediateData/downregulated_spleen_genes_list.csv", 
                                                                sep = ",", col.names = NA, append = TRUE))


#### To extract up and downregulated genes in each category in liver ####


# Initializing lists to store upregulated and downregulated genes for each category in liver
upregulated_liver_genes_list <- list()
downregulated_liver_genes_list <- list()

# Iterating through each category
for (i in 2:length(list_of_results_liver)) {
  result_liver <- list_of_results_liver[[i]]  # Gets results for the current category
  
  # Extracting significant rows based on adjusted p-value threshold
  significant_liver_rows <- result_liver[!is.na(result_liver$padj) & result_liver$padj < 0.01, ]
  
  # Extracting upregulated and downregulated genes
  upregulated_liver_genes <- significant_liver_rows[significant_liver_rows$log2FoldChange > 0, ]
  downregulated_liver_genes <- significant_liver_rows[significant_liver_rows$log2FoldChange < 0, ]
  
  # Storing the genes in the respective lists
  upregulated_liver_genes_list[[i - 1]] <- upregulated_liver_genes
  downregulated_liver_genes_list[[i - 1]] <- downregulated_liver_genes
}

# Printing or further processing the lists of genes
print("Upregulated liver Genes List:")
print(upregulated_liver_genes_list)

print("Downregulated liver Genes List:")
print(downregulated_liver_genes_list)


# Writing the lists to a CSV file
lapply(upregulated_liver_genes_list, function(x) write.table(x, file = "/home/brenda/BatHepatoTransc/intermediateData/upregulated_liver_genes_list.csv", 
                                                             sep = ",", col.names = NA, append = TRUE))


lapply(downregulated_liver_genes_list, function(x) write.table(x, file = "/home/brenda/BatHepatoTransc/intermediateData/downregulated_liver_genes_list.csv", 
                                                               sep = ",", col.names = NA, append = TRUE))


# To summarize the counts of differentially expressed transcripts (DETs) for each condition separately

# Category names
category_names <- c("Season_Rainy_vs_Dry", "Age_2category_Young_vs_Adult", 
                    "Sex_Male_vs_Female", "rpmh_scaled")


# DETs counts summary for liver per category


# Initializing a vector to store the counts for each condition
liver_DETs_summary_counts <- numeric(length(category_names))

# Iterating through each condition and count the number of DETs
for (i in seq_along(category_names)) {
  liver_DETs_summary_counts[i] <- nrow(significant_liver_rows_list[[i + 1]])
}

# Creating a summary data frame
liver_DETs_counts_summary_df <- data.frame(Category = category_names, Count = liver_DETs_summary_counts)

# Print the summary data frame
print(liver_DETs_counts_summary_df)


# Saving the table to a CSV file
write.table(liver_DETs_counts_summary_df, file = "/home/brenda/BatHepatoTransc/intermediateData/liver_DETs_counts_summary.csv", 
            sep = ",", col.names = NA)



# DETs counts summary for spleen per category

# Initializing a vector to store the counts for each condition
spleen_DETs_summary_counts <- numeric(length(category_names))

# Iterating through each condition and count the number of DETs
for (i in seq_along(category_names)) {
  spleen_DETs_summary_counts[i] <- nrow(significant_spleen_rows_list[[i + 1]])
}

# Creating a summary data frame
spleen_DETs_counts_summary_df <- data.frame(Category = category_names, Count = spleen_DETs_summary_counts)

# Print the summary data frame
print(spleen_DETs_counts_summary_df)

# Saving the table to a CSV file
write.table(spleen_DETs_counts_summary_df, file = "/home/brenda/BatHepatoTransc/intermediateData/spleen_DETs_counts_summary.csv", 
            sep = ",", col.names = NA)



#### To extract the significant results in the combined categories ####

# Extracting significant liver rows based on adjusted p-value threshold for combined categories
significant_liver_rows <- result_liver[!is.na(result_liver$padj) & result_liver$padj < 0.01, ]

# significant liver results
print(significant_liver_rows)

# Saving the significant transcripts to a CSV file
write.csv(significant_liver_rows, file = "/home/brenda/BatHepatoTransc/intermediateData/significant_liver_transcripts.csv",
          row.names = TRUE)


# Extracting significant spleen rows based on adjusted p-value threshold
significant_spleen_rows <- result_spleen[!is.na(result_spleen$padj) & result_spleen$padj < 0.01, ]

# significant spleen results
print(significant_spleen_rows)

# Saving the significant transcripts to a CSV file
write.csv(significant_spleen_rows, file = "/home/brenda/BatHepatoTransc/intermediateData/significant_spleen_transcripts.csv",
          row.names = TRUE)



#### To extract the significant results in each category ####

# Significant spleen transcripts per category

# Initializing lists to store significant rows for each category
significant_spleen_rows_list <- list()

# Iterating through each category
for (i in 1:length(list_of_results_spleen)) {
  result_spleen <- list_of_results_spleen[[i]]  # Gets results for the current category
  
  # Extracting significant DETs based on adjusted p-value and log2FoldChange thresholds
  significant_spleen_DETs <- result_spleen[!is.na(result_spleen$padj) & result_spleen$padj < 0.01, ] #& abs(result_spleen$log2FoldChange) > 1, ]
  
  # Storing significant rows in the list
  significant_spleen_rows_list[[i]] <- significant_spleen_DETs
}

#str(all_significant_spleen_rows)
str(significant_spleen_rows_list)


# Printing the length of each element in significant_spleen_rows_list
for (i in 1:length(significant_spleen_rows_list)) {
  cat(names(list_of_results_spleen)[i], "has", 
      nrow(significant_spleen_rows_list[[i]]), "significant rows.\n")
}

# Creating a list to store the categories
spleen_sig_DETs_categories_list <- list()

# Add each category's data frame to the list
for (i in 2:5) {
  spleen_sig_DETs_category_name <- names(list_of_results_spleen)[i]
  spleen_sig_DETs_categories_list[[spleen_sig_DETs_category_name]] <- as.data.frame(significant_spleen_rows_list[[i]])
}

# Printing the list
print(spleen_sig_DETs_categories_list)

# Saving to a csv file
lapply(spleen_sig_DETs_categories_list, 
       function(x) write.table(x, file = "/home/brenda/BatHepatoTransc/intermediateData/spleen_significant_DETs_per_category_list.csv", 
                               sep = ",", col.names = NA, append = TRUE))





# Significant liver transcripts per category

# Initializing lists to store significant rows for each category
significant_liver_rows_list <- list()

# Iterating through each category
for (i in 1:length(list_of_results_liver)) {
  result_liver <- list_of_results_liver[[i]]  # Gets results for the current category
  
  # Extracting significant DETs based on adjusted p-value and log2FoldChange thresholds
  significant_liver_DETs <- result_liver[!is.na(result_liver$padj) & result_liver$padj < 0.01, ] #& abs(result_liver$log2FoldChange) > 1, ]
  
  # Storing significant rows in the list
  significant_liver_rows_list[[i]] <- significant_liver_DETs
}

#str(all_significant_liver_rows)
str(significant_liver_rows_list)


# Printing the length of each element in significant_liver_rows_list
for (i in 1:length(significant_liver_rows_list)) {
  cat(names(list_of_results_liver)[i], "has", 
      nrow(significant_liver_rows_list[[i]]), "significant rows.\n")
}

# Creating a list to store the categories
liver_sig_DETs_categories_list <- list()

# Add each category's data frame to the list
for (i in 2:5) {
  liver_sig_DETs_category_name <- names(list_of_results_liver)[i]
  liver_sig_DETs_categories_list[[liver_sig_DETs_category_name]] <- as.data.frame(significant_liver_rows_list[[i]])
}

# Printing the list
print(liver_sig_DETs_categories_list)

# Saving to a csv file
lapply(liver_sig_DETs_categories_list, 
       function(x) write.table(x, file = "/home/brenda/BatHepatoTransc/intermediateData/liver_significant_DETs_per_category_list.csv", 
                               sep = ",", col.names = NA, append = TRUE))



#### For condition B vs condition A #### 


# Counting upregulated and downregulated genes separately for condition B (upregulated in condition B compared to condition A)

# For liver DETs in condition B vs A

# Initializing variables to store counts
upregulated_liver_counts_BvsA <- numeric(length(list_of_results_liver) - 1)  # excluding the first result which is the baseMean
downregulated_liver_counts_BvsA <- numeric(length(list_of_results_liver) - 1)

# Iterating through each category except the first one
for (i in 2:length(list_of_results_liver)) {
  result_liver <- list_of_results_liver[[i]]  # Gets results for the current category
  
  # Inverting the log2FoldChange to switch the comparison from A vs B to B vs A
  #result_liver$log2FoldChange <- -result_liver$log2FoldChange
  
  # Extracting significant rows based on adjusted p-value threshold
  significant_liver_rows <- result_liver[!is.na(result_liver$padj) & result_liver$padj < 0.01, ]
  
  # Counting upregulated and downregulated transcripts for condition B vs A
  upregulated_liver_counts_BvsA[i - 1] <- sum(significant_liver_rows$log2FoldChange < 0)  # Reverse the comparison
  downregulated_liver_counts_BvsA[i - 1] <- sum(significant_liver_rows$log2FoldChange > 0)  # Reverse the comparison
}

# Printing the counts
print(upregulated_liver_counts_BvsA)
print(downregulated_liver_counts_BvsA)

# Combining the counts vectors into a data frame
liver_counts_BvsA_df <- data.frame(Upregulated = upregulated_liver_counts_BvsA, Downregulated = downregulated_liver_counts_BvsA)

# Setting the row names for the counts data frame
rownames(liver_counts_BvsA_df) <- resultsNames(dds_liver)[-1]  # Exclude the first category

# Transposing the data frame to have categories as headers
counts_table_liver_BvsA <- t(liver_counts_BvsA_df)

# Printing the table
print(counts_table_liver_BvsA)


# Saving the table to a CSV file
write.table(counts_table_liver_BvsA, file = "/home/brenda/BatHepatoTransc/intermediateData/liverDETS_table_by_category_BvsA.csv", 
            sep = ",", col.names = NA)



# For spleen DETs counts in condition B vs A


# Initializing variables to store counts
upregulated_spleen_counts_BvsA <- numeric(length(list_of_results_spleen) - 1)  # excluding the first result which is the baseMean
downregulated_spleen_counts_BvsA <- numeric(length(list_of_results_spleen) - 1)

# Iterating through each category except the first one
for (i in 2:length(list_of_results_spleen)) {
  result_spleen <- list_of_results_spleen[[i]]  # Gets results for the current category
  
  # Extracting significant rows based on adjusted p-value threshold
  significant_spleen_rows <- result_spleen[!is.na(result_spleen$padj) & result_spleen$padj < 0.01, ]
  
  # Counting upregulated and downregulated transcripts for condition B vs A
  upregulated_spleen_counts_BvsA[i - 1] <- sum(significant_spleen_rows$log2FoldChange < 0)  # Reverse the comparison
  downregulated_spleen_counts_BvsA[i - 1] <- sum(significant_spleen_rows$log2FoldChange > 0)  # Reverse the comparison
}

# Printing the counts
print(upregulated_spleen_counts_BvsA)
print(downregulated_spleen_counts_BvsA)

# Combining the counts vectors into a data frame
spleen_counts_BvsA_df <- data.frame(Upregulated = upregulated_spleen_counts_BvsA, Downregulated = downregulated_spleen_counts_BvsA)

# Setting the row names for the counts data frame
rownames(spleen_counts_BvsA_df) <- resultsNames(dds_spleen)[-1]  # Exclude the first category

# Transposing the data frame to have categories as headers
counts_table_spleen_BvsA <- t(spleen_counts_BvsA_df)

# Printing the table
print(counts_table_spleen_BvsA)


# Saving the table to a CSV file
write.table(counts_table_spleen_BvsA, file = "/home/brenda/BatHepatoTransc/intermediateData/spleenDETS_table_by_category_BvsA.csv", 
            sep = ",", col.names = NA)



# Summarizing the counts of differentially expressed transcripts (DETs) for condition B vs A

# DETs counts summary for spleen condition B vs A per category

# Initializing a vector to store the counts for each condition
spleen_DETs_summary_counts_BvsA <- numeric(length(category_names))

# Iterating through each category and count the number of DETs for condition B vs A
for (i in seq_along(category_names)) {
  spleen_DETs_summary_counts_BvsA[i] <- upregulated_spleen_counts_BvsA[i] + downregulated_spleen_counts_BvsA[i]
}

# Creating a summary data frame for spleen DETs counts for condition B vs A
spleen_DETs_counts_summary_BvsA_df <- data.frame(Category = category_names, Count = spleen_DETs_summary_counts_BvsA)

# Print the summary data frame
print(spleen_DETs_counts_summary_BvsA_df)


# Saving the table to a CSV file
write.table(spleen_DETs_counts_summary_BvsA_df, file = "/home/brenda/BatHepatoTransc/intermediateData/summary_counts_spleen_BvsA.csv", 
            sep = ",", col.names = NA)


# DETs counts summary for liver condition B vs A per category

# Initializing a vector to store the counts for each condition
liver_DETs_summary_counts_BvsA <- numeric(length(category_names))

# Iterating through each category and count the number of DETs for condition B vs A
for (i in seq_along(category_names)) {
  liver_DETs_summary_counts_BvsA[i] <- upregulated_liver_counts_BvsA[i] + downregulated_liver_counts_BvsA[i]
}

# Creating a summary data frame for liver DETs counts for condition B vs A
liver_DETs_counts_summary_BvsA_df <- data.frame(Category = category_names, Count = liver_DETs_summary_counts_BvsA)

# Print the summary data frame
print(liver_DETs_counts_summary_BvsA_df)


# Saving the table to a CSV file
write.table(liver_DETs_counts_summary_BvsA_df, file = "/home/brenda/BatHepatoTransc/intermediateData/summary_counts_liver_BvsA.csv", 
            sep = ",", col.names = NA)

