# DETs analysis for host spleen transcriptome from 3' Tag RNAseqs
# differential expression tests are based on a negative binomial 
#generalized linear model

# analysis for host spleen transcripts
BiocManager::install("apeglm")
install.packages("VennDiagram")
install.packages("ggVennDiagram")
library(DESeq2)
library(ggplot2)
library(dplyr)
library(apeglm)
library(VennDiagram)
library(ggVennDiagram)


# Filtering metadata to keep only rows for spleen samples
metadata_spleen <- metadata %>%
                   filter(Organ == 'Spleen')

# separating host liver from host spleen counts
# spleen host counts
host_counts_spleen <- host_counts[, metadata_spleen$ID]


####  DETs analysis for spleen ####

# constructing the DESeqdataset object with rpmh_scaled as condition of test
dds_spleen <- DESeqDataSetFromMatrix(countData = host_counts_spleen,
                              colData = metadata_spleen,
                              design = ~Season+Age_2category+Sex+
                              rpmh_scaled,
                              tidy = FALSE)

# viewing the object
dds_spleen

# differential expression analysis (uses wald test statistics)
dds_spleen <- DESeq(dds_spleen)
        # does everything from normalization to linear modeling

        #estimateSizeFactors
        #This calculates the relative library depth of each sample 
        
        #estimateDispersions
        #estimates the dispersion of counts for each gene 
        
        #nbinomWaldTest
        #calculates the significance of coefficients in a 
        #Negative Binomial GLM using the size and dispersion outputs


# getting the results table
#res_spleen <- results(dds_spleen)
#res_spleen

#alternatively
list_of_spleen_results <- lapply(resultsNames(dds_spleen), 
                          function(n){
                          results(dds_spleen, name = n)
                          })

lapply(list_of_spleen_results, head)

# head(results(dds_spleen, tidy = TRUE))

# summary of DTE
#summary(res_spleen)

# sorting table by p values
#res_spleen <- res_spleen[order(res_spleen$padj),]
#head(res_spleen)


# list of transcripts with significant p value for all the conditions
list_of_spleen_DETs <- lapply(list_of_spleen_results, 
                      function(rdf){
                      rdf <- rdf[!is.na(rdf$padj),]
                      rownames(rdf[rdf$padj< 0.1,])
                      })

# number of DETs for the categories
lapply(list_of_spleen_DETs, length)
lapply(list_of_spleen_DETs, head)



          #### plotting the results #### 

# Calculating overlaps between the 5 categories
overlap_matrix_spleen <- matrix(0, nrow = length(list_of_spleen_DETs), 
                         ncol = length(list_of_spleen_DETs))
for (i in 1:length(list_of_spleen_DETs)) {
  for (j in 1:length(list_of_spleen_DETs)) {
    if (i != j) {
      common_genes <- length(intersect(list_of_spleen_DETs[[i]], 
                                       list_of_spleen_DETs[[j]]))
      overlap_matrix[i, j] <- common_genes
    }
  }
}

# Create a summary table for the DE transcripts overlaps
colnames(overlap_matrix_spleen) <- c("Intercept", "Season", "Age_2category", "Sex", "rpmh_scaled")
rownames(overlap_matrix_spleen) <- colnames(overlap_matrix_spleen)
overlap_matrix_spleen


# venn diagram for all the DETs in the 4 categories

# Defining category names 
category_names <- c("Season(Rainy vs Dry)", 
                    "Age 2category(Young vs Adult)", 
                    "Sex (Male vs Female)", "rpmh scaled")

# Defining colors for each category
category_colors <- c("yellow", "green", "red", "purple") 


# Creating Venn diagrams for all transcripts in each category
for (i in 1:length(category_names)) {
  cat <- category_names[i]
  
  # Select all DETs for the current category
  selected_spleen_DETs <- list_of_spleen_DETs[[i]]
  
  # Create Venn diagram for the category
  venn.plotspleen <- venn.diagram(
    x = selected_spleen_DETs,
    category.names = cat,
    filename = NULL,
    output = FALSE, # Prevents the creation of log files
    col = category_colors[i],
    fill = category_colors[i],
    annotation.cex = 1.2  
  )
  
  # Saving the Venn diagram as a PDF file
  pdf(paste0("plots/Venn_all_categories_spleen.pdf"), width = 4, height = 4)  
  grid.draw(venn.plots)
  dev.off()  # Close the PDF device
}

# Define category names (excluding "Intercept")
category_names <- c("Season(Rainy vs Dry)", 
                    "Age 2category(Young vs Adult)", 
                    "Sex (Male vs Female)", "rpmh scaled")

# Define colors for each category
category_colors <- c("yellow", "green", "red", "purple")


# Create Venn diagram for the other four categories
venn.plotspleen <- venn.diagram(
  x = list_of_spleen_DETs[-1],  # Exclude the first category ("Intercept")
  category.names = category_names,
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = category_colors,
  fill = category_colors,
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of all spleen DETs in all categories"
)

# Saving the Venn diagram as a PDF file
pdf("plots/Venn_all_DETs_all_categories_spleen.pdf", width = 8, height = 8)  # Adjust width and height as needed
grid.draw(venn.plotspleen)
dev.off()  # Close the PDF device


# loop for all 5 categories venn diagrams compared pairwise

# Defining category names
category_names <- c("Intercept", "Season(Rainy vs Dry)", 
                    "Age 2category(Young vs Adult)", 
                    "Sex (Male vs Female)", "rpmh scaled")

# Defining colors for each category
category_colors <- c("dodgerblue", "yellow", "green", "red", "purple")

# Creating Venn diagrams for all pairs of categories using all DETs
for (i in 1:(length(category_names) - 1)) {
  for (j in (i + 1):length(category_names)) {
    cat1 <- category_names[i]
    cat2 <- category_names[j]
    
    # Select all DETs for the two categories
    selected_spleen_DETs <- list_of_spleen_DETs[category_names %in% c(cat1, cat2)]
    
    # Creating Venn diagram for the pair of categories
    venn.plotspleen <- venn.diagram(
      x = selected_spleen_DETs,
      category.names = c(cat1, cat2),
      filename = NULL,
      output = FALSE,  # Prevents the creation of log files
      col = category_colors[category_names %in% c(cat1, cat2)],
      fill = category_colors[category_names %in% c(cat1, cat2)],
      annotation.cex = 1.2
    )
    
    # Saving the Venn diagram as a PDF file
    pdf(paste0("plots/Venn_", cat1, "_vs_", cat2, "_spleen.pdf"), width = 8, height = 8)
    grid.draw(venn.plotspleen)
    dev.off()  # Close the PDF device
  }
}



#MA plots
        #shows the log2 fold changes attributable to a given variable 
        #over the mean of normalized counts for all the samples in the
        #DESeqDataSet. Points are colored blue if the adjusted p value 
        #is less than 0.1 (statistically significant points). 
        #Points which fall out of the window are plotted as open 
        #triangles pointing either up or down.

# Function to create and save MA plots in DESeq2 style
create_MA_plot_spleen <- function(result, category_name,plots) {
  pdf(file.path("plots/", paste(category_name, "_spleen_MA_plot.pdf")))
  plotMA(result, ylim = c(-2, 2), main = paste("MA Plot for", category_name))
  dev.off()
}

# Opening a PDF device to save multiple MA plots
pdf("plots/Host_DETs_plots/Host_spleen_DETs_plots/spleen_MA_plots.pdf", width = 12, height = 12)

# Create and save MA plots for each category using a loop
for (i in 2:length(list_of_spleen_results)) { #starts from second category
  category_name <- resultsNames(dds_spleen)[i]
  plotMA(list_of_spleen_results[[i]], 
         ylim = c(-2, 2), 
         main = paste("MA Plot for", category_name))
}

# Closing the PDF device
dev.off()

# getting the resLFC - Log fold change shrinkage for visualization and ranking

# name of coefficient to shrink
resultsNames(dds_spleen)

res_spleen_LFC <- lfcShrink(dds_spleen, coef = "Sex_Male_vs_Female",
                           type = "apeglm")

res_spleen_LFC

      # It is more useful visualize the MA-plot for the shrunken log2 fold 
      # changes, which remove the noise associated with log2 fold changes 
      # from low count genes without requiring arbitrary filtering thresholds.

# plotting the LFC
pdf("plots/spleen_LFC_DETs_sex_MA_plot.pdf", width = 12, height = 12)
plotMA(res_spleen_LFC, ylim = c(-2, 2))
dev.off()

# plotting the LFC for all my 4 categories

# Mapping between category names and coefficient names
category_coefficients <- c(
  "Season(Rainy vs Dry)" = "Season_Rainy_vs_Dry",
  "Age 2category(Young vs Adult)" = "Age_2category_Young_vs_Adult",
  "Sex (Male vs Female)" = "Sex_Male_vs_Female",
  "rpmh scaled" = "rpmh_scaled"
)

# Looping through each category to create and save MA plots
for (category_name in names(category_coefficients)) {
  coef_name <- category_coefficients[category_name]
  
  # Check if the coefficient name exists in resultsNamesDDS
  if (coef_name %in% resultsNames(dds_spleen)) {
    # Calculating resLFC for the given category
    res_spleen_LFC <- lfcShrink(dds_spleen, coef = coef_name, 
                               type = "apeglm")
    
    # Creating the MA plot
    pdf(file.path("plots/", paste("spleen_LFC_", 
                                  category_name, "_MA_plot.pdf")), 
        width = 12, height = 12)
    plotMA(res_spleen_LFC, ylim = c(-2, 2), 
           main = paste("Spleen LFC MA Plot for category", category_name))
    dev.off()
  } else {
    cat("Coefficient", coef_name, "not found in DESeq results.\n")
  }
}


# Plot counts

  # plotCounts function compares the normalized counts between the 
  # comparison groups, sex Male vs female, for the genes

# plotting transcript with lowest p adjusted value from results table
pdf("plots/spleen_DTEs_plotcounts_for_min_padjvalue_transcript.pdf", width = 12, height = 12)
plotCounts(dds_spleen, gene=which.min(res_spleen$padj), 
           intgroup="Sex")

dev.off()

#alternatively plotting using ggplot2
pdf("plots/spleen_DTEs_plotcounts_for_min_padjvalue2_transcript.pdf", width = 12, height = 12)

ggplot(plotCounts(dds_spleen, gene=which.min(res_spleen$padj), 
                  intgroup="Sex", returnData = TRUE), 
       aes(x=Sex, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

dev.off()

# alternatively plotting the counts for individual transcripts
pdf("plots/spleen_DETs_plotcounts_for_LOC107508184_transcript_category_sex.pdf", width = 12, height = 12)

plotCounts(dds_spleen, gene = "LOC107508184", intgroup = "Sex" )

dev.off()

# getting information on variables and tests used 
mcols(res_spleen)$description


# volcano plots

# Looping through each category to create and save volcano plots with padj value < 0.01
for (i in seq_along(list_of_spleen_results)) {
  category_name <- resultsNames(dds_spleen)[i]  # Getting the category name
  result <- list_of_spleen_results[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
    pdf(file.path("plots/", paste("spleen_volcano_", category_name, "_padjvalue<0.01.pdf")), 
        width = 12, height = 12)  # Opens PDF for the current category
    
    plot(log2FoldChange, -log10(padj), 
         pch = 20, main = paste("Volcano plot for spleen", category_name), 
         xlim = c(-3, 3))
    
    # Adding colored points: blue if padj < 0.01, red if log2FC > 1 and padj < 0.05)
    points(log2FoldChange, -log10(padj), pch = 20, 
           col = ifelse(padj < 0.01, "blue", 
                        ifelse(abs(log2FoldChange) > 1 & padj < 0.05, "red", "black")))
    
    # Customize the legend
    legend("topright", legend = c("padj < 0.01", "log2FC > 1 & padj < 0.05", "Other"), 
           col = c("blue", "red", "black"), pch = 20)
    
    dev.off()  # Close the PDF for the current category
  })
}

# Looping through each category to create and save volcano plots with padj value < 0.1
for (i in seq_along(list_of_spleen_results)) {
  category_name <- resultsNames(dds_spleen)[i]  # Getting the category name
  result <- list_of_spleen_results[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
    pdf(file.path("plots/", paste("spleen_volcano_", category_name, "_padjvalue<0.1.pdf")), 
        width = 12, height = 12)  # Opens PDF for the current category
    
    plot(log2FoldChange, -log10(padj), 
         pch = 20, main = paste("Volcano plot for spleen", category_name), 
         xlim = c(-3, 3))
    
    # Adding colored points: blue if padj < 0.1, red if log2FC > 1 and padj < 0.05)
    points(log2FoldChange, -log10(padj), pch = 20, 
           col = ifelse(padj < 0.1, "blue", 
                        ifelse(abs(log2FoldChange) > 1 & padj < 0.05, "red", "black")))
    
    # Customize the legend
    legend("topright", legend = c("padj < 0.1", "log2FC > 1 & padj < 0.05", "Other"), 
           col = c("blue", "red", "black"), pch = 20)
    
    dev.off()  # Close the PDF for the current category
  })
}


# PCA plots
#First transform the raw count data
#vst function performs variance stabilizing transformation


# PCAs for the 4 categories
for (category_name in c("Season", "Age_2category", "Sex", "rpmh_scaled")) {
  # Performing variance stabilizing transformation (VST) for the current category
  vs_data_spleen <- vst(dds_spleen, blind = TRUE)
  
  # Create the PCA plot
  pdf_file_name <- paste("plots/", "Spleen_", category_name, "_PCA_plot.pdf", sep = "")
  pdf(pdf_file_name, width = 12, height = 12)
  
  # Plot PCA using ggplot2
  pca_data_spleen <- plotPCA(vs_data_spleen, intgroup = category_name, returnData = TRUE)
  
  # Create the ggplot2 PCA plot and print it
  pca_plot_spleen <- ggplot(pca_data_spleen, 
    aes(x = PC1, y = PC2, color = as.factor(get(category_name)))) +
    geom_point() +
    labs(title = paste("PCA plot for spleen", category_name), color = category_name) # sets legend name
  
  print(pca_plot_spleen)
  
  # Close the current PDF file
  dev.off()
}

# Close the main PDF file
dev.off()
