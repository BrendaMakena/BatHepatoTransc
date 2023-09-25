# DETs analysis for host liver transcriptome from 3' Tag RNAseqs
# differential expression tests are based on a negative binomial 
#generalized linear model

# analysis for host liver transcripts  
BiocManager::install("apeglm")
install.packages("VennDiagram")
install.packages("ggVennDiagram")
library(DESeq2)
library(ggplot2)
library(dplyr)
library(apeglm)
library(VennDiagram)
library(ggVennDiagram)
library(gridExtra)


# first filter the counts data to keep only host read counts 
#and > 500 counts across all samples
host_counts <-tagseqRNAfeatureCounts %>% 
              filter(!grepl("HEP_",rownames(tagseqRNAfeatureCounts),
              ignore.case = TRUE),
              rowSums(tagseqRNAfeatureCounts)>500)

# Filtering metadata to keep only rows for liver samples
metadata_liver <- metadata %>%
                  filter(Organ == 'Liver')

# separating host liver from host spleen counts
# liver host counts
host_counts_liver <- host_counts[, metadata_liver$ID]


    #### DETs analysis for liver ####

# constructing the DESeqdataset object with rpmh_scaled as condition of test
dds_liver <- DESeqDataSetFromMatrix(countData = host_counts_liver,
                              colData = metadata_liver,
                              design = ~Season+Age_2category+Sex+
                                rpmh_scaled,
                              tidy = FALSE)

# viewing the object
dds_liver

# differential expression analysis (uses wald test statistics)
dds_liver <- DESeq(dds_liver)
             # does everything from normalization to linear modeling
                
                #estimateSizeFactors
                #This calculates the relative library depth of each sample 
                
                #estimateDispersions
                #estimates the dispersion of counts for each gene 
                
                #nbinomWaldTest
                #calculates the significance of coefficients in a 
                #Negative Binomial GLM using the size and dispersion outputs

# getting the results table
list_of_results <- lapply(resultsNames(dds_liver), function(n){
      results(dds_liver, name = n)
  })

lapply(list_of_results, head)

# list of transcripts with significant p value for all the conditions
list_of_DETs <- lapply(list_of_results, function(rdf){
          rdf <- rdf[!is.na(rdf$padj),]
         rownames(rdf[rdf$padj< 0.1,])
           })

# number of DETs for the categories
lapply(list_of_DETs, length)
lapply(list_of_DETs, head)


  
          #### plotting the results ####  


# Calculating overlaps between the 5 categories
overlap_matrix <- matrix(0, nrow = length(list_of_DETs), 
                  ncol = length(list_of_DETs))
for (i in 1:length(list_of_DETs)) {
  for (j in 1:length(list_of_DETs)) {
    if (i != j) {
      common_genes <- length(intersect(list_of_DETs[[i]], 
                      list_of_DETs[[j]]))
      overlap_matrix[i, j] <- common_genes
    }
  }
}

# Create a summary table for the DE transcripts overlaps
colnames(overlap_matrix) <- c("Intercept", "Season", "Age_2category", "Sex", "rpmh_scaled")
rownames(overlap_matrix) <- colnames(overlap_matrix)
overlap_matrix


# Creating Venn diagrams for all liver and spleen DETs in each category
# side by side

# Defining category names (excluding "Intercept")
category_names <- c("Season(Rainy vs Dry)", 
                    "Age 2category(Young vs Adult)", 
                    "Sex (Male vs Female)", "rpmh scaled")

# Defining colors for each category
category_colors <- c("yellow", "green", "red", "purple")


# Creating Venn diagram for the liver categories
venn.plotliver <- venn.diagram(
  x = list_of_DETs[-1],  # Exclude the first category ("Intercept")
  category.names = category_names,
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = category_colors,
  fill = category_colors,
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of all liver DETs in all categories"
)

# Creating Venn diagram for the spleen categories
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
pdf("plots/Venn_all_DETs_all_categories_liver_and_spleen_with_vertical_grid_line.pdf", width = 16, height = 8)

# Creating a layout with two columns
layout_matrix <- rbind(
  c(1, 2, 3),  # 1st column, plot 1, 2nd column
  c(1, 2, 3)   # 1st column, line, 2nd column
)

# Creating a grid with the line and the plots
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3, widths = unit(c(5, 1, 5), "null"))))

# Positioning the plots and the line in the grid
grid.draw(
  arrangeGrob(venn.plotliver, venn.plotspleen, ncol = 2)
)

# Drawing a vertical line to separate the plots
grid.lines(x = unit(0.5, "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2))


dev.off()  # Closing the PDF device


# venn diagrams for the top 250 DETs in each of the five categories

# Extracting the top 250 DETs from each category
top_250_DETs_list <- lapply(list_of_DETs, 
                     function(det_list) head(det_list, 250))

# Creating Venn diagrams
venn.plot <- venn.diagram(
  x = top_250_DETs_list,
  category.names = c("Intercept", "Season(Rainy vs Dry) ", 
                     "Age 2category(Young vs Adult)", 
                     "Sex (Male vs Female)", "rpmh scaled"),
  filename = NULL,  
  output = TRUE,
  col = c("dodgerblue", "yellow", "green", "red", "purple"), # Specifies the edge colors
  fill = c("dodgerblue", "yellow", "green", "red", "purple"), # Specifies the fill colors
  annotation.cex = 1.2,
  main = "Top 250 liver DETs"
)

# Saving the Venn diagram as a pdf file
pdf("plots/Top_250_DETs_Venn_diagram_liver.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()  # Close the pdf file


# venn diagram of eg only two of the 5 categories
# after line 115 choose the categories to include in the Venn diagram
selected_categories <- c("Season(Rainy vs Dry)", "rpmh scaled")

# Creating Venn diagrams for the selected categories
selected_DETs <- top_250_DETs_list[c("Intercept", "Season(Rainy vs Dry)", 
                                     "Age 2category(Young vs Adult)", 
                                     "Sex (Male vs Female)", "rpmh scaled") %in% selected_categories]

venn.plot <- venn.diagram(
  x = selected_DETs,
  category.names = selected_categories,
  filename = NULL,
  output = TRUE,
  col = c("dodgerblue", "yellow", "green", "red", "purple")[c("Intercept", "Season(Rainy vs Dry)", 
                                                              "Age 2category(Young vs Adult)", 
                                                              "Sex (Male vs Female)", "rpmh scaled") %in% selected_categories],
  fill = c("dodgerblue", "yellow", "green", "red", "purple")[c("Intercept", "Season(Rainy vs Dry)", 
                                                               "Age 2category(Young vs Adult)", 
                                                               "Sex (Male vs Female)", "rpmh scaled") %in% selected_categories],
  annotation.cex = 1.2,
  main = "liver season vs rpmh scaled"
)

# Saving the Venn diagram as a PDF file
pdf("plots/Season_vs_rpmh_scaled_Venn_diagram_liver.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()  # Close the PDF file


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
    selected_DETs <- list_of_DETs[category_names %in% c(cat1, cat2)]
    
    # Creating Venn diagram for the pair of categories
    venn.plot <- venn.diagram(
      x = selected_DETs,
      category.names = c(cat1, cat2),
      filename = NULL,
      output = FALSE,  # Prevents the creation of log files
      col = category_colors[category_names %in% c(cat1, cat2)],
      fill = category_colors[category_names %in% c(cat1, cat2)],
      annotation.cex = 1.2  # Adjust the font size as needed
    )
    
    # Saving the Venn diagram as a PDF file
    pdf(paste0("plots/Venn_", cat1, "_vs_", cat2, "_liver.pdf"), width = 8, height = 8)  # Adjust width and height as needed
    grid.draw(venn.plot)
    dev.off()  # Close the PDF device
  }
}



# MA plots
        #shows the log2 fold changes attributable to a given variable 
        #over the mean of normalized counts for all the samples in the
        #DESeqDataSet. Points are colored blue if the adjusted p value 
        #is less than 0.1 (statistically significant points). 
        #Points which fall out of the window are plotted as open 
        #triangles pointing either up or down.


# Function to create and save MA plots in DESeq2 style
create_MA_plot_liver <- function(result, category_name,plots) {
  pdf(file.path("plots/", paste(category_name, "_liver_MA_plot.pdf")))
  plotMA(result, ylim = c(-2, 2), main = paste("MA Plot for", category_name))
  dev.off()
}

# Opening a PDF device to save multiple MA plots
pdf("plots/Host_DETs_plots/Host_liver_DETs_plots/liver_MA_plots.pdf", width = 12, height = 12)


# Creating and saving MA plots for each category using a loop
for (i in 2:length(list_of_results)) {
  category_name <- resultsNames(dds_liver)[i]
  plotMA(list_of_results[[i]], 
      ylim = c(-2, 2), 
      main = paste("MA Plot for", category_name))
}

# Closing the PDF device
dev.off()

# getting the resLFC - Log fold change shrinkage for visualization and ranking

# name of coefficient to shrink for only one category
resultsNames(dds_liver)

res_liver_LFC <- lfcShrink(dds_liver, coef = "Sex_Male_vs_Female",
                           type = "apeglm")
res_liver_LFC

       # It is more useful visualize the MA-plot for the shrunken log2 fold 
       # changes, which remove the noise associated with log2 fold changes 
       # from low count genes without requiring arbitrary filtering thresholds.

# plotting the LFC for the one category
pdf("plots/liver_LFC_DTEs_sex_MA_plot.pdf", width = 12, height = 12)
plotMA(res_liver_LFC, ylim = c(-2, 2))
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
  if (coef_name %in% resultsNames(dds_liver)) {
    # Calculating resLFC for the given category
    res_liver_LFC <- lfcShrink(dds_liver, coef = coef_name, 
                     type = "apeglm")
    
    # Creating the MA plot
    pdf(file.path("plots/", paste("liver_LFC_", 
                                  category_name, "_MA_plot.pdf")), 
        width = 12, height = 12)
    plotMA(res_liver_LFC, ylim = c(-2, 2), 
           main = paste("Liver LFC MA Plot for category", category_name))
    dev.off()
  } else {
    cat("Coefficient", coef_name, "not found in DESeq results.\n")
  }
}


# Plot counts
        # plotCounts function compares the normalized counts between the comparison groups,
        # sex Male vs female, for the genes

# plotting transcript with lowest padjusted value from results table
pdf("plots/liver_DTEs_plotcounts_for_min_padjvalue_transcript.pdf", 
    width = 12, height = 12)

plotCounts(dds_liver, gene=which.min(res_liver$padj), 
           intgroup="Sex")

dev.off()

#alternatively plotting using ggplot2
pdf("plots/liver_DTEs_plotcounts_for_min_padjvalue2_transcript.pdf", width = 12, height = 12)

ggplot(plotCounts(dds_liver, gene=which.min(res_liver$padj), 
                  intgroup="Sex", returnData = TRUE), 
       aes(x=Sex, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

dev.off()

# plot counts for my 5 categories in a loop
# Looping through each category and creating and saving ggplot2-based count plots
for (category_name in c("Season", "Age_2category", "Sex", "rpmh_scaled")) {
  # Create a list to store data for all genes in the category
  gene_data_list <- list()
  
  # Loop through genes and create ggplot2 plots for each
  all_genes <- rownames(dds_liver)
  for (gene in all_genes) {
    # Subset the data for the specific gene and create the ggplot2 plot
    gene_data <- plotCounts(dds_liver, gene = gene, 
                 intgroup =  category_name, returnData = TRUE)
    
    gene_data_list[[gene]] <- gene_data
  }
  
  # Save ggplot2-based count plots for all genes in the category
  pdf(file.path("plots/", paste("liver_DTEs_plotcounts_for_", category_name, "_transcripts.pdf")), 
      width = 12, height = 12)
  
  # Plot all genes in the category
  for (gene in all_genes) {
    ggplot(gene_data_list[[gene]], aes(x = .data[[category_name]], y = count)) + 
      geom_point(position = position_jitter(w = 0.1, h = 0)) + 
      scale_y_log10(breaks = c(25, 100, 400))
  }
  
  dev.off()
}


# Iterate through category names
for (category_name in c("Season", "Age_2category", "Sex", "rpmh_scaled")) {
  # Create a list to store data for all genes in the category
  gene_data_list <- list()
  
  # Loop through DETs and create ggplot2 plots for each
  for (gene in list_of_DETs[[i]]) {
    # Subset the data for the specific DET and create the ggplot2 plot
    gene_data <- plotCounts(dds_liver, gene = gene, 
                            intgroup = category_name, returnData = TRUE)
    
    gene_data_list[[gene]] <- gene_data
  }
  
  # Save ggplot2-based count plots for DETs in the category
  pdf(file.path("plots/", paste("liver_DETs_plotcounts_for_", category_name, "_transcripts.pdf")), 
      width = 12, height = 12)
  
  # Plot all DETs in the category
  for (gene in list_of_DETs[[i]]) {
    ggplot(gene_data_list[[gene]], aes(x = category_name, y = count)) + 
      geom_point(position = position_jitter(w = 0.1, h = 0)) + 
      scale_y_log10(breaks = c(25, 100, 400))
  }
  
  dev.off()
}

# Iterate through category names
for (category_name in c("Season", "Age_2category", "Sex", "rpmh_scaled")) {
  # Create a list to store data for all genes in the category
  gene_data_list <- list()
  
  # Loop through DETs and create ggplot2 plots for each
  for (gene in list_of_DETs[[category_name]]) {
    # Subset the data for the specific DET and create the ggplot2 plot
    gene_data <- plotCounts(dds_liver, gene = gene, 
                            intgroup = category_name, returnData = TRUE)
    
    gene_data_list[[gene]] <- gene_data
  }
  
  # Save ggplot2-based count plots for DETs in the category
  pdf(file.path("plots/", paste("liver_DETs_plot_counts_for_", category_name, "_transcripts.pdf")), 
      width = 12, height = 12)
  
  # Plot all DETs in the category
  for (gene in list_of_DETs[[category_name]]) {
    p <- ggplot(gene_data_list[[gene]], aes(x = .data[[category_name]], y = count)) + 
      geom_point(position = position_jitter(w = 0.1, h = 0)) + 
      scale_y_log10(breaks = c(25, 100, 400))
    
    print(p)  # Print the plot to generate and save it
  }
  
  dev.off()
}


# alternatively plotting the counts for individual transcripts
pdf("plots/liver_DTEs_plotcounts_for_LOC107508184_transcript.pdf", width = 12, height = 12)

plotCounts(dds_liver, gene = "LOC107508184", intgroup = "Sex" )

dev.off()

# getting information on variables and tests used 
mcols(res_liver)$description


# volcano plots

# volcano plots for all categories with padjvalue<0.01

for (i in seq_along(list_of_results)) {
  category_name <- resultsNames(dds_liver)[i]  # Getting the category name
  result <- list_of_results[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
    pdf(file.path("plots/", paste("liver_volcano_", category_name, "_padjvalue<0.01.pdf")), 
        width = 12, height = 12)  # Opens PDF for the current category
    
    plot(log2FoldChange, -log10(padj), 
         pch = 20, main = paste("Volcano plot for liver", category_name), 
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

# Save the PDF file
dev.off()


# plotting the 5 categories with padjvalue of 0.1 instead of 0.01
for (i in seq_along(list_of_results)) {
  category_name <- resultsNames(dds_liver)[i]  # Getting the category name
  result <- list_of_results[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
    pdf(file.path("plots/", paste("liver_volcano_", category_name, "_padjvalue<0.1.pdf")), 
        width = 12, height = 12)  # Opens PDF for the current category
    
    plot(log2FoldChange, -log10(padj), 
         pch = 20, main = paste("Volcano plot for liver", category_name), 
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
  vs_data_liver <- vst(dds_liver, blind = TRUE)
  
  # Create the PCA plot
  pdf_file_name <- paste("plots/", "Liver_", category_name, "_PCA_plot.pdf", sep = "")
  pdf(pdf_file_name, width = 12, height = 12)
  
  # Plot PCA using ggplot2
  pca_data_liver <- plotPCA(vs_data_liver, intgroup = category_name, returnData = TRUE)
  
  # Create the ggplot2 PCA plot and print it
  pca_plot_liver <- ggplot(pca_data_liver, 
    aes(x = PC1, y = PC2, color = as.factor(get(category_name)))) +
    geom_point() +
    labs(title = paste("PCA plot for liver", category_name), color = category_name) # sets legend name
  
  print(pca_plot_liver)
  
  # Close the current PDF file
  dev.off()
}

# Close the main PDF file
dev.off()
