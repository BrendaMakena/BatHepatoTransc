## DETs analysis for host liver transcriptome from 3' Tag RNAseqs
## differential expression tests are based on a negative binomial
## generalized linear model

## The script will ouput a list of DE genes (only the gene names for
## now), this could obviously be changed to output/transfer the whole
## of the results tables to the next script

## analysis for host liver and spleen transcripts  
library(DESeq2)
library(ggplot2)
library(dplyr)

library(VennDiagram)
library(ggVennDiagram)
library(gridExtra)
library(magrittr)


## We need the count data and the metadata for the DE analysis

redoCounting <- FALSE
redoMetadata <- FALSE

if(redoCounting){
    source("scripts/3_featurecounts.R")
}else{
    tagseqRNAfeatureCounts <- readRDS("intermediateData/countTable.RDS")
}

if(redoMetadata){
    source("scripts/4_metadata.R")
}else{
    metadata <- read.csv("intermediateData/metadata_expanded.csv")
}
rownames(metadata) <- metadata$ID
metadata$rpmh_scaled <- scale(metadata$rpmh)

# first filter the counts data to keep only host read counts 
#and > 500 counts across all samples
host_counts <-tagseqRNAfeatureCounts[!grepl("HEP_",rownames(tagseqRNAfeatureCounts)) &
                                     rowSums(tagseqRNAfeatureCounts)>500
                                    ,]


# Fltering metadata to separate rows for liver and spleen samples
liverIDs <- metadata$ID[metadata$Organ%in%"Liver"]
spleenIDs <- metadata$ID[metadata$Organ%in%"Spleen"]

# constructing the liver DESeqdataset object with rpmh_scaled as condition of test
dds_liver <- DESeqDataSetFromMatrix(countData = host_counts[,liverIDs],
                              colData = metadata[liverIDs,],
                              design = ~Season+Age_2category+Sex+
                                  rpmh_scaled,
                              tidy = FALSE)

# differential expression analysis (uses wald test statistics)
dds_liver <- DESeq(dds_liver)


# constructing the spleen DESeqdataset object with rpmh_scaled as condition of test
dds_spleen <- DESeqDataSetFromMatrix(countData = host_counts[,spleenIDs],
                              colData = metadata[spleenIDs,],
                              design = ~Season+Age_2category+Sex+
                                  rpmh_scaled,
                              tidy = FALSE)

# differential expression analysis (uses wald test statistics)
dds_spleen <- DESeq(dds_spleen)


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


# list of transcripts with significant p value for all the conditions
list_of_DETs  <- lapply(list_of_results, function(rdf){
    rdf <- rdf[!is.na(rdf$padj),]
    rownames(rdf[rdf$padj< 0.1,])
})

# list of liver transcripts with significant p value for all the conditions
list_of_liver_DETs  <- lapply(list_of_results_liver, function(rdf){
  rdf <- rdf[!is.na(rdf$padj),]
  rownames(rdf[rdf$padj< 0.1,])
})

# list of transcripts with significant p value for all the conditions
list_of_spleen_DETs  <- lapply(list_of_results_spleen, function(rdf){
  rdf <- rdf[!is.na(rdf$padj),]
  rownames(rdf[rdf$padj< 0.1,])
})


## here the output of the analysis for the pipeline
DETs_ALL <- c(list_of_DETs, overall=list(rownames(list_of_results[[1]])))

DETs_liver <- c(list_of_liver_DETs, overall=list(rownames(list_of_results_liver[[1]])))

DETs_spleen <- c(list_of_spleen_DETs, overall=list(rownames(list_of_results_spleen[[1]])))


saveRDS(DETs_ALL, "intermediateData/DETs_ALL.RDS")


### FROM HERE ONLY VISUALISATION (PLOTS) AND TABLES OUPUT
  
### Calculating overlaps between the 5 categories
## doing this to get a table for pairwise comparison of the categories
## could be a substitute for venn diagrams

getOverlapMatrix <- function(list_of_DETs){
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
    ## Create a summary table for the DE transcripts overlaps
    colnames(overlap_matrix) <- names(list_of_DETs)
    rownames(overlap_matrix) <- colnames(overlap_matrix)
    overlap_matrix
}

getOverlapMatrix(list_of_DETs[grep("liver", names(list_of_DETs))])

getOverlapMatrix(list_of_DETs[grep("spleen", names(list_of_DETs))])

getOverlapMatrix(list_of_DETs)

## overlap table of liver vs spleen rpmh scaled categories
table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:rpmh_scaled"]],
      liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:rpmh_scaled"]])

## the 133 DETs overlapping between liver vs spleen rpmh scaled category are statistically significant
chisq.test(table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:rpmh_scaled"]],
                 liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:rpmh_scaled"]]))
                      
                      # X-squared = 106.06, df = 1, p-value < 2.2e-16


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
  x = list_of_liver_DETs[-1],  # Exclude the first category ("Intercept")
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
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3, 
                              widths = unit(c(5, 1, 5), "null"))))

# Positioning the plots and the line in the grid
grid.draw(
  arrangeGrob(venn.plotliver, venn.plotspleen, ncol = 2)
)

# Drawing a vertical line to separate the plots
grid.lines(x = unit(0.5, "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2))


dev.off()  # Closing the PDF device


# venn diagram of eg only two of the 5 categories

# loop for all 5 categories venn diagrams compared pairwise

# Defining category names
category_names <- c("Season(Rainy vs Dry)", 
                    "Age 2category(Young vs Adult)", 
                    "Sex (Male vs Female)", "rpmh scaled")

# Defining colors for each category
category_colors <- c("yellow", "green", "red", "purple")


## loop for liver categories 

# Saving the Venn diagram as a PDF file
pdf("plots/Host_DETs_plots/Host_liver_DETs_plots/Venn_pairwise_categories_liver.pdf", 
    width = 8, height = 8)

## Creating liver Venn diagrams for all pairs of categories using all DETs
for (i in 1:(length(category_names))) {
  for (j in 1:length(category_names)) {
    cat1 <- category_names[i]
    cat2 <- category_names[j]
    
    # Select all DETs for the two categories
    selected_liver_DETs <- list(cat1 = list_of_liver_DETs[[i]],
                          cat2 = list_of_liver_DETs[[j]]) 
    
    # Creating Venn diagram for the pair of categories
    venn.plotliver <- venn.diagram(
      x = selected_liver_DETs,
      category.names = c(cat1, cat2),
      filename = NULL,
      output = FALSE,  # Prevents the creation of log files
      col = category_colors[category_names %in% c(cat1, cat2)],
      fill = category_colors[category_names %in% c(cat1, cat2)],
      annotation.cex = 1.2,
      main = paste("Venn Diagram liver:", cat1, "vs", cat2)
    )
    
    # Creating a new page in the PDF for each Venn diagram
    if (i != 1 || j != 1) {
      cat("PageBreak\n", file = "plots/Host_DETs_plots/Host_liver_DETs_plots/Venn_pairwise_categories_liver.pdf", append = TRUE)
    }
    
    # Drawing the Venn diagram on the PDF device
    grid.newpage()
    grid.draw(venn.plotliver)
    
    }
}

## closing the PDF and saving the combined Venn diagrams PDF
dev.off()

## loop for spleen categories

# Saving the Venn diagram as a PDF file
pdf("plots/Host_DETs_plots/Host_spleen_DETs_plots/Venn_pairwise_categories_spleen.pdf", 
    width = 8, height = 8)


## loop for spleen Venn diagrams for all pairs of categories using all DETs
for (i in 1:(length(category_names))) {
  for (j in 1:length(category_names)) {
    cat1 <- category_names[i]
    cat2 <- category_names[j]
    
    # Select all DETs for the two categories
    selected_spleen_DETs <- list(cat1 = list_of_spleen_DETs[[i]],
                          cat2 = list_of_spleen_DETs[[j]]) 
    
    # Creating Venn diagram for the pair of categories
    venn.plotspleen <- venn.diagram(
      x = selected_spleen_DETs,
      category.names = c(cat1, cat2),
      filename = NULL,
      output = FALSE,  # Prevents the creation of log files
      col = category_colors[category_names %in% c(cat1, cat2)],
      fill = category_colors[category_names %in% c(cat1, cat2)],
      annotation.cex = 1.2,
      main = paste("Venn Diagram spleen:", cat1, "vs", cat2)
    )
    
    # Creating a new page in the PDF for each Venn diagram
    if (i != 1 || j != 1) {
      cat("PageBreak\n", file = "plots/Host_DETs_plots/Host_spleen_DETs_plots/Venn_pairwise_categories_spleen.pdf", append = TRUE)
    }
    
    # Drawing the Venn diagram on the PDF device
    grid.newpage()
    grid.draw(venn.plotspleen)
  }
}

# Closing the PDF device to save the combined Venn diagrams PDF
dev.off()



# Creating Venn diagram for the liver vs spleen rpmh scaled category
venn.plotliver_vs_spleen_rpmh_scaled <- venn.diagram(
  x = list(liver = list_of_liver_DETs[["liver:rpmh_scaled"]],
           spleen = list_of_spleen_DETs[["spleen:rpmh_scaled"]]),
  category.names = c("liver", "spleen"),
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = "gold",
  fill = "gold",
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of liver vs spleen DETs rpmh scaled"
)

# Saving the Venn diagram as a PDF file
pdf(paste0("plots/Venn_liver_vs_spleen_rpmh_scaled.pdf"), width = 8, height = 8)
grid.draw(venn.plotliver_vs_spleen_rpmh_scaled)
dev.off()  # Close the PDF device




# MA plots
        #shows the log2 fold changes attributable to a given variable 
        #over the mean of normalized counts for all the samples in the
        #DESeqDataSet. Points are colored blue if the adjusted p value 
        #is less than 0.1 (statistically significant points). 
        #Points which fall out of the window are plotted as open 
        #triangles pointing either up or down.

## MA plots for liver DETs

# Opening a PDF device to save multiple MA plots
pdf("plots/Host_DETs_plots/Host_liver_DETs_plots/liver_MA_plots.pdf", width = 12, height = 12)

# Creating and saving MA plots for each category using a loop
for (i in 2:length(list_of_results_liver)) {
  category_name <- resultsNames(dds_liver)[i]
  plotMA(list_of_results_liver[[i]], 
      ylim = c(-2, 2), 
      main = paste("MA Plot for liver", category_name))
}

# Closing the PDF device
dev.off()

## MA plots for spleen DETs

# Opening a PDF device to save multiple MA plots
pdf("plots/Host_DETs_plots/Host_spleen_DETs_plots/spleen_MA_plots.pdf", width = 12, height = 12)

# Creating and saving MA plots for each category using a loop
for (i in 2:length(list_of_results_spleen)) {
  category_name <- resultsNames(dds_spleen)[i]
  plotMA(list_of_results_spleen[[i]], 
         ylim = c(-2, 2), 
         main = paste("MA Plot for spleen", category_name))
}

# Closing the PDF device
dev.off()


## getting the liver resLFC - Log fold change shrinkage for visualization and ranking

# name of coefficient to shrink for the categories
resultsNames(dds_liver)
resultsNames(dds_spleen)
       # It is more useful visualize the MA-plot for the shrunken log2 fold 
       # changes, which remove the noise associated with log2 fold changes 
       # from low count genes without requiring arbitrary filtering thresholds.

# plotting the LFC for all my 4 liver categories

# Mapping between category names and coefficient names
category_coefficients <- c(
  "Season(Rainy vs Dry)" = "Season_Rainy_vs_Dry",
  "Age 2category(Young vs Adult)" = "Age_2category_Young_vs_Adult",
  "Sex (Male vs Female)" = "Sex_Male_vs_Female",
  "rpmh scaled" = "rpmh_scaled"
)

## opening pdf file and setting path
pdf(file.path("plots/Host_DETs_plots/Host_liver_DETs_plots/All_liver_LFC_MA_plot.pdf"), 
    width = 12, height = 12)

# Looping through each category to create and save MA plots
for (category_name in names(category_coefficients)) {
  coef_name <- category_coefficients[category_name]
  
  # Check if the coefficient name exists in resultsNamesDDS
  if (coef_name %in% resultsNames(dds_liver)) {
    # Calculating resLFC for the given category
    res_liver_LFC <- lfcShrink(dds_liver, coef = coef_name, 
                     type = "apeglm")
    
    # Creating the MA plot
    plotMA(res_liver_LFC, ylim = c(-2, 2), 
           main = paste("Liver LFC MA Plot for category", category_name))
   } else {
    cat("Coefficient", coef_name, "not found in DESeq results.\n")
  }
}
  
## Closing the PDF device
dev.off()


## plotting the LFC for all my 4 spleen categories

## opening pdf file and setting path
pdf(file.path("plots/Host_DETs_plots/Host_spleen_DETs_plots/All_spleen_LFC_MA_plot.pdf"), 
    width = 12, height = 12)

# Looping through each category to create and save MA plots
for (category_name in names(category_coefficients)) {
  coef_name <- category_coefficients[category_name]
  
  # Check if the coefficient name exists in resultsNamesDDS
  if (coef_name %in% resultsNames(dds_spleen)) {
    # Calculating resLFC for the given category
    res_spleen_LFC <- lfcShrink(dds_spleen, coef = coef_name, 
                               type = "apeglm")
    
    # Creating the MA plot
    plotMA(res_spleen_LFC, ylim = c(-2, 2), 
           main = paste("Spleen LFC MA Plot for category", category_name))
  } else {
    cat("Coefficient", coef_name, "not found in DESeq results.\n")
  }
}

## Closing the PDF device
dev.off()

      ###using 'apeglm' for LFC shrinkage. If used in published research, please cite:
        #Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
        #sequence count data: removing the noise and preserving large differences.
        #Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895



### volcano plots

## liver volcano plots
# Opening PDF and defining path
pdf(file.path("plots/Host_DETs_plots/Host_liver_DETs_plots/liver_volcano_padjvalue<0.01.pdf"), 
    width = 12, height = 12)  

# volcano plots for all liver categories with padjvalue<0.01
for (i in 2:length(list_of_results_liver)) {
  category_name <- resultsNames(dds_liver)[i]  # Getting the category name
  result <- list_of_results_liver[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
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
    })
}

# Save the PDF file
dev.off()


# plotting the 4 categories with padjvalue of 0.1 instead of 0.01

# Opening PDF and defining path
pdf(file.path("plots/Host_DETs_plots/Host_liver_DETs_plots/liver_volcano_padjvalue<0.1.pdf"), 
    width = 12, height = 12)  


for (i in 2:length(list_of_results_liver)) {
  category_name <- resultsNames(dds_liver)[i]  # Getting the category name
  result <- list_of_results_liver[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
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
  })
}

# Closing the PDF
dev.off() 


## spleen volcano plots

# Opening PDF and defining path
pdf(file.path("plots/Host_DETs_plots/Host_spleen_DETs_plots/spleen_volcano_padjvalue<0.01.pdf"), 
    width = 12, height = 12)  

# volcano plots for all liver categories with padjvalue<0.01
for (i in 2:length(list_of_results_spleen)) {
  category_name <- resultsNames(dds_spleen)[i]  # Getting the category name
  result <- list_of_results_spleen[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
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
  })
}

# Save the PDF file
dev.off()


# plotting the 4 categories with padjvalue of 0.1 instead of 0.01

# Opening PDF and defining path
pdf(file.path("plots/Host_DETs_plots/Host_spleen_DETs_plots/spleen_volcano_padjvalue<0.1.pdf"), 
    width = 12, height = 12)  


for (i in 2:length(list_of_results_spleen)) {
  category_name <- resultsNames(dds_spleen)[i]  # Getting the category name
  result <- list_of_results_spleen[[i]]  # Get results for the current category
  
  # Create the volcano plot
  with(result, {
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
  })
}

# Closing the PDF
dev.off() 



### PCA plots

#First transform the raw count data
#vst function performs variance stabilizing transformation


## PCAs for the 4 liver categories

## Opening PDF and setting path
pdf(file.path("plots/Host_DETs_plots/Host_liver_DETs_plots/Liver_PCA_plots.pdf"), 
              width = 12, height = 12)

## looping through the categories
for (category_name in c("Season", "Age_2category", "Sex", "rpmh_scaled")) {
  # Performing variance stabilizing transformation (VST) for the current category
  vs_data_liver <- vst(dds_liver, blind = TRUE)
  
  # Creating the PCA plot
  
  # Plot PCA using ggplot2
  pca_data_liver <- plotPCA(vs_data_liver, intgroup = category_name, returnData = TRUE)
  
  # Create the ggplot2 PCA plot and print it
  pca_plot_liver <- ggplot(pca_data_liver, 
    aes(x = PC1, y = PC2, color = as.factor(get(category_name)))) +
    geom_point() +
    labs(title = paste("PCA plot for liver", category_name), color = category_name) # sets legend name
  
  print(pca_plot_liver)
  
}

# Closing the main PDF file
dev.off()


## PCAs for the 4 spleen categories

## Opening PDF and setting path
pdf(file.path("plots/Host_DETs_plots/Host_spleen_DETs_plots/Spleen_PCA_plots.pdf"), 
    width = 12, height = 12)

## looping through the categories
for (category_name in c("Season", "Age_2category", "Sex", "rpmh_scaled")) {
  # Performing variance stabilizing transformation (VST) for the current category
  vs_data_spleen <- vst(dds_spleen, blind = TRUE)
  
  # Creating the PCA plot
  
  # Plot PCA using ggplot2
  pca_data_spleen <- plotPCA(vs_data_spleen, intgroup = category_name, returnData = TRUE)
  
  # Create the ggplot2 PCA plot and print it
  pca_plot_spleen <- ggplot(pca_data_spleen, 
                           aes(x = PC1, y = PC2, color = as.factor(get(category_name)))) +
    geom_point() +
    labs(title = paste("PCA plot for spleen", category_name), color = category_name) # sets legend name
  
  print(pca_plot_spleen)
  
}

# Closing the main PDF file
dev.off()
