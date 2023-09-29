# Annotating the DETs for the host and liver

#install.packages("topGO")
library(biomaRt)
library(topGO)
library(org.Hs.eg.db)  # Loads the appropriate organism annotation package (e.g., human)


redoDE <- FALSE
redoAnnotation <- FALSE


if(redoAnnotation){
    source("scripts/6_annotation.R")
} else{
   allLocusGO <- readRDS("intermediateData/GOtermAnnot.RDS")
}

if(redoDE){
    source("scripts/DE_analysis.R")
} else {
    "intermediateData/DETs_ALL.RDS"
}



# getting annotation database from ensembl (human genes)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")



### loop to get GO enrichment terms for each of the 4 categories

# Creating a list to store the GO enrichment results for each condition
liver_GO_results_list <- list()

# Defining a function to select significant genes based on our criteria
significant_liver_genes <- function(condition_results) {
  # Extracting the row names (transcript IDs) from the condition_results
  transcript_ids <- rownames(condition_results)
  return(transcript_ids)
}

# Defining a custom gene selection function for topGO
custom_liver_gene_selection <- function(TopGOdata, node) {
  allGenes <- TopGOdata$allGenes
  # Filter genes based on our criteria here
  # For example, selecting significant genes by comparing with a threshold
  # Here, we'll select all genes as an example
  selected_genes <- allGenes
  return(selected_genes)
}

# Defining a custom annotation function for topGO
custom_annotation_fun_liver <- function(ontology, allGenes) {
  # Extracting GO terms for each gene and creating a named vector
  gene_GO_terms <- unique(ontology$GOdata$allLocusGO$GOID[ontology$GOdata$allLocusGO$entrezgene_accession %in% allGenes])
  return(gene_GO_terms)
}


# Looping through the list of results for each condition, excluding the intercept
for (i in 2:length(list_of_results)) {  # Starts from the second element
  # Extracting the results for the current condition
  condition_results <- list_of_results[[i]]
  
  # Using the gene selection function to get significant genes
  significant_liver_genes <- significant_liver_genes(condition_results)
  
  # Creating a named vector for allGenes using transcript_ids
  significant_liver_genes_named <- setNames(rep(1, length(significant_liver_genes)), 
                                   significant_liver_genes)
  
  # Subset the gene-GO mapping data for the current condition's DETs
  #allLocusGO_liver_subset <- allLocusGO[allLocusGO$entrezgene_accession %in% 
  #                            significant_liver_genes, ]
  
  # Creating a topGOdata object for the current condition
  liver_godata <- new("topGOdata", allLocusGO = allLocusGO, 
                ontology = "BP",   # for biological process GO branch
                allGenes = significant_liver_genes_named, nodeSize = 10,
                geneSelectionFun = custom_liver_gene_selection,
                annotationFun = custom_annotation_fun_liver)
                        # Using the custom gene selection function
                        # and annotation function
                        # allGenes is a vector containing all the gene IDs.
                        # nodeSize is the minimum number of genes required to 
                        # define a GO term as significant.
  
  # Performing GO enrichment analysis
  liver_GO_result <- runTest(liver_godata, algorithm = "classic", 
                    statistic = "fisher")
  
  # Storing the results in the list
  liver_GO_results_list[[i - 1]] <- liver_GO_result  # Use (i - 1) as an index to skip the intercept
}

# liver DETs GO and enrichment analysis
#liver DETs <- list_of_results  #DETs in the 4 categories
#liver significant DETs <- list_of_DETs_liver
# gene-GO mapping data <- allLocusGO

# names of our DETs (not significant DETs)
geneNames_liver <- featureNames(list_of_results_liver)
length(geneNames_liver)
    # 13317
str(geneNames_liver)

geneList_liver <- factor(as.integer(geneNames_liver %in% 
                                      list_of_DETs_liver))

# Ensuring there are two levels in the factor
levels(geneList_liver) <- c("0", "1")

names(geneList_liver) <- geneNames_liver  
str(geneList_liver) 

# creating the custom annotation function
custom_annotation_fun_liver <- function(GOdata, allGenes) {
  # Initializing an empty list to store gene-to-GO-term mappings
  gene_to_GO_terms_liver <- list()
  
  # Looping through each gene in allGenes
  for (gene in allGenes) {
    # Extracting the GO terms associated with the current gene
    terms <- GOdata$go_id[GOdata$entrezgene_accession == gene]
    
    # Checking if there are any GO terms for this gene
    if (length(terms) > 0) {
      gene_to_GO_terms_liver[[gene]] <- terms
    }
  }
  
  # Converting the list to a named list
  gene_to_GO_terms_liver <- lapply(names(gene_to_GO_terms_liver), 
                              function(gene) {
    setNames(gene_to_GO_terms_liver[[gene]], gene)
  })
  
  return(gene_to_GO_terms_liver)
}

# Creating a named vector with gene identifiers and binary indicators
allGenes_liver <- rep(0, length(geneNames_liver))
allGenes_liver[geneNames_liver %in% list_of_DETs_liver] <- 1
names(allGenes_liver) <- geneNames_liver
str(allGenes_liver)

# Defining a simple gene selection function
custom_gene_selection_fun <- function(TopGOdata, node) {
  allGenes <- TopGOdata$allGenes
  return(allGenes)
}

# Creating a topGOdata object for the current condition
liver_godata <- new("topGOdata", 
                    ontology = "BP",   # for biological process GO branch
                    allGenes = allGenes_liver, 
                    nodeSize = 10,
                    annotationFun = custom_annotation_fun_liver,
                    geneSelectionFun = custom_gene_selection_fun)
                    #annot = allLocusGO) # Specifying GO data frame
                     
                    
# Using the custom gene selection function
# and annotation function
# allGenes is a vector containing all the gene IDs.
# nodeSize is the minimum number of genes required to 
# define a GO term as significant.

# Performing GO enrichment analysis
liver_GO_result <- runTest(liver_godata, algorithm = "classic", 
                           statistic = "fisher")

# Storing the results in the list
liver_GO_results_list[[i - 1]] <- liver_GO_result  # Use (i - 1) as an index to skip the intercept


str(allLocusGO)


# Creating a custom annotation function
custom_annotation_fun_liver <- function(GOdata, allGenes) {
  # Initialize an empty list to store gene-to-GO-term mappings
  gene_to_GO_terms_liver <- list()
  
  # Loop through each gene in allGenes
  for (gene in names(allGenes)) {
    # Extract the GO terms associated with the current gene
    terms <- GOdata$go_id[GOdata$entrezgene_accession == gene]
    
    # Checking if there are any GO terms for this gene
    if (length(terms) > 0) {
      gene_to_GO_terms_liver[[gene]] <- terms
    }
  }
  
  # Converting the list to a named list
  gene_to_GO_terms_liver <- lapply(names(gene_to_GO_terms_liver), 
                             function(gene) {
    setNames(gene_to_GO_terms_liver[[gene]], gene)
  })
  
  return(gene_to_GO_terms_liver)
}

# Defining a simple gene selection function
custom_gene_selection_fun <- function(TopGOdata, node) {
  allGenes <- TopGOdata$allGenes
  return(allGenes)
}

# Creating a topGOdata object for the current condition
liver_GOdata <- new("topGOdata", 
                    ontology = "BP",   # for biological process GO branch
                    allGenes = allGenes_liver, 
                    nodeSize = 10,
                    annotationFun = custom_annotation_fun_liver,
                    #annot = allLocusGO,
                    geneSelectionFun = custom_gene_selection_fun)
library(topGO)
str(allLocusGO)
str(allGenes_liver)
head(allLocusGO)


# Load the topGO library if not already loaded
library(topGO)

# Define your custom annotation function
custom_annotation_fun_liver <- function(allLocusGO, allGenes) {
  # Initialize an empty list to store gene-to-GO-term mappings
  gene_to_GO_terms_liver <- list()
  
  # Loop through each gene in allGenes
  for (gene in names(allGenes)) {
    # Extract the GO terms associated with the current gene
    terms <- allLocusGO$go_id[allLocusGO$entrezgene_accession == gene]
    
    # Checking if there are any GO terms for this gene
    if (length(terms) > 0) {
      gene_to_GO_terms_liver[[gene]] <- terms
    }
  }
  
  # Converting the list to a named list
  gene_to_GO_terms_liver <- lapply(names(gene_to_GO_terms_liver), 
                                   function(gene) {
                                     setNames(gene_to_GO_terms_liver[[gene]], gene)
                                   })
  
  return(gene_to_GO_terms_liver)
}


custom_liver_gene_selection_fun <- function(allGenes) {
  return(allGenes)  # Returns all genes without any selection
}


# Creating a topGOdata object for the current condition
liver_GOdata <- new("topGOdata", 
                    ontology = "BP",   # for biological process GO branch
                    allGenes = allGenes_liver, 
                    nodeSize = 10,
                    annotationFun = custom_annotation_fun_liver,
                    geneSelectionFun = custom_liver_gene_selection_fun)
                    #annot = allLocusGO)

# Check the topGOdata object
liver_GOdata




# Defining a custom gene selection function
custom_gene_selection_fun <- function(allGenes, selectedGenes) {
  selectedGenes <- selectedGenes[selectedGenes %in% names(allGenes)]
  return(selectedGenes)
}


# Creating a vector of 0s and 1s where 1 indicates the gene is detected
detected_genes_liver <- as.numeric(names(list_of_results_liver) %in% list_of_DETs_liver)

# Assigning names to the vector using gene names or IDs
names(detected_genes_liver) <- names(list_of_results_liver)


str(detected_genes_liver)
# Creating the topGOdata object
liver_GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = detected_genes_liver,
                    geneSel = detected_genes_liver,
                    annotationFun = allLocusGO,
                    geneSelectionFun = custom_gene_selection_fun,
                    mapping = "org.Hs.eg.db"  # Use the appropriate organism annotation package
)

# Define a custom annotation function
custom_annotation_fun_liver <- function(allLocusGO, allGenes) {
  # Initialize an empty list to store gene-to-GO-term mappings
  gene_to_GO_terms_liver <- list()
  
  # Loop through each gene in allGenes
  for (gene in names(allGenes)) {
    # Extract the GO terms associated with the current gene
    terms <- allLocusGO$go_id[allLocusGO$entrezgene_accession == gene]
    
    # Checking if there are any GO terms for this gene
    if (length(terms) > 0) {
      gene_to_GO_terms_liver[[gene]] <- terms
    }
  }
  
  # Converting the list to a named list
  gene_to_GO_terms_liver <- lapply(names(gene_to_GO_terms_liver), 
                                   function(gene) {
                                     setNames(gene_to_GO_terms_liver[[gene]], gene)
                                   })
  
  return(gene_to_GO_terms_liver)
}

# Create a vector of 0s and 1s where 1 indicates the gene is detected
detected_genes_liver <- as.numeric(names(list_of_results_liver) %in% list_of_DETs_liver)

# Assign names to the vector using gene names or IDs
names(detected_genes_liver) <- names(list_of_results_liver)

# Create the topGOdata object with the custom annotation function
liver_GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = detected_genes_liver,
                    geneSel = detected_genes_liver,
                    annotationFun = custom_annotation_fun_liver,
                    geneSelectionFun = function(allGenes) allGenes,  # Use all genes without selection
                    mapping = "org.Hs.eg.db"  # Use the appropriate organism annotation package
)

# Load the topGO library if not already loaded
library(topGO)

# Create a named vector of detected genes where 1 indicates detection
detected_genes_liver <- ifelse(names(list_of_results_liver) %in% list_of_DETs_liver, 1, 0)

# Create the topGOdata object
liver_GOdata <- new("topGOdata",
                    ontology = "BP",  # for biological process GO branch
                    allGenes = detected_genes_liver,
                    nodeSize = 10,
                    annot = allLocusGO,
                    mapping = "org.Hs.eg.db",  # Use the appropriate organism annotation package
                    ID = "entrezgene_accession"  # Assuming "entrezgene" as the gene identifier
)

# Creating a named vector of detected genes
detected_genes_liver <- numeric(length(allLocusGO$entrezgene_accession))
detected_genes_liver[names(list_of_results_liver) %in% list_of_DETs_liver] <- 1
names(detected_genes_liver) <- allLocusGO$entrezgene_accession


# Defining a custom gene selection function
custom_gene_selection_fun <- function(allGenes) {
  return(allGenes)
}


# Creating the topGOdata object
liver_GOdata <- new("topGOdata",
                    ontology = "BP",  # for biological process GO branch
                    allGenes = detected_genes_liver,
                    nodeSize = 10,
                    #annot = allLocusGO,
                    annot = org.Hs.eg.db,  # Use the appropriate organism annotation package
                    ID = "entrezgene_accession",  # Assuming "entrezgene" as the gene identifier
                    geneSelectionFun = custom_gene_selection_fun
)


