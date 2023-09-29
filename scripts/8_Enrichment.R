# Annotating the DETs for the host and liver

#install.packages("topGO")
library(biomaRt)
library(topGO)
library(org.Hs.eg.db)  # Loads the appropriate organism annotation package (e.g., human)
library(GO.db)

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
  DETs_ALL  <- readRDS("intermediateData/DETs_ALL.RDS")
}



# getting annotation database from ensembl (human genes)
#mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                dataset="hsapiens_gene_ensembl")




####### try 1

# Removing duplicates from entrezgene_accession
unique_entrezgene_accession <- unique(allLocusGO$entrezgene_accession)

# Finding the intersection of DETs_ALL and unique_entrezgene_accession
common_genes <- intersect(DETs_ALL$overall, unique_entrezgene_accession)

# Creating a numeric vector of detected genes where 1 indicates detection
detected_genes <- numeric(length(unique_entrezgene_accession))

# Marking the detected genes with 1
detected_genes[common_genes %in% DETs_ALL$overall] <- 1

# Creating the named vector with names as unique_entrezgene_accession
names(detected_genes) <- unique_entrezgene_accession

########## try 2

# Creating a named vector from the combined list of significant DETs
detected_genes_named <- as.character(DETs_ALL$overall)
names(detected_genes_named) <- detected_genes_named

# Creating a named vector with 1 indicating detection
detected_genes_named[detected_genes_named %in% 
                       rownames(allLocusGO$entrezgene_accession)] <- "1"


######## try 3

# Extracting all significant DETs from the 'overall' list in DETs_ALL
all_DETs <- DETs_ALL$overall

# Create a named vector with 1 indicating detection
allGenes <- rep(0, length(allLocusGO$entrezgene_accession))
names(allGenes) <- allLocusGO$entrezgene_accession

# Mark significant DETs as detected
allGenes[all_DETs %in% names(allGenes)] <- 1


str(allGenes)

# Defining a custom gene selection function
custom_gene_selection_fun <- function(allGenes) {
  return(allGenes == 1)
}


########### try 4 

# DETs_ALL is the list of significant genes
detected_genes <- rep(0, length(allLocusGO$entrezgene_accession))
detected_genes[names(allLocusGO$entrezgene_accession) %in% 
                 DETs_ALL$overall] <- 1
names(detected_genes) <- allLocusGO$entrezgene_accession

str(detected_genes)

# Creating the topGOdata object
GOdata <- new("topGOdata",
                    ontology = "BP",  # for biological process GO branch
                    allGenes = detected_genes,
                    nodeSize = 10,
                    annot = allLocusGO,
                    #annot = org.Hs.eg.db,  # Use the appropriate organism annotation package
                    #ID = "entrezgene_accession",  # Assuming "entrezgene" as the gene identifier
                    geneSelectionFun = function(k) k == 1)  # Custom gene selection function


str(DETs_ALL)


if ("topGO" %in% installed.packages()) {
  detach("package:topGO", unload = TRUE)
  remove.packages("topGO")
}

install.packages("topGO")


library(AnnotationDbi)
library(topGO)





