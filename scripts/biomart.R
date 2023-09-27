# Annotating the DETs for the host and liver

library(biomaRt)
install.packages("topGO")
library(topGO)


# getting annotation database from ensembl (human genes)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl")


## where did this table acually come from MAKE A NOTE, so that we are
## able to cite this!!!
proteins <- read.csv("inputdata/rousettus_proteins.csv")

    # table from ncbi

# filtering necessary data (columns) from entrez file
#getBM(mart=mart,
#      attributes = c("entrezgene_accession", "entrezgene_id",
 #                    "go_id"),
  #    filters =  "entrezgene_accession",                       
   #   values = "PPP6R3") ### value in "Locus" testing one

### This works only for the short abreviated "gene names", not for the
### accessions starting with LOCxxxxxxx (x are numbers)


### But let's see whether this is a problem and how big...

### Before you use this you might want to subset your "proteins" to
### only those really tested in the DE analysis!

allLocusGO <- getBM(mart=mart,
                    attributes = c("entrezgene_accession",
                                   "entrezgene_id",
                                   "entrezgene_description",
                                   "go_id"),
                    filters =  "entrezgene_accession",
                    values = proteins$Locus) ###all values in "Locus"

table(unique(proteins$Locus)%in%allLocusGO$entrezgene_accession)
    ## FALSE  TRUE 
    ##  3540 15529 

## we miss (only?) 3540 genes in the overall set
## how is this for the expression-analysed genes?

table(found=unique(proteins$Locus)%in%allLocusGO$entrezgene_accession,
      LOC=grepl("LOC", unique(proteins$Locus)))

                    # LOC
        # found   FALSE  TRUE
        # FALSE   353  3187   # 15520 protein identifiers not LOC have
        # TRUE  15520     9   # annotations, 9 don't have 

## even better: almost all the ones we're are really those with "LOC"

## how is this for the expression-analysed genes?

# for our >500 read counts host transcripts, how many are annotated 
table(inTransCounts=unique(proteins$Locus)%in%rownames(host_counts), 
  inAnnotation=unique(proteins$Locus)%in%allLocusGO$entrezgene_accession)

                    # inAnnotation
    # inTransCounts FALSE  TRUE
            # FALSE  2197  4175
            # TRUE   1343 11354

# chi square test to show statistical significance for the annotated transcripts
chisq.test(table(inTransCounts=unique(proteins$Locus)%in%rownames(host_counts), 
          inAnnotation=unique(proteins$Locus)%in%allLocusGO$entrezgene_accession))
  
  # X-squared = 1601.7, df = 1, p-value < 2.2e-16

# table of transcripts ids in unfiltered host feature counts file 
# that are in the protein annotation file
table(rownames(tagseqRNAfeatureCounts)%in%proteins$Locus)

    # FALSE  TRUE 
    # 4835 18611 

# table of annotated unfiltered host transcripts in our DETs 
table(inTransCounts=proteins$Locus%in%rownames(tagseqRNAfeatureCounts), 
      inAnnotation=proteins$Locus%in%allLocusGO$entrezgene_accession)

              # inAnnotation
# inTransCounts FALSE  TRUE
        # FALSE   317   192
        # TRUE   6089 49793

# table of unique annotated transcripts in our DETs
table(inTransCounts=unique(proteins$Locus)%in%rownames(tagseqRNAfeatureCounts), 
      inAnnotation=unique(proteins$Locus)%in%allLocusGO$entrezgene_accession)
 
                    # inAnnotation
      # inTransCounts FALSE  TRUE
              # FALSE   299   159
              # TRUE   3241 15370

# table of unique annotated transcripts in our DETs 
#from the >500 counts host transcripts
table(inTransCounts=unique(proteins$Locus)%in%rownames(host_counts), 
      inAnnotation=unique(proteins$Locus)%in%allLocusGO$entrezgene_accession)

                  # inAnnotation
    # inTransCounts FALSE  TRUE
           # FALSE  2197  4175
           # TRUE   1343 11354



#### annotation of the 11354 transcripts with identifiers using TopGO #### 

# separate analysis for liver and spleen DETs

# liver DETs GO and enrichment analysis
#liver DETs <- list_of_results  #DETs in the 4 categories
# gene-GO mapping data <- allLocusGO

# subsetting the gene-GO mapping data to only contain transcripts in 
# the DETs list
#allLocusGO_subset <- allLocusGO[allLocusGO[, 1] %in% list_of_results, ]

# the transcript ids for the individual categories are the rownames
rownames(list_of_results$Season_Rainy_vs_Dry)
head(list_of_results$Season_Rainy_vs_Dry)

?topGO

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







