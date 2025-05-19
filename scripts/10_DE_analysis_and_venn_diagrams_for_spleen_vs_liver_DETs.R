## DETs analysis for host liver transcriptome from 3' Tag RNAseqs
## differential expression tests are based on a negative binomial
## generalized linear model


## Then plot venn diagrams of the liver vs spleen DETs and downstream enrichment analysis


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

# list of spleen transcripts with significant p value for all the conditions
list_of_spleen_DETs  <- lapply(list_of_results_spleen, function(rdf){
  rdf <- rdf[!is.na(rdf$padj),]
  rownames(rdf[rdf$padj< 0.1,])
})


## here the output of the analysis for the pipeline
DETs_ALL <- c(list_of_DETs, overall=list(rownames(list_of_results[[1]])))

DETs_liver <- c(list_of_liver_DETs, overall=list(rownames(list_of_results_liver[[1]])))

DETs_spleen <- c(list_of_spleen_DETs, overall=list(rownames(list_of_results_spleen[[1]])))


#saveRDS(DETs_ALL, "intermediateData/DETs_ALL.RDS")


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


#        liver
# spleen  FALSE  TRUE
# FALSE 11051  1866
# TRUE    267   133

## the 133 DETs overlapping between liver vs spleen rpmh scaled category are statistically significant
chisq.test(table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:rpmh_scaled"]],
                 liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:rpmh_scaled"]]))
                      
                      # X-squared = 106.06, df = 1, p-value < 2.2e-16


## overlap table of liver vs spleen season rainy vs dry categories
table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:Season_Rainy_vs_Dry"]],
      liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:Season_Rainy_vs_Dry"]])

#        liver
# spleen  FALSE  TRUE
#   FALSE 10556  1979
#   TRUE    531   251


## the 251 DETs overlapping between liver vs spleen Season_Rainy_vs_Dry category are statistically significant
chisq.test(table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:Season_Rainy_vs_Dry"]],
                 liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:Season_Rainy_vs_Dry"]]))

# X-squared = 139.27, df = 1, p-value < 2.2e-16


## overlap table of liver vs spleen age young vs adult categories
table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:Age_2category_Young_vs_Adult"]],
      liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:Age_2category_Young_vs_Adult"]])

#         liver
# spleen  FALSE TRUE
# FALSE  4847 1906
# TRUE   4427 2137

## the 2137 DETs overlapping between liver vs spleen age young vs adult category are statistically significant
chisq.test(table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:Age_2category_Young_vs_Adult"]],
                 liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:Age_2category_Young_vs_Adult"]]))

# X-squared = 29.338, df = 1, p-value = 6.078e-08


## overlap table of liver vs spleen sex male vs female categories
table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:Sex_Male_vs_Female"]],
      liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:Sex_Male_vs_Female"]])

#          liver
#  spleen  FALSE  TRUE
#  FALSE 12963   236
#  TRUE     70    48

## the 48 DETs overlapping between liver vs spleen sex male vs female category are statistically significant
chisq.test(table(spleen=rownames(list_of_results[[1]])%in%list_of_DETs[["spleen:Sex_Male_vs_Female"]],
                 liver=rownames(list_of_results[[1]])%in%list_of_DETs[["liver:Sex_Male_vs_Female"]]))

# Chi-squared approximation may be incorrect


# the expected counts are small
chisq.test(table(spleen=rownames(list_of_results[[1]]) %in% list_of_DETs[["spleen:Sex_Male_vs_Female"]],
                 liver=rownames(list_of_results[[1]]) %in% list_of_DETs[["liver:Sex_Male_vs_Female"]]))$expected

#             liver
# spleen       FALSE       TRUE
# FALSE 12917.5165 281.483517
# TRUE    115.4835   2.516483


# running Fisher's test instead because the expected counts are small (<5)
fisher.test(table(spleen=rownames(list_of_results[[1]]) %in% list_of_DETs[["spleen:Sex_Male_vs_Female"]],
                  liver=rownames(list_of_results[[1]]) %in% list_of_DETs[["liver:Sex_Male_vs_Female"]]))

# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval: 24.90346 56.45337
# sample estimates:   odds ratio 37.61018

# the p-value < 2.2e-16 is extremely small showing the DETs in the sex males vs females category
# in the liver and spleen are statistically significantly associated






#### visualising the DETs between the liver and spleen tissues ####


# Creating Venn diagram for the liver vs spleen infection intensity category
venn.plotliver_vs_spleen_rpmh_scaled <- venn.diagram(
  x = list(liver = list_of_liver_DETs[["liver:rpmh_scaled"]],
           spleen = list_of_spleen_DETs[["spleen:rpmh_scaled"]]),
  category.names = c("liver", "spleen"),
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = "blue",
  fill = "blue",
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of liver vs spleen DETs infection intensity category"
)

# Saving the Venn diagram as a PDF file
pdf(paste0("plots/Venn_liver_vs_spleen_infection_intensity.pdf"), width = 8, height = 8)
grid.draw(venn.plotliver_vs_spleen_rpmh_scaled)
dev.off()  # Close the PDF device



# Creating Venn diagram for the liver vs spleen season rainy vs dry category
venn.plotliver_vs_spleen_season_rainy_vs_dry <- venn.diagram(
  x = list(liver = list_of_liver_DETs[["liver:Season_Rainy_vs_Dry"]],
           spleen = list_of_spleen_DETs[["spleen:Season_Rainy_vs_Dry"]]),
  category.names = c("liver", "spleen"),
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = "orange",
  fill = "orange",
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of liver vs spleen DETs season rainy vs dry"
)

# Saving the Venn diagram as a PDF file
pdf(paste0("plots/Venn_liver_vs_spleen_season_rainy_vs_dry.pdf"), width = 8, height = 8)
grid.draw(venn.plotliver_vs_spleen_season_rainy_vs_dry)
dev.off()  # Close the PDF device



# Creating Venn diagram for the liver vs spleen age category young vs adult
venn.plotliver_vs_spleen_age_category_young_vs_adult <- venn.diagram(
  x = list(liver = list_of_liver_DETs[["liver:Age_2category_Young_vs_Adult"]],
           spleen = list_of_spleen_DETs[["spleen:Age_2category_Young_vs_Adult"]]),
  category.names = c("liver", "spleen"),
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = "purple",
  fill = "purple",
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of liver vs spleen DETs age category youngs vs adults"
)

# Saving the Venn diagram as a PDF file
pdf(paste0("plots/Venn_liver_vs_spleen_age_category_young_vs_adult.pdf"), width = 8, height = 8)
grid.draw(venn.plotliver_vs_spleen_age_category_young_vs_adult)
dev.off()  # Close the PDF device



# Creating Venn diagram for the liver vs spleen sex male vs female
venn.plotliver_vs_spleen_sex_male_vs_female <- venn.diagram(
  x = list(liver = list_of_liver_DETs[["liver:Sex_Male_vs_Female"]],
           spleen = list_of_spleen_DETs[["spleen:Sex_Male_vs_Female"]]),
  category.names = c("liver", "spleen"),
  filename = NULL,
  output = FALSE, # Prevents the creation of log files
  col = "green",
  fill = "green",
  annotation.cex = 1.2,  # Adjust the font size as needed
  main = "Venn diagram of liver vs spleen DETs sex category males vs females"
)

# Saving the Venn diagram as a PDF file
pdf(paste0("plots/Venn_liver_vs_spleen_sex_male_vs_female.pdf"), width = 8, height = 8)
grid.draw(venn.plotliver_vs_spleen_sex_male_vs_female)
dev.off()  # Close the PDF device

