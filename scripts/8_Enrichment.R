# Annotating the DETs for the host and liver

#install.packages("topGO")
library(biomaRt)
library(topGO)
library(org.Hs.eg.db)  # Loads the appropriate organism annotation package (e.g., human)
library(GO.db)
library(magrittr)

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

#### we need a list associating each gene ID with its (potentially
#### multiple) GO terms. We transform the GO data frame to such a list
gene2GO <- by(allLocusGO, allLocusGO$entrezgene_accession,
              function(x) as.character(x$go_id))


### This function will compute the enrichment tests for the GO
### ontology ("MF" = molecular function, "BP" = biological process and
### "CC" = cellular compartment. allgenes is always our complete set
### of tested genes and annot always our annotation list
TOGO.all.onto <- function (ontology, allgenes, gene.set, annot) {
  g <- factor(as.integer( allgenes%in%gene.set ))
  names(g) <- allgenes
  toGO <-  new("topGOdata", ontology = ontology, allGenes = g, annot = annFUN.gene2GO,
               gene2GO = annot)
  resultFis <- runTest(toGO, algorithm = "classic", statistic = "fisher")
  list(toGO, resultFis) ## returns a list first data then result
}

### now we can loop over the DE gene sets we don't want to run the
### tests for the intercepts or the
MF_enrichment <- lapply(DETs_ALL[!grepl("Intercept|overall", names(DETs_ALL))], function (mySet){
    TOGO.all.onto("MF", DETs_ALL[["overall"]], mySet, gene2GO)
})

## above I did it only for "MF", but we could also loop over "MF",
## "BP" and "CC" ;-)


### now a function to perform correction for multiple testing and to
### extract a table from the resulst
gene.table.topGO <- function(TOGO.list, fdr=0.1){
    all <- GenTable(TOGO.list[[1]], TOGO.list[[2]], topNodes=100)
    names(all)[names(all)%in%"result1"] <- "p.value"
    all$fdr <- p.adjust(all$p.value, method="BH")
    return(all[all$fdr<pval,])
}


### and run it along the list of Enrichment results
MF_enrichment_tables <- lapply(MF_enrichment, gene.table.topGO)


lapply(MF_enrichment_tables, head)

MF_enrichment_tables["spleen:rpmh_scaled"]

MF_enrichment_tables["liver:rpmh_scaled"]



