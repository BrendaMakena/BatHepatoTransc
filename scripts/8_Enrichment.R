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


## let's do one test for the intersection of liver and spleen rpmh_scaled

TOGO.all.onto("MF", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:rpmh_scaled"]], DETs_ALL[["spleen:rpmh_scaled"]]),
              gene2GO) %>% gene.table.topGO()

## ##       GO.ID                                        Term Annotated Significant Expected p.value        fdr
## ## 1  GO:0070325       lipoprotein particle receptor binding        25           4     0.25 0.00011 0.01050000
## ## 2  GO:0005198                structural molecule activity       459          14     4.63 0.00021 0.01050000
## ## 3  GO:0050750 low-density lipoprotein particle recepto...        21           3     0.21 0.00117 0.03228571
## ## 4  GO:0005518                            collagen binding        50           4     0.50 0.00159 0.03228571
## ## 5  GO:0043394                        proteoglycan binding        24           3     0.24 0.00174 0.03228571
## ## 6  GO:0008307            structural constituent of muscle        26           3     0.26 0.00220 0.03228571
## ## 7  GO:0001540                        amyloid-beta binding        55           4     0.56 0.00226 0.03228571
## ## 8  GO:0004385                   guanylate kinase activity         8           2     0.08 0.00272 0.03400000
## ## 9  GO:0042277                             peptide binding       162           6     1.64 0.00592 0.04011538
## ## 10 GO:1901681                     sulfur compound binding       168           6     1.70 0.00704 0.04011538
## ## > 
