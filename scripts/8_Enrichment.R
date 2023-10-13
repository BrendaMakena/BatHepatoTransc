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
  toGO <-  new("topGOdata", ontology = ontology, allGenes = g, 
               annot = annFUN.gene2GO,
               gene2GO = annot)
  resultFis <- runTest(toGO, algorithm = "classic", statistic = "fisher")
  list(toGO, resultFis) ## returns a list first data then result
}

### now we can loop over the DE gene sets, we don't want to run the
### tests for the intercepts or the overall
# running the molecular function enrichment analysis
MF_enrichment <- lapply(DETs_ALL[!grepl("Intercept|overall", names(DETs_ALL))], 
                        function (mySet){
    TOGO.all.onto("MF", DETs_ALL[["overall"]], mySet, gene2GO)
})

## above I did it only for "MF", but we could also loop over "MF",
## "BP" and "CC" ;-)


### now a function to perform correction for multiple testing and to
### extract a table from the results
gene.table.topGO <- function(TOGO.list, fdr=0.1){
    all <- GenTable(TOGO.list[[1]], TOGO.list[[2]], topNodes=100)
    names(all)[names(all)%in%"result1"] <- "p.value"
    all$fdr <- p.adjust(all$p.value, method="BH")
    return(all[all$fdr<fdr,])
}



### and run it along the list of Enrichment results
MF_enrichment_tables <- lapply(MF_enrichment, gene.table.topGO)


lapply(MF_enrichment_tables, head)

MF_enrichment_tables["spleen:rpmh_scaled"]

MF_enrichment_tables["liver:rpmh_scaled"]

names(MF_enrichment_tables)



## let's do one test for the intersection of liver and spleen rpmh_scaled

TOGO.all.onto("MF", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:rpmh_scaled"]], ### overlaping genes (133)
                        DETs_ALL[["spleen:rpmh_scaled"]]),
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


## for the genes enriched in collagen binding MF, how many of the 133 DETs overlapping between
## spleen and liver rpmh scaled category, are in our DETs and also annotated as collagen  binding
table(inGO = allLocusGO$go_id%in%"GO:0005518", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
     DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:rpmh_scaled"]], ### overlaping genes (133)
                DETs_ALL[["spleen:rpmh_scaled"]]))

          ## 3 DETs
          ## this can be applied to other interesting enrichment clusters


### looping over the DE gene sets we don't want to run the
### tests for the intercepts or the overall
## running the biological process ontology analysis
BP_enrichment <- lapply(DETs_ALL[!grepl("Intercept|overall", names(DETs_ALL))], 
                        function (mySet){
  TOGO.all.onto("BP", DETs_ALL[["overall"]], mySet, gene2GO)
})


### now a function to perform correction for multiple testing and to
### extract a table from the results
gene.table.topGO <- function(TOGO.list, fdr=0.1){
  all <- GenTable(TOGO.list[[1]], TOGO.list[[2]], topNodes=100)
  names(all)[names(all)%in%"result1"] <- "p.value"
  all$fdr <- p.adjust(all$p.value, method="BH")
  return(all[all$fdr<fdr,])
}



### and run it along the list of Enrichment results
BP_enrichment_tables <- lapply(BP_enrichment, gene.table.topGO)

names(BP_enrichment_tables)

lapply(BP_enrichment_tables, head)

BP_enrichment_tables["spleen:rpmh_scaled"]

BP_enrichment_tables["liver:rpmh_scaled"]

## and for season category as QC step
BP_enrichment_tables["liver:Season_Rainy_vs_Dry"]

BP_enrichment_tables["spleen:Season_Rainy_vs_Dry"]

## doing one test for the intersection of liver and spleen rpmh_scaled

TOGO.all.onto("BP", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:rpmh_scaled"]], 
                        DETs_ALL[["spleen:rpmh_scaled"]]),
              gene2GO) %>% gene.table.topGO()


## sample output

## ##         GO.ID                                        Term Annotated Significant Expected p.value    fdr
## ##   GO:0071504                cellular response to heparin         4           2     0.04 0.00058 0.01548437
## ##   GO:0061448               connective tissue development       193           8     1.91 0.00065 0.01548437
## ##   GO:0071503                         response to heparin         5           2     0.05 0.00095 0.01548437
## ##   GO:0030198           extracellular matrix organization       207           8     2.05 0.00102 0.01548437
## ##   GO:0043062        extracellular structure organization       207           8     2.05 0.00102 0.01548437
## ##   GO:0045229 external encapsulating structure organiz...       208           8     2.06 0.00105 0.01548437
## ##   GO:0030199                collagen fibril organization        47           4     0.47 0.00117 0.01548437
## ##   GO:0010043                        response to zinc ion        25           3     0.25 0.00185 0.01548437
## ##   GO:0120178        steroid hormone biosynthetic process        25           3     0.25 0.00185 0.01548437
## ##  GO:2000586 regulation of platelet-derived growth fa...         7           2     0.07 0.00198 0.01548437
## ##  GO:0098758                   response to interleukin-8         1           1     0.01 0.00991 0.01548437
## ##  GO:0098759          cellular response to interleukin-8         1           1     0.01 0.00991 0.01548437
## ##  GO:0099502 calcium-dependent activation of synaptic...         1           1     0.01 0.00991 0.01548437
## ##  GO:1901634 positive regulation of synaptic vesicle ...         1           1     0.01 0.00991 0.01548437
## ##  GO:1901875 positive regulation of post-translationa...         1           1     0.01 0.00991 0.01548437
## ##  GO:2000607 negative regulation of cell proliferatio...         1           1     0.01 0.00991 0.01548437
## ##  GO:2000683    regulation of cellular response to X-ray         1           1     0.01 0.00991 0.01548437
## ##  GO:2000699 fibroblast growth factor receptor signal...         1           1     0.01 0.00991 0.01548437
## ## >


## for the genes enriched in cellular response to heparin BP, how many of the 133 DETs overlapping between
## spleen and liver rpmh scaled category, are in our DETs and also annotated as responsible for cellular 
## response to heparin
table(inGO = allLocusGO$go_id%in%"GO:0071504", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:rpmh_scaled"]], ### overlaping genes (133)
                  DETs_ALL[["spleen:rpmh_scaled"]]))

                      ## 2 DETs
                      ## this can be applied to other interesting enrichment clusters


### looping over the DE gene sets we don't want to run the
### tests for the intercepts or the overall
## running the cellular compartment ontology analysis
CC_enrichment <- lapply(DETs_ALL[!grepl("Intercept|overall", names(DETs_ALL))], 
                        function (mySet){
                          TOGO.all.onto("CC", DETs_ALL[["overall"]], mySet, gene2GO)
                        })


### now a function to perform correction for multiple testing and to
### extract a table from the results
gene.table.topGO <- function(TOGO.list, fdr=0.1){
  all <- GenTable(TOGO.list[[1]], TOGO.list[[2]], topNodes=100)
  names(all)[names(all)%in%"result1"] <- "p.value"
  all$fdr <- p.adjust(all$p.value, method="BH")
  return(all[all$fdr<fdr,])
}



### and running it along the list of Enrichment results
CC_enrichment_tables <- lapply(CC_enrichment, gene.table.topGO)


lapply(CC_enrichment_tables, head)

CC_enrichment_tables["spleen:rpmh_scaled"]

CC_enrichment_tables["liver:rpmh_scaled"]


## doing one test for the intersection of liver and spleen rpmh_scaled

TOGO.all.onto("CC", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:rpmh_scaled"]], ## 133 overlapping genes
                        DETs_ALL[["spleen:rpmh_scaled"]]),
              gene2GO) %>% gene.table.topGO()


## output

## ##       GO.ID                             Term Annotated Significant Expected p.value     fdr
## ## GO:0005614              interstitial matrix         6           2     0.06  0.0014 0.06666667
## ## GO:0005581                  collagen trimer        53           4     0.53  0.0019 0.06666667
## ## GO:0005883                    neurofilament         7           2     0.07  0.0020 0.06666667
## ## GO:0031012             extracellular matrix       318           9     3.16  0.0043 0.08800000
## ## GO:0030312 external encapsulating structure       319           9     3.17  0.0044 0.08800000
## ## GO:0005840                         ribosome       217           7     2.15  0.0058 0.09666667
## ## > 


## for the genes enriched in interstitial matrix CC, how many of the 133 DETs overlapping between
## spleen and liver rpmh scaled category, are in our DETs and also annotated as being in the cellular 
## compartment interstitial matrix
table(inGO = allLocusGO$go_id%in%"GO:0005614", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:rpmh_scaled"]], ### overlaping genes (133)
                  DETs_ALL[["spleen:rpmh_scaled"]]))

                ## 2 DETs
                ## this can be applied to other interesting enrichment clusters


