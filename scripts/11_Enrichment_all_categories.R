# Annotating the DETs for the host and liver in all the 4 categories

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
    source("scripts/7_DE_analysis.R")
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
## "BP" and "CC" 


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



## let's do one test for the intersection of liver and spleen infection intensity

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


## listing the names of the 3 collagen binding DETs that are in 
##the 133 overlapping liver and spleen DETs in rpmH scaled category

cat("collagen binding DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0005518" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:rpmh_scaled", "spleen:rpmh_scaled")])
          ]
        ), 
        decreasing = TRUE
      ), 3
    )), "\n")


## collagen binding DETs: COL14A1, HSD17B12 and LOX genes

## protein IDs on NCBI: NP_066933.1, NP_057226.1 and NP_002308.2 
              


##Doing the test for the intersection of liver and spleen season category

TOGO.all.onto("MF", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Season_Rainy_vs_Dry"]], ### overlaping genes (251)
                        DETs_ALL[["spleen:Season_Rainy_vs_Dry"]]),
              gene2GO) %>% gene.table.topGO()



#     GO.ID                                        Term       Annotated   Significant Expected p.value  fdr
# 1   GO:0016491                     oxidoreductase activity       463          24     8.59   5e-06 0.000500000
# 2   GO:0004857                   enzyme inhibitor activity       255          14     4.73 0.00028 0.006800000
# 3   GO:0055102                   lipase inhibitor activity         8           3     0.15 0.00033 0.006800000
# 4   GO:0047747                 cholate-CoA ligase activity         2           2     0.04 0.00034 0.006800000
# 5   GO:0070653 high-density lipoprotein particle recept...         2           2     0.04 0.00034 0.006800000
# 6   GO:0016878                  acid-thiol ligase activity        20           4     0.37 0.00044 0.007333333
# 7   GO:0016616 oxidoreductase activity, acting on the C...        78           7     1.45 0.00059 0.007466667
# 8   GO:0033293                 monocarboxylic acid binding        39           5     0.72 0.00072 0.007466667
# 9   GO:0031406                     carboxylic acid binding       108           8     2.00 0.00088 0.007466667
# 10  GO:0016614 oxidoreductase activity, acting on CH-OH...        84           7     1.56 0.00092 0.007466667



## for the genes enriched in oxidoreductase activity MF, how many of the 251 DETs overlapping between
## spleen and liver season category, are in our DETs and also annotated as oxidoreductase activity
table(inGO = allLocusGO$go_id%in%"GO:0016491", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Season_Rainy_vs_Dry"]], ### overlaping genes (251)
                  DETs_ALL[["spleen:Season_Rainy_vs_Dry"]]))

## 19 DETs



## listing the names of the 19 collagen binding DETs that are in 
##the 251 overlapping liver and spleen DETs in season category

cat("oxidoreductase activity DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0016491" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Season_Rainy_vs_Dry", "spleen:Season_Rainy_vs_Dry")])
          ]
        ), 
        decreasing = TRUE
      ), 19
    )), "\n")

## oxidoreductase activity DETs: ACAD11 ACOX2 AKR1A1 ALDH1L1 CAT GFOD2 GRHPR HGD HPD HSD17B10 
##                               JMJD1C ME3 MIOX MOXD1 PAH SELENOM SLC27A5 SORD XDH genes

## protein IDs on NCBI: NP_115545.3  NP_003491.1 NP_001189342.1 NP_001257293.1  NP_110446.3 




## The test for the intersection of liver and spleen age category

TOGO.all.onto("MF", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Age_2category_Young_vs_Adult"]], ### overlaping genes (2137)
                        DETs_ALL[["spleen:Age_2category_Young_vs_Adult"]]),
              gene2GO) %>% gene.table.topGO()



#       GO.ID                                        Term       Annotated   Significant Expected p.value  fdr
#   1   GO:0004879                   nuclear receptor activity        39          19     6.63 4.7e-06 0.000235000
#   2   GO:0098531 ligand-activated transcription factor ac...        39          19     6.63 4.7e-06 0.000235000
#   3   GO:0003723                                 RNA binding      1387         289   235.67 3.8e-05 0.001266667
#   4   GO:0044877          protein-containing complex binding      1303         269   221.39 0.00014 0.003000000
#   5   GO:0001221           transcription coregulator binding        96          31    16.31 0.00017 0.003000000
#   6   GO:0003676                        nucleic acid binding      2760         531   468.95 0.00018 0.003000000
#   7   GO:0003682                           chromatin binding       452         105    76.80 0.00032 0.004500000
#   8   GO:0046965         nuclear retinoid X receptor binding        16           9     2.72 0.00042 0.004500000
#   9   GO:0140657                      ATP-dependent activity       437         101    74.25 0.00050 0.004500000
#   10  GO:0140662     ATP-dependent protein folding chaperone        26          12     4.42 0.00051 0.004500000
#   11  GO:0016887                     ATP hydrolysis activity       293          72    49.78 0.00052 0.004500000
#   12  GO:0050839              cell adhesion molecule binding       433         100    73.57 0.00054 0.004500000
#   13  GO:0097159             organic cyclic compound binding      4181         773   710.40 0.00060 0.004615385



## for the genes enriched in RNA binding MF, how many of the 2137 DETs overlapping between
## spleen and liver age category, are in our DETs and also annotated as RNA binding
table(inGO = allLocusGO$go_id%in%"GO:0003723", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Age_2category_Young_vs_Adult"]], ### overlaping genes (2137)
                  DETs_ALL[["spleen:Age_2category_Young_vs_Adult"]]))

## 276 DETs



## listing the names of some of the RNA activity DETs that are in 
##the 2137 overlapping liver and spleen DETs in age category

cat("RNA binding DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0003723" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Age_2category_Young_vs_Adult", 
                                                                         "spleen:Age_2category_Young_vs_Adult")])
          ]
        ), 
        decreasing = TRUE
      ), 20
    )), "\n")

## RNA binding DETs: A1CF AARS1 AARS2 ACO1 ADAD1 ADAR ADD1 AGO1 AGO3 AKAP1 ALDH18A1 ALDOA ALKBH5 ANP32A 
#                    ARID5A BCLAF1 BUD23 BZW1 C1D C4BPA 

## protein IDs on NCBI: NP_001185747.1  NP_001596.2  NP_065796.2 NP_001265281.1  NP_001152757.1



## The test for the intersection of liver and spleen sex category

TOGO.all.onto("MF", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Sex_Male_vs_Female"]], ### overlaping genes (48)
                        DETs_ALL[["spleen:Sex_Male_vs_Female"]]),
              gene2GO) %>% gene.table.topGO()


#         GO.ID                                 Term            Annotated   Significant Expected p.value  fdr
#   1   GO:0016705 oxidoreductase activity, acting on paire...        91           4     0.27 0.00015 0.01300000
#   2   GO:0016491                     oxidoreductase activity       463           7     1.39 0.00037 0.01300000
#   3   GO:0016706 2-oxoglutarate-dependent dioxygenase act...        48           3     0.14 0.00039 0.01300000
#   4   GO:0051213                        dioxygenase activity        70           3     0.21 0.00118 0.02950000
#   5   GO:0141052             histone H3 demethylase activity        21           2     0.06 0.00177 0.03000000
#   6   GO:0032452                histone demethylase activity        24           2     0.07 0.00231 0.03000000
#   7   GO:0140457                protein demethylase activity        24           2     0.07 0.00231 0.03000000
#   8   GO:0050251                  retinol isomerase activity         1           1     0.00 0.00300 0.03000000
#   9   GO:0051908 double-stranded DNA 5'-3' DNA exonucleas...         1           1     0.00 0.00300 0.03000000
#   10  GO:0070080                      titin Z domain binding         1           1     0.00 0.00300 0.03000000
#   11  GO:0032451                        demethylase activity        31           2     0.09 0.00385 0.03266667
#   12  GO:0004175                      endopeptidase activity       218           4     0.65 0.00392 0.03266667


## for the genes enriched in dioxygenase activity MF, how many of the 48 DETs overlapping between
## spleen and liver sex category, are in our DETs and also annotated as dioxygenase activity
table(inGO = allLocusGO$go_id%in%"GO:0051213", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Sex_Male_vs_Female"]], ### overlaping genes (48)
                  DETs_ALL[["spleen:Sex_Male_vs_Female"]]))

## 3 DETs



## listing the names of some of the dioxygenase activity DETs that are in 
##the 48 overlapping liver and spleen DETs in sex category

cat("dioxygenase activity DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0051213" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Sex_Male_vs_Female", 
                                                                         "spleen:Sex_Male_vs_Female")])
          ]
        ), 
        decreasing = TRUE
      ), 3
    )), "\n")

## dioxygenase activity DETs: EGLN1 KDM5C KDM6A
## protein IDs on NCBI: NP_001364189.1  NP_001140174.1  NP_001278344.1




#### Enrichment analysis for the DETs for the BP category #### 


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



## doing the test for the intersection of liver and spleen infection intensity category

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
## spleen and liver infection intensity category, are in our DETs and also annotated as responsible for cellular 
## response to heparin
table(inGO = allLocusGO$go_id%in%"GO:0071504", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:rpmh_scaled"]], ### overlaping genes (133)
                  DETs_ALL[["spleen:rpmh_scaled"]]))

                      ## 2 DETs
                      


## listing the names of the 2 cellular response to heparin DETs that are in 
##the 133 overlapping liver and spleen DETs in infection intensity category

cat("cellular response to heparin DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0071504" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:rpmh_scaled", 
                                                                         "spleen:rpmh_scaled")])
          ]
        ), 
        decreasing = TRUE
      ), 2
    )), "\n")

## cellular response to heparin DETs: EGR1 SLIT2 
## protein IDs on NCBI: NP_001955.1 NP_001276064.1 




## doing the test for the intersection of liver and spleen season category

TOGO.all.onto("BP", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Season_Rainy_vs_Dry"]], 
                        DETs_ALL[["spleen:Season_Rainy_vs_Dry"]]),
              gene2GO) %>% gene.table.topGO()


#           GO.ID                                 Term          Annotated   Significant Expected p.value    fdr
#   1   GO:0044282            small molecule catabolic process       264          25     4.91 1.8e-11 1.125000e-09
#   2   GO:0044281            small molecule metabolic process      1272          58    23.66 3.0e-11 1.125000e-09
#   3   GO:0016054              organic acid catabolic process       190          21     3.53 4.5e-11 1.125000e-09
#   4   GO:0046395           carboxylic acid catabolic process       190          21     3.53 4.5e-11 1.125000e-09
#   5   GO:0019752           carboxylic acid metabolic process       624          37    11.61 2.3e-10 4.600000e-09
#   6   GO:0043436                   oxoacid metabolic process       638          37    11.87 4.3e-10 7.166667e-09
#   7   GO:0006082              organic acid metabolic process       642          37    11.94 5.1e-10 7.285714e-09
#   8   GO:0072329       monocarboxylic acid catabolic process       105          15     1.95 8.3e-10 1.037500e-08
#   9   GO:0032787       monocarboxylic acid metabolic process       424          27     7.89 1.9e-08 1.900000e-07
#   10  GO:0044283         small molecule biosynthetic process       424          27     7.89 1.9e-08 1.900000e-07
#   11  GO:0046394        carboxylic acid biosynthetic process       223          19     4.15 3.1e-08 2.818182e-07
#   12  GO:0016053           organic acid biosynthetic process       226          19     4.20 3.9e-08 3.250000e-07
#   13  GO:1901615 organic hydroxy compound metabolic proce...       346          23     6.44 1.1e-07 8.461538e-07
#   14  GO:0044248                  cellular catabolic process      1222          47    22.73 7.4e-07 5.285714e-06
#   15  GO:0006699              bile acid biosynthetic process        22           6     0.41 2.2e-06 1.466667e-05



## for the genes enriched in bile acid biosynthetic process BP, how many of the 251 DETs overlapping between
## spleen and liver season category, are in our DETs and also annotated as responsible for  
## bile acid biosynthetic process
table(inGO = allLocusGO$go_id%in%"GO:0006699", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Season_Rainy_vs_Dry"]], ### overlaping genes (251)
                  DETs_ALL[["spleen:Season_Rainy_vs_Dry"]]))

## 6 DETs


## listing the names of the 6 bile acid biosynthetic process DETs that are in 
##the 251 overlapping liver and spleen DETs in season category

cat("bile acid biosynthetic process DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0006699" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Season_Rainy_vs_Dry", 
                                                                         "spleen:Season_Rainy_vs_Dry")])
          ]
        ), 
        decreasing = TRUE
      ), 6
    )), "\n")

## bile acid biosynthetic process DETs: ABCB11 ACOX2 HSD17B10 SCP2 SLC27A2 SLC27A5 
## protein IDs on NCBI: NP_003733.2  NP_003491.1 NP_001032900.1 NP_001007099.1 NP_001153101.1 NP_001308125.1



## doing the test for the intersection of liver and spleen age category BP DETs

TOGO.all.onto("BP", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Age_2category_Young_vs_Adult"]], 
                        DETs_ALL[["spleen:Age_2category_Young_vs_Adult"]]),
              gene2GO) %>% gene.table.topGO()



#           GO.ID                                 Term          Annotated   Significant Expected p.value    fdr
#   1   GO:0051276                     chromosome organization       498         135    85.25 7.3e-09 0.0000007300
#   2   GO:0140014                    mitotic nuclear division       225          68    38.52 7.1e-07 0.0000355000
#   3   GO:1903047                  mitotic cell cycle process       621         150   106.30 2.7e-06 0.0000500000
#   4   GO:0000280                            nuclear division       321          87    54.95 3.7e-06 0.0000500000
#   5   GO:0048285                           organelle fission       363          96    62.14 3.7e-06 0.0000500000
#   6   GO:1904872 regulation of telomerase RNA localizatio...        18          12     3.08 4.1e-06 0.0000500000
#   7   GO:0051301                               cell division       519         128    88.84 4.9e-06 0.0000500000
#   8   GO:0000070        mitotic sister chromatid segregation       158          50    27.05 5.0e-06 0.0000500000
#   9   GO:0000819                sister chromatid segregation       193          58    33.04 5.6e-06 0.0000500000
#   10  GO:0070887      cellular response to chemical stimulus      1786         371   305.73 6.0e-06 0.0000500000
#   11  GO:0000278                          mitotic cell cycle       732         170   125.31 7.2e-06 0.0000500000
#   12  GO:0071310      cellular response to organic substance      1348         288   230.75 9.0e-06 0.0000500000
#   13  GO:0014070         response to organic cyclic compound       609         145   104.25 9.0e-06 0.0000500000




## for the genes enriched in cell division BP, how many of the 2137 DETs overlapping between
## spleen and liver age category, are in our DETs and also annotated as responsible for  
## cell division
table(inGO = allLocusGO$go_id%in%"GO:0051301", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Age_2category_Young_vs_Adult"]], ### overlaping genes (2137)
                  DETs_ALL[["spleen:Age_2category_Young_vs_Adult"]]))

## 99 DETs


## listing the names of ten of the cell division DETs that are in 
##the 2137 overlapping liver and spleen DETs in age category

cat("cell division DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0051301" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Age_2category_Young_vs_Adult", 
                                                                         "spleen:Age_2category_Young_vs_Adult")])
          ]
        ), 
        decreasing = TRUE
      ), 10
    )), "\n")

## cell division DETs: ANAPC10 ANAPC7 ASPM AURKB BABAM1 BUB1 BUB1B CCNA2 CCNB1 CCNB2 
## protein IDs on NCBI: NP_001243635.1 NP_001131136.2 NP_001193775.1 NP_001243763.1  NP_001028721.1  




## doing the test for the intersection of liver and spleen sex category BP DETs

TOGO.all.onto("BP", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Sex_Male_vs_Female"]], 
                        DETs_ALL[["spleen:Sex_Male_vs_Female"]]),
              gene2GO) %>% gene.table.topGO()



#   GO.ID                                        Term Annotated Significant Expected p.value        fdr
#   1   GO:0006999                   nuclear pore organization        15           2     0.04 0.00085 0.02742105
#   2   GO:0044249               cellular biosynthetic process      5193          24    15.20 0.00116 0.02742105
#   3   GO:1901576      organic substance biosynthetic process      5252          24    15.37 0.00142 0.02742105
#   4   GO:0009058                        biosynthetic process      5283          24    15.46 0.00157 0.02742105
#   5   GO:1904385            cellular response to angiotensin        24           2     0.07 0.00220 0.02742105
#   6   GO:1990776                     response to angiotensin        27           2     0.08 0.00278 0.02742105
#   7   GO:0008057            eye pigment granule organization         1           1     0.00 0.00293 0.02742105
#   8   GO:0033386 geranylgeranyl diphosphate biosynthetic ...         1           1     0.00 0.00293 0.02742105
#   9   GO:1900736 regulation of phospholipase C-activating...         1           1     0.00 0.00293 0.02742105
#   10  GO:1900738 positive regulation of phospholipase C-a...         1           1     0.00 0.00293 0.02742105
#   11  GO:0007200 phospholipase C-activating G protein-cou...        37           2     0.11 0.00518 0.02742105
#   12  GO:0033384    geranyl diphosphate biosynthetic process         2           1     0.01 0.00585 0.02742105
#   13  GO:0033385 geranylgeranyl diphosphate metabolic pro...         2           1     0.01 0.00585 0.02742105
#   14  GO:0034159 regulation of toll-like receptor 8 signa...         2           1     0.01 0.00585 0.02742105
#   15  GO:0034161 positive regulation of toll-like recepto...         2           1     0.01 0.00585 0.02742105



## for the genes enriched in cellular biosynthetic process BP, how many of the 48 DETs overlapping between
## spleen and liver sex category, are in our DETs and also annotated as responsible for  
## cellular biosynthetic process
table(inGO = allLocusGO$go_id%in%"GO:0044249", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Sex_Male_vs_Female"]], ### overlaping genes (48)
                  DETs_ALL[["spleen:Sex_Male_vs_Female"]]))

## 1 DET


## listing the name of the cellular biosynthetic process DET that is in 
##the overlap of liver and spleen DETs in sex category

cat("cellular biosynthetic process DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0044249" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Sex_Male_vs_Female", 
                                                                         "spleen:Sex_Male_vs_Female")])
          ]
        ), 
        decreasing = TRUE
      ), 1
    )), "\n")

##  cellular biosynthetic process DETs: DEGS1
## protein IDs on NCBI: NP_001308470.1





#### looping over the DE gene sets we don't want to run the ####
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

names(CC_enrichment_tables)

lapply(CC_enrichment_tables, head)

CC_enrichment_tables["spleen:rpmh_scaled"]

CC_enrichment_tables["liver:rpmh_scaled"]


## doing one test for the intersection of liver and spleen infection intensity CC DETs

TOGO.all.onto("CC", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:rpmh_scaled"]], ## 133 overlapping genes
                        DETs_ALL[["spleen:rpmh_scaled"]]),
              gene2GO) %>% gene.table.topGO()



## ##       GO.ID                             Term Annotated Significant Expected p.value     fdr
## ## GO:0005614              interstitial matrix         6           2     0.06  0.0014 0.06666667
## ## GO:0005581                  collagen trimer        53           4     0.53  0.0019 0.06666667
## ## GO:0005883                    neurofilament         7           2     0.07  0.0020 0.06666667
## ## GO:0031012             extracellular matrix       318           9     3.16  0.0043 0.08800000
## ## GO:0030312 external encapsulating structure       319           9     3.17  0.0044 0.08800000
## ## GO:0005840                         ribosome       217           7     2.15  0.0058 0.09666667
## ## > 


## for the genes enriched in interstitial matrix CC, how many of the 133 DETs overlapping between
## spleen and liver infection intensity category, are in our DETs and also annotated as being in the cellular 
## compartment interstitial matrix
table(inGO = allLocusGO$go_id%in%"GO:0005614", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:rpmh_scaled"]], ### overlaping genes (133)
                  DETs_ALL[["spleen:rpmh_scaled"]]))

                ## 2 DETs
                

## listing the name of the 2 interstitial matrix DETs that is in 
##the 133 overlap of liver and spleen DETs in infection intensity category CC DETs

cat("interstitial matrix DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0005614" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:rpmh_scaled", 
                                                                         "spleen:rpmh_scaled")])
          ]
        ), 
        decreasing = TRUE
      ), 2
    )), "\n")

##  interstitial matrix DETs: CCDC80 COL14A1
## protein IDs on NCBI: NP_955805.1  NP_001371876.1



## doing one test for the intersection of liver and spleen season category CC DETs

TOGO.all.onto("CC", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Season_Rainy_vs_Dry"]], ##  overlapping genes (251)
                        DETs_ALL[["spleen:Season_Rainy_vs_Dry"]]),
              gene2GO) %>% gene.table.topGO()



#         GO.ID                                    Term         Annotated Significant Expected p.value    fdr
#   1  GO:0072562                         blood microparticle        65           8     1.20 2.3e-05 0.000975000
#   2  GO:0005576                        extracellular region      2308          67    42.47 3.5e-05 0.000975000
#   3  GO:0005777                                  peroxisome       112          10     2.06 3.9e-05 0.000975000
#   4  GO:0042579                                   microbody       112          10     2.06 3.9e-05 0.000975000
#   5  GO:0005782                          peroxisomal matrix        38           6     0.70 6.1e-05 0.001016667
#   6  GO:0031907                             microbody lumen        38           6     0.70 6.1e-05 0.001016667
#   7  GO:0005615                         extracellular space      1952          58    35.92 7.8e-05 0.001114286
#   8  GO:0034366 spherical high-density lipoprotein parti...         6           3     0.11 0.00012 0.001444444
#   9  GO:0034364           high-density lipoprotein particle        15           4     0.28 0.00013 0.001444444
#   10 GO:0034361       very-low-density lipoprotein particle        17           4     0.31 0.00022 0.002000000



## for the genes enriched in extracellular region CC, how many of the 251 DETs overlapping between
## spleen and liver season category, are in our DETs and also annotated as being in the cellular 
## compartment extracellular region
table(inGO = allLocusGO$go_id%in%"GO:0005576", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Season_Rainy_vs_Dry"]], ### overlaping genes ()
                  DETs_ALL[["spleen:Season_Rainy_vs_Dry"]]))

## 40 DETs


## listing the name of 10 of the 40 extracellular region DETs that is in 
##the 251 overlap of liver and spleen DETs in season category CC DETs

cat("extracellular region DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0005576" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Season_Rainy_vs_Dry", 
                                                                         "spleen:Season_Rainy_vs_Dry")])
          ]
        ), 
        decreasing = TRUE
      ), 10
    )), "\n")

##  extracellular region DETs: ADAMTS1 APOA1 APOC1 APOC2 APOC3 APOLD1 C6 CAT COL4A1 CPN1 
## protein IDs on NCBI: NP_008919.3  NP_000030.1 NP_001307994.1 NP_000474.2  NP_000031.1  NP_001123887.1  



## doing one test for the intersection of liver and spleen age CC DETs

TOGO.all.onto("CC", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Age_2category_Young_vs_Adult"]], ##  overlapping genes (2137)
                        DETs_ALL[["spleen:Age_2category_Young_vs_Adult"]]),
              gene2GO) %>% gene.table.topGO()



#         GO.ID                                     Term         Annotated Significant Expected p.value    fdr
#   1   GO:0000793                        condensed chromosome       219          64    37.28 4.4e-06 0.0004400000
#   2   GO:0032991                  protein-containing complex      4419         832   752.33 2.2e-05 0.0007666667
#   3   GO:0005694                                  chromosome      1306         276   222.35 2.3e-05 0.0007666667
#   4   GO:0043228              non-membrane-bounded organelle      3707         703   631.12 7.0e-05 0.0014000000
#   5   GO:0043232 intracellular non-membrane-bounded organ...      3707         703   631.12 7.0e-05 0.0014000000
#   6   GO:0072562                         blood microparticle        65          24    11.07 9.4e-05 0.0015666667
#   7   GO:0000775              chromosome, centromeric region       208          57    35.41 0.00011 0.0015714286
#   8   GO:0031974                     membrane-enclosed lumen      4440         825   755.91 0.00020 0.0019090909
#   9   GO:0043233                             organelle lumen      4440         825   755.91 0.00020 0.0019090909
#   10  GO:0070013               intracellular organelle lumen      4440         825   755.91 0.00020 0.0019090909
#   11  GO:0031982                                     vesicle      2730         526   464.78 0.00021 0.0019090909
#   12  GO:0005615                         extracellular space      1952         386   332.33 0.00025 0.0019285714
#   13  GO:0098687                          chromosomal region       325          80    55.33 0.00026 0.0019285714
#   14  GO:0072686                             mitotic spindle       155          44    26.39 0.00027 0.0019285714
#   15  GO:0000779    condensed chromosome, centromeric region       151          43    25.71 0.00029 0.0019333333



## for the genes enriched in extracellular space CC, how many of the 2137 DETs overlapping between
## spleen and liver age category, are in our DETs and also annotated as being in the cellular 
## compartment extracellular space
table(inGO = allLocusGO$go_id%in%"GO:0005615", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Age_2category_Young_vs_Adult"]], ### overlaping genes (2137)
                  DETs_ALL[["spleen:Age_2category_Young_vs_Adult"]]))

## 151 DETs


## listing the name of 10 of the 151 extracellular space DETs that is in 
##the 2137 overlap of liver and spleen DETs in age category CC DETs

cat("extracellular space DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0005615" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Age_2category_Young_vs_Adult", 
                                                                         "spleen:Age_2category_Young_vs_Adult")])
          ]
        ), 
        decreasing = TRUE
      ), 10
    )), "\n")

##  extracellular space DETs: A1BG ACTB ADAM9 AFM AGT AHSG AKR1A1 ALDOA ANGPT1 ANGPTL4
## protein IDs on NCBI:  NP_570602.2  NP_001092.1  NP_003807.1  NP_001124.1 NP_001369746.2 





## doing one test for the intersection of liver and spleen sex CC DETs

TOGO.all.onto("CC", DETs_ALL[["overall"]],
              intersect(DETs_ALL[["liver:Sex_Male_vs_Female"]], ##  overlapping genes (48)
                        DETs_ALL[["spleen:Sex_Male_vs_Female"]]),
              gene2GO) %>% gene.table.topGO()


#         GO.ID                    Term         Annotated Significant Expected p.value  fdr
#   1 GO:0031080     nuclear pore outer ring        10           2     0.03 0.00039 0.03000
#   2 GO:0030864 cortical actin cytoskeleton        62           3     0.19 0.00083 0.03000
#   3 GO:0031143                pseudopodium        15           2     0.05 0.00090 0.03000
#   4 GO:0030863       cortical cytoskeleton        87           3     0.26 0.00221 0.05525


## for the genes enriched in nuclear pore outer ring CC, how many of the 48 DETs overlapping between
## spleen and liver sex category, are in our DETs and also annotated as being in the cellular 
## compartment nuclear pore outer ring
table(inGO = allLocusGO$go_id%in%"GO:0031080", 
      tested = allLocusGO$entrezgene_accession%in%DETs_ALL[["overall"]],
      DE = allLocusGO$entrezgene_accession%in%
        intersect(DETs_ALL[["liver:Sex_Male_vs_Female"]], ### overlaping genes (48)
                  DETs_ALL[["spleen:Sex_Male_vs_Female"]]))

## 2 DETs


## listing the name of the 2 nuclear pore outer ring DETs that is in 
##the 48 overlap of liver and spleen DETs in sex category CC DETs

cat("nuclear pore outer ring DETs:", 
    names(head(
      sort(
        table(
          allLocusGO$entrezgene_accession[allLocusGO$go_id %in% "GO:0005615" &
                                            allLocusGO$entrezgene_accession %in% DETs_ALL$overall &
                                            allLocusGO$entrezgene_accession %in% 
                                            Reduce(intersect, DETs_ALL[c("liver:Sex_Male_vs_Female", 
                                                                         "spleen:Sex_Male_vs_Female")])
          ]
        ), 
        decreasing = TRUE
      ), 2
    )), "\n")

##  nuclear pore outer ring DETs: F2 IGFBP2 
## protein IDs on NCBI:  NP_000497.1  NP_000588.3 
