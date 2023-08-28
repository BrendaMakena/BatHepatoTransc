# DETs analysis for host liver transcriptome from 3' Tag RNAseqs
# differential expression tests are based on a negative binomial 
#generalized linear model

# separate analysis for liver and spleen 
BiocManager::install("apeglm")
library(DESeq2)
library(ggplot2)
library(dplyr)
library(apeglm)

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

# constructing the DESeqdataset object with Sex as condition of test
dds_liver <- DESeqDataSetFromMatrix(countData = host_counts_liver,
                              colData = metadata_liver,
                              design = ~Season+Age_2category+Sex,
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
res_liver <- results(dds_liver)

# viewing the table
res_liver

# head(results(dds_liver, tidy = TRUE))

# summary of DTE
summary(res_liver)

# sorting table by p values
res_liver <- res_liver[order(res_liver$padj),]
head(res_liver)

write.csv(res_liver, file = "intermediateData/liver_DETs.csv")

# adjusted p-values less than 0.1
sum(res_liver$padj < 0.1, na.rm=TRUE)  # =306 DETs


# plotting the results 

#MA plots
pdf("plots/liver_DTEs_MA_plot.pdf", width = 12, height = 12)
plotMA(res_liver, ylim = c(-2, 2))
dev.off()
      
      #shows the log2 fold changes attributable to a given variable 
      #over the mean of normalized counts for all the samples in the
      #DESeqDataSet. Points are colored blue if the adjusted p value 
      #is less than 0.1 (statistically significant points). 
      #Points which fall out of the window are plotted as open 
      #triangles pointing either up or down.

# getting the resLFC - Log fold change shrinkage for visualization and ranking

# name of coefficient to shrink
resultsNames(dds_liver)

res_liver_LFC <- lfcShrink(dds_liver, coef = "Sex_Male_vs_Female",
                           type = "apeglm")

res_liver_LFC

                # It is more useful visualize the MA-plot for the shrunken log2 fold 
                # changes, which remove the noise associated with log2 fold changes 
                # from low count genes without requiring arbitrary filtering thresholds.

# plotting the LFC
pdf("plots/liver_LFC_DTEs_MA_plot.pdf", width = 12, height = 12)
plotMA(res_liver_LFC, ylim = c(-2, 2))
dev.off()


# Plot counts
        # plotCounts function compares the normalized counts between the comparison groups,
        # sex Male vs female, for the genes

# plotting transcript with lowest padjusted value from results table
pdf("plots/liver_DTEs_plotcounts_for_min_padjvalue_transcript.pdf", width = 12, height = 12)
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

# alternatively plotting the counts for individual transcripts
pdf("plots/liver_DTEs_plotcounts_for_LOC107508184_transcript.pdf", width = 12, height = 12)

plotCounts(dds_liver, gene = "LOC107508184", intgroup = "Sex" )

dev.off()

# getting information on variables and tests used 
mcols(res_liver)$description


# volcano plots
pdf("plots/liver_DTEs_volcano_plot.pdf", width = 12, height = 12)

with(res_liver, plot(log2FoldChange, -log10(pvalue), 
                     pch=20, main="Volcano plot", xlim=c(-3,3)))

# Adding colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_liver, padj<.01 ), 
     points(log2FoldChange, -log10(pvalue), 
            pch=20, col="blue"))

with(subset(res_liver, (padj<.01 & abs(log2FoldChange)>2)),
     points(log2FoldChange, -log10(pvalue), 
            pch=20, col="red"))

dev.off()

# PCA plots
#First we need to transform the raw count data
#vst function performs variance stabilizing transformation

vs_liver_counts_data <- vst(dds_liver, blind=FALSE)

pdf("plots/liver_DTEs_PCA_plot.pdf", width = 12, height = 12)

plotPCA(vs_liver_counts_data, intgroup = "Sex")

dev.off()

      # shows how samples group by sex