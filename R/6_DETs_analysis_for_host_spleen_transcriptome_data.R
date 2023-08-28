# DETs analysis for host spleen transcriptome from 3' Tag RNAseqs
# differential expression tests are based on a negative binomial 
#generalized linear model

# separate analysis for liver and spleen 
library(DESeq2)
library(ggplot2)
library(dplyr)
library(apeglm)


# Filtering metadata to keep only rows for spleen samples
metadata_spleen <- metadata %>%
                   filter(Organ == 'Spleen')

# separating host liver from host spleen counts
# spleen host counts
host_counts_spleen <- host_counts[, metadata_spleen$ID]


####  DETs analysis for spleen ####

# constructing the DESeqdataset object with Sex as condition of test
dds_spleen <- DESeqDataSetFromMatrix(countData = host_counts_spleen,
                                    colData = metadata_spleen,
                                    design = ~Season+Age_2category+Sex,
                                    tidy = FALSE)

# viewing the object
dds_spleen

# differential expression analysis (uses wald test statistics)
dds_spleen <- DESeq(dds_spleen)
        # does everything from normalization to linear modeling

        #estimateSizeFactors
        #This calculates the relative library depth of each sample 
        
        #estimateDispersions
        #estimates the dispersion of counts for each gene 
        
        #nbinomWaldTest
        #calculates the significance of coefficients in a 
        #Negative Binomial GLM using the size and dispersion outputs


# getting the results table
res_spleen <- results(dds_spleen)

# viewing the table
res_spleen

# head(results(dds_spleen, tidy = TRUE))

# summary of DTE
summary(res_spleen)

# sorting table by p values
res_spleen <- res_spleen[order(res_spleen$padj),]
head(res_spleen)

write.csv(res_spleen, file = "intermediateData/spleen_DETs.csv")

# adjusted p-values less than 0.1
sum(res_spleen$padj < 0.1, na.rm=TRUE)  # =68 DETs


# plotting the results 

#MA plots
pdf("plots/spleen_DTEs_MA_plot.pdf", width = 12, height = 12)
plotMA(res_spleen, ylim = c(-2, 2))
dev.off()

        #shows the log2 fold changes attributable to a given variable 
        #over the mean of normalized counts for all the samples in the
        #DESeqDataSet. Points are colored blue if the adjusted p value 
        #is less than 0.1 (statistically significant points). 
        #Points which fall out of the window are plotted as open 
        #triangles pointing either up or down.

# getting the resLFC - Log fold change shrinkage for visualization and ranking

# name of coefficient to shrink
resultsNames(dds_spleen)

res_spleen_LFC <- lfcShrink(dds_spleen, coef = "Sex_Male_vs_Female",
                           type = "apeglm")

res_spleen_LFC

      # It is more useful visualize the MA-plot for the shrunken log2 fold 
      # changes, which remove the noise associated with log2 fold changes 
      # from low count genes without requiring arbitrary filtering thresholds.

# plotting the LFC
pdf("plots/spleen_LFC_DTEs_MA_plot.pdf", width = 12, height = 12)
plotMA(res_spleen_LFC, ylim = c(-2, 2))
dev.off()


# Plot counts
# plotCounts function compares the normalized counts between the comparison groups,
# sex Male vs female, for the genes

# plotting transcript with lowest p adjusted value from results table
pdf("plots/spleen_DTEs_plotcounts_for_min_padjvalue_transcript.pdf", width = 12, height = 12)
plotCounts(dds_spleen, gene=which.min(res_spleen$padj), 
           intgroup="Sex")

dev.off()

#alternatively plotting using ggplot2
pdf("plots/spleen_DTEs_plotcounts_for_min_padjvalue2_transcript.pdf", width = 12, height = 12)

ggplot(plotCounts(dds_spleen, gene=which.min(res_spleen$padj), 
                  intgroup="Sex", returnData = TRUE), 
       aes(x=Sex, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

dev.off()

# alternatively plotting the counts for individual transcripts
pdf("plots/spleen_DTEs_plotcounts_for_LOC107508184_transcript.pdf", width = 12, height = 12)

plotCounts(dds_spleen, gene = "LOC107508184", intgroup = "Sex" )

dev.off()

# getting information on variables and tests used 
mcols(res_spleen)$description


# volcano plots
pdf("plots/spleen_DTEs_volcano_plot.pdf", width = 12, height = 12)

with(res_spleen, plot(log2FoldChange, -log10(pvalue), 
                     pch=20, main="Volcano plot", xlim=c(-3,3)))

# Adding colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_spleen, padj<.01 ), 
     points(log2FoldChange, -log10(pvalue), 
            pch=20, col="blue"))

with(subset(res_spleen, (padj<.01 & abs(log2FoldChange)>2)),
     points(log2FoldChange, -log10(pvalue), 
            pch=20, col="red"))

dev.off()

# PCA plots
#First we need to transform the raw count data
#vst function performs variance stabilizing transformation

vs_spleen_counts_data <- vst(dds_spleen, blind=FALSE)

pdf("plots/spleen_DTEs_PCA_plot.pdf", width = 12, height = 12)

plotPCA(vs_spleen_counts_data, intgroup = "Sex")

dev.off()

        # shows how samples group by sex