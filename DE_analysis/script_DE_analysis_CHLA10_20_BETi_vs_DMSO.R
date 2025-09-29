
library("pasilla")
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")

library("DEGreport")
library("pasilla")
library("ggrepel")

sampleData <- read.table("/active/taylor_s/people/nkatiy/CP-E_sarcoma_BETi/lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/resources/paired_end_samplesheet_CHLA10_20days.csv", sep=",", header=T)
rawCounts <- read.delim("/active/taylor_s/people/nkatiy/CP-E_sarcoma_BETi/lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Read_count_matrix_cell_line_CHLA10_20days.txt", sep="\t", header=T, row.names=1)

# Convert count data to a matrix of appropriate form that DEseq2 can read
rawCounts <- as.matrix(rawCounts)
cts <- rawCounts

rownames(sampleData) <- sampleData$id

coldata <- sampleData[,c("cell_line", "condition", "id")]
rownames(coldata) <- sampleData$id

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

dds

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

#Prefiltering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Note on factor levels
dds$condition <- factor(dds$condition, levels = c("DMSO", "BETi"))

#Genes differentially expressed between conditions
## probably the most important command :)
dds <- DESeq(dds)

res <- results(dds)
res

#res <- results(dds, name="condition_1uM_vs_0uM")
#res <- results(dds, contrast=c("drug_amt","1uM","0uM"))

#Log fold change shrinkage for visualization and ranking
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_BETi_vs_DMSO", type="apeglm")
# noticed any differences to the values in lfcSE compared to res ?
resLFC

#p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]

summary(res)

#res05 <- results(dds, alpha=0.05)
#summary(res05)

#sum(res05$padj < 0.05, na.rm=TRUE)

#MA plot
pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Plot_MA.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Plot_MA_LFC.pdf")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Plot_Counts.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="cell_line", 
                returnData=TRUE)

library("ggplot2")
pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Plot_ggplot_count.pdf")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
dev.off()

mcols(res)$description

#Exporting results to CSV files
write.csv(as.data.frame(resOrdered), 
          file="./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/CHLA10_condition_treated_results.csv")

#resSig <- subset(resOrdered, padj < 0.1)
#resSig

#Data transformations and visualization

vsd <- vst(dds)
head(assay(vsd), 3)

#Data quality assessment by sample clustering and visualization
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("cell_line", "condition")])

pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Heatmap_sample_clustering.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Heatmap_sample_clustering_with_labels.pdf")
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$cell_line, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Heatmap_sample_distance.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#Principal component plot of the samples
pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/PCA_plot.pdf")
plotPCA(vsd, intgroup=c("cell_line", "condition"))
dev.off()

#Contrasts
degs <- degComps(dds, combs = "condition",
                 contrast = list( c("condition", "DMSO", "BETi")))
names(degs)

significants(degs, fc = 0, fdr = 0.01, full = TRUE)

#Gene Plots
# example
pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/DEG_plot.pdf")
degPlot(dds = dds, res = res, n = 6, xs = "condition")
dev.off()

df_res <- rownames_to_column(as.data.frame(res), var ="Gene")
write.table(df_res, "./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/DEG_results.txt", sep="\t", quote=FALSE)

## Obtain logical vector where TRUE values denote padj values < 0.01

res_tableOE_tb <- df_res %>% 
                  mutate(threshold_OE = padj < 0.01)

## Create a column to indicate which genes to label
res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")

res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$Gene[1:10]

#View(res_tableOE_tb)
#Next, we plot it as before with an additiona layer for geom_text_repel() wherein we can specify the column of gene labels we just created.

pdf("./additional_results/Results_CHLA10_20days_BETI_vs_DMSO/Volcano_plot_CHLA10.pdf")
ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(colour = threshold_OE)) +
        geom_text_repel(aes(label = genelabels)) +
        ggtitle("CHLA10 overexpression") +
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25))) 
dev.off()


