
#Script to create venn diagram to check overlap between upregulated and downregulated genes between the 3 celllines for 3 days, BETi vs DMSO.
library(dplyr)
library(tidyverse)
library(ggVennDiagram)
library(ggplot2)

A673_DE_20days <- read.csv("../../lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Results_A673_20days_BETI_vs_DMSO/A673_condition_treated_results.csv", row.names=1)
CHLA10_DE_20days <- read.csv("../../lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Results_CHLA10_20days_BETI_vs_DMSO/CHLA10_condition_treated_results.csv",row.names=1)
TC32_DE_20days <- read.csv("../../lawlorlab_shireen_2023.07.13_bulk_rnaseq_tc32_count_nf/additional_results/Results_TC32_20days_BETI_vs_DMSO/TC32_condition_treated_results.csv",row.names=1)

A673_DE_20days_upreg <- A673_DE_20days %>% filter(log2FoldChange > 0, padj < 0.01)
CHLA10_DE_20days_upreg <- CHLA10_DE_20days %>% filter(log2FoldChange > 0, padj < 0.01)
TC32_DE_20days_upreg <- TC32_DE_20days %>% filter(log2FoldChange > 0, padj < 0.01)

A673_DE_20days_downreg <- A673_DE_20days %>% filter(log2FoldChange < 0, padj < 0.01)
CHLA10_DE_20days_downreg <- CHLA10_DE_20days %>% filter(log2FoldChange < 0, padj < 0.01)
TC32_DE_20days_downreg <- TC32_DE_20days %>% filter(log2FoldChange < 0, padj < 0.01)

write.table(A673_DE_20days_upreg, "A673_20days_BETi_DMSO_upreg.txt", sep="\t", quote=FALSE) 
write.table(CHLA10_DE_20days_upreg, "CHLA10_20days_BETi_DMSO_upreg.txt", sep="\t", quote=FALSE) 
write.table(TC32_DE_20days_upreg, "TC32_20days_BETi_DMSO_upreg.txt", sep="\t", quote=FALSE)

write.table(A673_DE_20days_downreg, "A673_20days_BETi_DMSO_downreg.txt", sep="\t", quote=FALSE) 
write.table(CHLA10_DE_20days_downreg, "CHLA10_20days_BETi_DMSO_downreg.txt", sep="\t", quote=FALSE) 
write.table(TC32_DE_20days_downreg, "TC32_20days_BETi_DMSO_downreg.txt", sep="\t", quote=FALSE)

#Make the plot
upreg <- list(A673 = rownames(A673_DE_20days_upreg), CHLA10 = rownames(CHLA10_DE_20days_upreg), TC32 = rownames(TC32_DE_20days_upreg))
downreg <- list(A673 = rownames(A673_DE_20days_downreg), CHLA10 = rownames(CHLA10_DE_20days_downreg), TC32 = rownames(TC32_DE_20days_downreg))

pdf("Venndiag_upreg_20days_BETi_DMSO.pdf")
ggVennDiagram(upreg, label = "count") + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()

pdf("Venndiag_downreg_20days_BETi_DMSO.pdf")
ggVennDiagram(downreg, label = "count") + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()

Upreg_intersect_genes <- Reduce(intersect, list(A673 = rownames(A673_DE_20days_upreg), CHLA10 = rownames(CHLA10_DE_20days_upreg), TC32 = rownames(TC32_DE_20days_upreg)))
Downreg_intersect_genes <- Reduce(intersect, list(A673 = rownames(A673_DE_20days_upreg), CHLA10 = rownames(CHLA10_DE_20days_upreg), TC32 = rownames(TC32_DE_20days_upreg)))

write.table(Upreg_intersect_genes, "Upreg_intersect_genes.txt", sep="\n", quote=FALSE)
write.table(Downreg_intersect_genes, "Downreg_intersect_genes.txt", sep="\n", quote=FALSE)

##################----------------------######################----------------#############



