
#Script to create venn diagram to check overlap between upregulated and downregulated genes between the 3 celllines for 3 days, BETi vs DMSO.
#Author : Neerja Katiyar

library(dplyr)
library(tidyverse)
library(ggVennDiagram)
library(ggplot2)

A673_DE_3days <- read.csv("../../../lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Results_A673_3days_BETI_vs_DMSO/A673_condition_treated_results.csv", row.names=1)
CHLA10_DE_3days <- read.csv("../../../lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Results_CHLA10_3days_BETI_vs_DMSO/CHLA10_condition_treated_results.csv",row.names=1)
TC32_DE_3days <- read.csv("../../../lawlorlab_shireen_2023.07.13_bulk_rnaseq_tc32_count_nf/additional_results/Results_TC32_3days_BETI_vs_DMSO/TC32_condition_treated_results.csv",row.names=1)

A673_DE_3days_upreg <- A673_DE_3days %>% filter(log2FoldChange > 0, padj < 0.05)
CHLA10_DE_3days_upreg <- CHLA10_DE_3days %>% filter(log2FoldChange > 0, padj < 0.05)
TC32_DE_3days_upreg <- TC32_DE_3days %>% filter(log2FoldChange > 0, padj < 0.05)

A673_DE_3days_downreg <- A673_DE_3days %>% filter(log2FoldChange < 0, padj < 0.05)
CHLA10_DE_3days_downreg <- CHLA10_DE_3days %>% filter(log2FoldChange < 0, padj < 0.05)
TC32_DE_3days_downreg <- TC32_DE_3days %>% filter(log2FoldChange < 0, padj < 0.05)

write.table(A673_DE_3days_upreg, "A673_3days_BETi_DMSO_upreg.txt", sep="\t", quote=FALSE) 
write.table(CHLA10_DE_3days_upreg, "CHLA10_3days_BETi_DMSO_upreg.txt", sep="\t", quote=FALSE) 
write.table(TC32_DE_3days_upreg, "TC32_3days_BETi_DMSO_upreg.txt", sep="\t", quote=FALSE)

write.table(A673_DE_3days_downreg, "A673_3days_BETi_DMSO_downreg.txt", sep="\t", quote=FALSE) 
write.table(CHLA10_DE_3days_downreg, "CHLA10_3days_BETi_DMSO_downreg.txt", sep="\t", quote=FALSE) 
write.table(TC32_DE_3days_downreg, "TC32_3days_BETi_DMSO_downreg.txt", sep="\t", quote=FALSE)

#Make the plot
upreg <- list(A673 = rownames(A673_DE_3days_upreg), CHLA10 = rownames(CHLA10_DE_3days_upreg), TC32 = rownames(TC32_DE_3days_upreg))
downreg <- list(A673 = rownames(A673_DE_3days_downreg), CHLA10 = rownames(CHLA10_DE_3days_downreg), TC32 = rownames(TC32_DE_3days_downreg))

pdf("Venndiag_upreg_3days_BETi_DMSO.pdf")
ggVennDiagram(upreg, label = "count") + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()

pdf("Venndiag_downreg_3days_BETi_DMSO.pdf")
ggVennDiagram(downreg, label = "count") + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()

##############------------##################
library(eulerr)

euler_upreg_plot <- euler(upreg)

#Venn diagrams for paper.

pdf("Venn_prop_upreg_3days_BETi_DMSO.pdf")
plot(euler_upreg_plot, quantities = TRUE, labels = list(font=4))
dev.off()

euler_downreg_plot <- euler(downreg)
pdf("Venn_prop_downreg_3days_BETi_DMSO.pdf")
plot(euler_downreg_plot, quantities = TRUE, labels = list(font=4))
dev.off()


