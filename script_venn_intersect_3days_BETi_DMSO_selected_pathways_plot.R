
#Script to create venn diagram to check overlap between upregulated and downregulated genes between the 3 celllines for 3 days, BETi vs DMSO.
library(dplyr)
library(tidyverse)
library(ggVennDiagram)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)

A673_DE_3days <- read.csv("../tables/A673_condition_treated_results.csv", row.names=1)
CHLA10_DE_3days <- read.csv("../tables/CHLA10_condition_treated_results.csv",row.names=1)
TC32_DE_3days <- read.csv("../tables/TC32_condition_treated_results.csv",row.names=1)

A673_DE_3days_upreg <- A673_DE_3days %>% filter(log2FoldChange > 0, padj < 0.05)
CHLA10_DE_3days_upreg <- CHLA10_DE_3days %>% filter(log2FoldChange > 0, padj < 0.05)
TC32_DE_3days_upreg <- TC32_DE_3days %>% filter(log2FoldChange > 0, padj < 0.05)

A673_DE_3days_downreg <- A673_DE_3days %>% filter(log2FoldChange < 0, padj < 0.05)
CHLA10_DE_3days_downreg <- CHLA10_DE_3days %>% filter(log2FoldChange < 0, padj < 0.05)
TC32_DE_3days_downreg <- TC32_DE_3days %>% filter(log2FoldChange < 0, padj < 0.05)

#Make the plot
upreg <- list(A673 = rownames(A673_DE_3days_upreg), CHLA10 = rownames(CHLA10_DE_3days_upreg), TC32 = rownames(TC32_DE_3days_upreg))
downreg <- list(A673 = rownames(A673_DE_3days_downreg), CHLA10 = rownames(CHLA10_DE_3days_downreg), TC32 = rownames(TC32_DE_3days_downreg))

A673_upreg = rownames(A673_DE_3days_upreg)
CHLA10_upreg = rownames(CHLA10_DE_3days_upreg)
TC32_upreg = rownames(TC32_DE_3days_upreg)

A673_downreg = rownames(A673_DE_3days_downreg)
CHLA10_downreg = rownames(CHLA10_DE_3days_downreg)
TC32_downreg = rownames(TC32_DE_3days_downreg)

##################-------------------------#####################

upreg_intersect_genes <- intersect(intersect(A673_upreg, CHLA10_upreg),TC32_upreg)
downreg_intersect_genes <- intersect(intersect(A673_downreg, CHLA10_downreg),TC32_downreg)

##################--------------------#################3

dotplot_custom <- function(data, showCategory){

data_df <- as.data.frame(data)
n = showCategory
if(nrow(data_df) >=20)
{
	n=20
}
else
{
	n=nrow(data_df)
}

data_df <- data_df[1:n,]
data_df$log10_padj <- -log10(data_df$p.adjust)
generatio_transformed <- sapply(data_df$GeneRatio, function(x) eval(parse(text=x)))
data_df$GeneRatio <- generatio_transformed

ggplot(data = data_df, aes(x = log10_padj, y = reorder(Description,log10_padj),
                        color = GeneRatio, size = Count)) +
  geom_point() +
  scale_color_gradient(name = "GeneRatio", low = "red", high = "blue") +
  theme_bw() +
  ylab("") +
  xlab("-log10(p.adjust)") +
  ggtitle("Enrichment analysis") + expand_limits(x = 0, y = 0)

}

m_REACTOME_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") %>% dplyr::select(gs_name, ensembl_gene)

em_REACTOME_upreg <- enricher(upreg_intersect_genes, TERM2GENE=m_REACTOME_t2g)
em_REACTOME_downreg <- enricher(downreg_intersect_genes, TERM2GENE=m_REACTOME_t2g)

vals_reactome_up = c("REACTOME_RHO_GTPASE_CYCLE", "REACTOME_SIGNALING_BY_VEGF", "REACTOME_PI_METABOLISM", "REACTOME_CLATHRIN_MEDIATED_ENDOCYTOSIS", "REACTOME_RAC1_GTPASE_CYCLE")

y_reactome_up <- em_REACTOME_upreg[em_REACTOME_upreg$ID %in% vals_reactome_up, asis = T]

pdf("../Plots/DotPlot_MSigdb_REACTOME_3days_upreg_v3.pdf")
dotplot_custom(y_reactome_up, showCategory=10)
dev.off()

vals_reactome_down = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION", "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES", "REACTOME_COLLAGEN_DEGRADATION", "REACTOME_ELASTIC_FIBRE_FORMATION")

y_reactome_down <- em_REACTOME_downreg[em_REACTOME_downreg$ID %in% vals_reactome_down, asis = T]

pdf("../Plots/DotPlot_MSigdb_REACTOME_3days_downreg_v3.pdf", width=10)
dotplot_custom(y_reactome_down, showCategory=10)
dev.off()

####################-------------------#####################
m_GO_BP_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, ensembl_gene)

em_GO_BP_upreg <- enricher(upreg_intersect_genes, TERM2GENE=m_GO_BP_t2g)
em_GO_BP_downreg <- enricher(downreg_intersect_genes, TERM2GENE=m_GO_BP_t2g)

go_bp_vals_up <- c("GOBP_ENDOSOMAL_TRANSPORT", "GOBP_VESICLE_LOCALIZATION", "GOBP_PHOSPHOLIPID_BIOSYNTHETIC_PROCESS", "GOBP_REGULATION_OF_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION_OR_REMOVAL", "GOBP_PROTEIN_LOCALIZATION_TO_CELL_PERIPHERY")

y_up <- em_GO_BP_upreg[em_GO_BP_upreg$ID %in% go_bp_vals_up, asis = T]

pdf("../Plots/DotPlot_MSigdb_GO_BP_3days_upreg_v3.pdf", width=10)
dotplot_custom(y_up, showCategory=20)
dev.off()

go_bp_vals_down <- c("GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY", "GOBP_REGULATION_OF_SMALL_GTPASE_MEDIATED_SIGNAL_TRANSDUCTION", "GOBP_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS", "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION")

y_down <- em_GO_BP_downreg[em_GO_BP_downreg$ID %in% go_bp_vals_down, asis = T]

pdf("../Plots/DotPlot_MSigdb_GO_BP_3days_downreg_v3.pdf", width=10)
dotplot_custom(y_down, showCategory=20)
dev.off()

####################--------#############################
m_H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)

em_H_upreg <- enricher(upreg_intersect_genes, TERM2GENE=m_H_t2g)
em_H_downreg <- enricher(downreg_intersect_genes, TERM2GENE=m_H_t2g)

pdf("../Plots/DotPlot_MSigdb_hallmark_3days_upreg_v3.pdf", width=7)
dotplot_custom(em_H_upreg, showCategory=20)
dev.off()

pdf("../Plots/DotPlot_MSigdb_hallmark_3days_downreg_v3.pdf", width=7)
dotplot_custom(em_H_downreg, showCategory=20)
dev.off()

