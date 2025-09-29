
# Attach the library
library(clusterProfiler)

# Package that contains MSigDB gene sets in tidy format
library(msigdb)

# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)
library(dplyr)
library(ggplot2)

A673_3days <- read.csv("../tables/A673_condition_treated_results.csv", row.names=1)
CHLA10_3days <- read.csv("../tables/CHLA10_condition_treated_results.csv",row.names=1)

A673_3days$Gene <- rownames(A673_3days)
CHLA10_3days$Gene <- rownames(CHLA10_3days)

# create a mapped data frame we can join to the differential.
A673_3days_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = A673_3days$Gene,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(A673_3days, by = c("Ensembl" = "Gene"))

# create a mapped data frame we can join to the differential.
CHLA10_3days_mapped_df <- data.frame(
  gene_symbol = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = CHLA10_3days$Gene,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_symbol)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(CHLA10_3days, by = c("Ensembl" = "Gene"))

#duplicated gene symbols
A673_dup_gene_symbols <- A673_3days_mapped_df %>%
  dplyr::filter(duplicated(gene_symbol)) %>%
  dplyr::pull(gene_symbol)

#duplicated gene symbols
CHLA10_dup_gene_symbols <- CHLA10_3days_mapped_df %>%
  dplyr::filter(duplicated(gene_symbol)) %>%
  dplyr::pull(gene_symbol)

#Now let’s take a look at the rows associated with the duplicated gene symbols.
A673_3days_mapped_df %>%
  dplyr::filter(gene_symbol %in% A673_dup_gene_symbols) %>%
  dplyr::arrange(gene_symbol)

#Now let’s take a look at the rows associated with the duplicated gene symbols.
CHLA10_3days_mapped_df %>%
  dplyr::filter(gene_symbol %in% CHLA10_dup_gene_symbols) %>%
  dplyr::arrange(gene_symbol)

# Let's create a named vector ranked based on the log2 fold change values
A673_lfc_vector <- A673_3days_mapped_df$log2FoldChange
names(A673_lfc_vector) <- A673_3days_mapped_df$gene_symbol

# Let's create a named vector ranked based on the log2 fold change values
CHLA10_lfc_vector <- CHLA10_3days_mapped_df$log2FoldChange
names(CHLA10_lfc_vector) <- CHLA10_3days_mapped_df$gene_symbol

# We need to sort the log2 fold change values in descending order here
A673_genelist <- sort(A673_lfc_vector, decreasing = TRUE)
CHLA10_genelist <- sort(CHLA10_lfc_vector, decreasing = TRUE)

#
ggo_BP_A673_3days <- gseGO(gene = A673_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)
ggo_BP_CHLA10_3days <- gseGO(gene = CHLA10_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)

#write.table(ggo_BP_A673_3days, "./Results_20_vs_3days_BETi_DMSO/GO_BP_A673_3days_GSEA.txt", sep="\t", quote=FALSE)

dotplot_custom <- function(data, showCategory){

data_df <- as.data.frame(data)

data_df$num_core_genes <- strsplit(as.character(data_df$core_enrichment),'/',fixed=TRUE)
data_df$GeneRatio <- (lengths(gregexpr("\\W+", data_df$num_core_genes)) - 1)/ data_df$setSize

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

ggplot(data = data_df, aes(x = log10_padj, y = reorder(Description,log10_padj),
			                                                      color = GeneRatio, size = NES)) +
				       geom_point() +
				           scale_color_gradient(name = "GeneRatio", low = "red", high = "blue") +
					         theme_bw() +
						         ylab("") +
							           xlab("-log10(p.adjust)") +
								               ggtitle("Enrichment analysis") + expand_limits(x=0, y=0)
}


library(enrichplot)
vals_A673_GO_BP_up = c("positive regulation of cell migration", "epithelial cell migration", "collagen fibril organization", "positive regulation of cell motility", "integrin-mediated signaling pathway")

y_A673_GO_BP_up <- ggo_BP_A673_3days[ggo_BP_A673_3days$Description %in% vals_A673_GO_BP_up, asis = T]

pdf("/chk/CP-E_sarcoma_BETi/Cell_line_level_analyses/GSEA_analysis/20days_vs_3days_comparison/GSEA_all_genes_analysis/Results_20_vs_3days_BETi_DMSO/Plots/DotPlot_GO_BP_A673_upreg_select_v3.pdf")
dotplot_custom(y_A673_GO_BP_up, showCategory=10)
dev.off()

vals_CHLA10_GO_BP_up = c("axon development", "positive regulation of cell migration", "regulation of neuron projection development", "angiogenesis", "positive regulation of cell motility")

y_CHLA10_GO_BP_up <- ggo_BP_CHLA10_3days[ggo_BP_CHLA10_3days$Description %in% vals_CHLA10_GO_BP_up, asis = T]

pdf("/chk/CP-E_sarcoma_BETi/Cell_line_level_analyses/GSEA_analysis/20days_vs_3days_comparison/GSEA_all_genes_analysis/Results_20_vs_3days_BETi_DMSO/Plots/DotPlot_GO_BP_CHLA10_upreg_select_v3.pdf")
dotplot_custom(y_CHLA10_GO_BP_up, showCategory=10)
dev.off()

##################--------------------------##################################3

ggo_CC_A673_3days <- gseGO(gene = A673_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)
ggo_CC_CHLA10_3days <- gseGO(gene = CHLA10_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)

vals_A673_GO_CC_up = c("focal adhesion", "cell-substrate junction", "endoplasmic reticulum lumen", "collagen-containing extracellular matrix", "cell-cell junction")

y_A673_GO_CC_up <- ggo_CC_A673_3days[ggo_CC_A673_3days$Description %in% vals_A673_GO_CC_up, asis = T]

pdf("/chk/CP-E_sarcoma_BETi/Cell_line_level_analyses/GSEA_analysis/20days_vs_3days_comparison/GSEA_all_genes_analysis/Results_20_vs_3days_BETi_DMSO/Plots/DotPlot_GO_CC_A673_upreg_select_v3.pdf")
dotplot_custom(y_A673_GO_CC_up, showCategory=10)
dev.off()

vals_CHLA10_GO_CC_up = c("collagen-containing extracellular matrix", "basement membrane", "actin cytoskeleton", "cell-substrate junction", "endoplasmic reticulum lumen")

y_CHLA10_GO_CC_up <- ggo_CC_CHLA10_3days[ggo_CC_CHLA10_3days$Description %in% vals_CHLA10_GO_CC_up, asis = T]

pdf("/chk/CP-E_sarcoma_BETi/Cell_line_level_analyses/GSEA_analysis/20days_vs_3days_comparison/GSEA_all_genes_analysis/Results_20_vs_3days_BETi_DMSO/Plots/DotPlot_GO_CC_CHLA10_upreg_select_v3.pdf")
dotplot_custom(y_CHLA10_GO_CC_up, showCategory=100)
dev.off()

