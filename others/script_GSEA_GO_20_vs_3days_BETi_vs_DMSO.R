
# Attach the library
library(clusterProfiler)
# Package that contains MSigDB gene sets in tidy format
library(msigdb)
# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(dplyr)

dir.create("Results_20_vs_3days_BETi_DMSO")

A673_3days <- read.csv("/active/taylor_s/people/nkatiy/CP-E_sarcoma_BETi/lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Results_A673_BETi_20_vs_3days/A673_condition_treated_results.csv", row.names=1)
CHLA10_3days <- read.csv("/active/taylor_s/people/nkatiy/CP-E_sarcoma_BETi/lawlorlab_shireen_2023.07.12_bulk_rnaseq_count_nf/additional_results/Results_CHLA10_BETi_20_vs_3days/CHLA10_condition_treated_results.csv",row.names=1)

#A673_3days <- filter(A673_3days, padj < 0.05)
#CHLA10_3days <- filter(CHLA10_3days, padj < 0.05)

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

write.table(ggo_BP_A673_3days, "./Results_20_vs_3days_BETi_DMSO/GO_BP_A673_3days_GSEA.txt", sep="\t", quote=FALSE)
write.table(ggo_BP_CHLA10_3days, "./Results_20_vs_3days_BETi_DMSO/GO_BP_CHLA10_3days_GSEA.txt", sep="\t", quote=FALSE)

library(enrichplot)

pdf("./Results_20_vs_3days_BETi_DMSO/DotPlot_GO_BP_A673.pdf")
dotplot(ggo_BP_A673_3days, showCategory=20)
dev.off()

pdf("./Results_20_vs_3days_BETi_DMSO/DotPlot_GO_BP_CHLA10.pdf")
dotplot(ggo_BP_CHLA10_3days, showCategory=20)
dev.off()

ggo_CC_A673_3days <- gseGO(gene = A673_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)
ggo_CC_CHLA10_3days <- gseGO(gene = CHLA10_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)

write.table(ggo_CC_A673_3days, "./Results_20_vs_3days_BETi_DMSO/GO_CC_A673_3days_GSEA.txt", sep="\t", quote=FALSE)
write.table(ggo_CC_CHLA10_3days, "./Results_20_vs_3days_BETi_DMSO/GO_CC_CHLA10_3days_GSEA.txt", sep="\t", quote=FALSE)

pdf("./Results_20_vs_3days_BETi_DMSO/DotPlot_GO_CC_A673.pdf")
dotplot(ggo_CC_A673_3days, showCategory=20)
dev.off()

pdf("./Results_20_vs_3days_BETi_DMSO/DotPlot_GO_CC_CHLA10.pdf")
dotplot(ggo_CC_CHLA10_3days, showCategory=20)
dev.off()

ggo_MF_A673_3days <- gseGO(gene = A673_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)
ggo_MF_CHLA10_3days <- gseGO(gene = CHLA10_genelist, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", minGSSize = 10, pvalueCutoff = 0.05, verbose = FALSE)

write.table(ggo_MF_A673_3days, "./Results_20_vs_3days_BETi_DMSO/GO_MF_A673_3days_GSEA.txt", sep="\t", quote=FALSE)
write.table(ggo_MF_CHLA10_3days, "./Results_20_vs_3days_BETi_DMSO/GO_MF_CHLA10_3days_GSEA.txt", sep="\t", quote=FALSE)

pdf("./Results_20_vs_3days_BETi_DMSO/DotPlot_GO_MF_A673.pdf")
dotplot(ggo_MF_A673_3days, showCategory=20)
dev.off()

pdf("./Results_20_vs_3days_BETi_DMSO/DotPlot_GO_MF_CHLA10.pdf")
dotplot(ggo_MF_CHLA10_3days, showCategory=20)
dev.off()

