
# Attach the library
library(clusterProfiler)
# Package that contains MSigDB gene sets in tidy format
library(msigdbr)
# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)
# We will need this so we can use the pipe: %>%
library(magrittr)
library(dplyr)
library(ggplot2)

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

#generatio_transformed <- sapply(data_df$generatio, function(x) eval(parse(text=x)))
#data_df$GeneRatio <- generatio_transformed

ggplot(data = data_df, aes(x = log10_padj, y = reorder(Description,log10_padj),
			                           color = GeneRatio, size = NES)) +
  geom_point() +
    scale_color_gradient(name = "GeneRatio", low = "red", high = "blue") +
      theme_bw() +
        ylab("") +
	  xlab("-log10(p.adjust)") +
	    ggtitle("Enrichment analysis") + expand_limits(x=0, y=0)
}

A673_20_vs_3days <- read.csv("../tables/A673_condition_treated_results.csv", row.names=1)
CHLA10_20_vs_3days <- read.csv("/tables/CHLA10_condition_treated_results.csv",row.names=1)

A673_20_vs_3days <- A673_20_vs_3days[order(A673_20_vs_3days$log2FoldChange, decreasing = TRUE),]
CHLA10_20_vs_3days <- CHLA10_20_vs_3days[order(CHLA10_20_vs_3days$log2FoldChange, decreasing = TRUE),]

#A673_20_vs_3days <- filter(A673_20_vs_3days, padj < 0.05)
#CHLA10_20_vs_3days <- filter(CHLA10_20_vs_3days, padj < 0.05)

A673_20_vs_3days$Gene <- rownames(A673_20_vs_3days)
CHLA10_20_vs_3days$Gene <- rownames(CHLA10_20_vs_3days)

# create a mapped data frame we can join to the differential.
A673_20_vs_3days_mapped_df <- data.frame(
  gene_id = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = A673_20_vs_3days$Gene,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene id, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_id)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(A673_20_vs_3days, by = c("Ensembl" = "Gene"))

# create a mapped data frame we can join to the differential.
CHLA10_20_vs_3days_mapped_df <- data.frame(
  gene_id = mapIds(
    # Replace with annotation package for the organism relevant to your data
    org.Hs.eg.db,
    keys = CHLA10_20_vs_3days$Gene,
    # Replace with the type of gene identifiers in your data
    keytype = "ENSEMBL",
    # Replace with the type of gene identifiers you would like to map to
    column = "ENTREZID",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a gene id, drop that
  # from the data frame
  dplyr::filter(!is.na(gene_id)) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(CHLA10_20_vs_3days, by = c("Ensembl" = "Gene"))


A673_20_vs_3days_mapped_df <- A673_20_vs_3days_mapped_df %>%
  # Sort so that the highest absolute values of the log2 fold change are at the top
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(gene_id, .keep_all = TRUE)

CHLA10_20_vs_3days_mapped_df <- CHLA10_20_vs_3days_mapped_df %>%
  # Sort so that the highest absolute values of the log2 fold change are at the top
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(gene_id, .keep_all = TRUE)

# Let's create a named vector ranked based on the log2 fold change values
A673_lfc_vector <- A673_20_vs_3days_mapped_df$log2FoldChange
names(A673_lfc_vector) <- A673_20_vs_3days_mapped_df$gene_id

# Let's create a named vector ranked based on the log2 fold change values
CHLA10_lfc_vector <- CHLA10_20_vs_3days_mapped_df$log2FoldChange
names(CHLA10_lfc_vector) <- CHLA10_20_vs_3days_mapped_df$gene_id

# We need to sort the log2 fold change values in descending order here
A673_genelist <- sort(A673_lfc_vector, decreasing = TRUE)
CHLA10_genelist <- sort(CHLA10_lfc_vector, decreasing = TRUE)

library(enrichplot)

print(head(A673_genelist))
print(head(CHLA10_genelist))

A673_H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

A673_H_gsea <- GSEA(A673_genelist, TERM2GENE = A673_H_t2g)
print(head(A673_H_gsea))

CHLA10_H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

CHLA10_H_gsea <- GSEA(CHLA10_genelist, TERM2GENE = CHLA10_H_t2g)
print(head(CHLA10_H_gsea))

#write.table(A673_H_gsea, "MSigdb_A673_hallmark_GSEA.txt", sep="\t", quote=FALSE)
#write.table(CHLA10_H_gsea, "MSigdb_CHLA10_hallmark_GSEA.txt", sep="\t", quote=FALSE)

#pdf("./DotPlot_MSigdb_hallmark_A673.pdf")
#dotplot(A673_H_gsea, showCategory=20)
#dev.off()

#pdf("./DotPlot_MSigdb_hallmark_CHLA10.pdf")
#dotplot(CHLA10_H_gsea, showCategory=20)
#dev.off()

########################3---------------------###################

vals_A673_H_gsea_up = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_DN", "HALLMARK_APICAL_JUNCTION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL2_STAT5_SIGNALING")

y_A673_H_gsea_up <- A673_H_gsea[A673_H_gsea$Description %in% vals_A673_H_gsea_up, asis = T]

pdf("DotPlot_A673_H_gsea_up_upreg.pdf", width=10)
dotplot_custom(y_A673_H_gsea_up, showCategory=10)
dev.off()

vals_CHLA10_H_gsea_up = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_DN", "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_ESTROGEN_RESPONSE_EARLY")

y_CHLA10_H_gsea_up <- CHLA10_H_gsea[CHLA10_H_gsea$Description %in% vals_CHLA10_H_gsea_up, asis = T]

pdf("DotPlot_CHLA10_H_gsea_up.pdf", width=10)
dotplot_custom(y_CHLA10_H_gsea_up, showCategory=10)
dev.off()


