library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gridExtra)
library(grid)
library(VennDiagram)
library(GenomicFeatures)

# Set working directory
setwd("D:/Lncrna diff analysis")

# Load and preprocess count data
countData <- read.table("SRRWITHGSM - Copy.tsv", header = TRUE, row.names = 1)
countDataFiltered <- countData[, !(colnames(countData) %in% c("Chr", "Start", "End", "Strand", "Length"))]

# Load lncRNA annotation (GTF) and extract lncRNA gene IDs
gtf <- makeTxDbFromGFF("D:/Lncrna diff analysis/gencode.v47.long_noncoding_RNAs.gtf", format="gtf")
lncrna_genes <- unique(genes(gtf)$gene_id)

# Filter lncRNA matrix
lncrna_matrix <- countDataFiltered[rownames(countDataFiltered) %in% lncrna_genes, ]
head(lncrna_matrix)

# Define sample metadata with MS classification stages
metadata <- data.frame(
  row.names = colnames(lncrna_matrix),
  sample_id = c(
    " sample1", "sample2", "sample3", "sample4",
    "sample5", "sample6", "sample7", "sample8"
  ),
  condition = c(
    "control", "control", "control", "control","RRMS", "RRMS","SPMS", "SPMS"
  )
)

# Create DESeq dataset with MS stages
dds <- DESeqDataSetFromMatrix(lncrna_matrix, metadata, design = ~condition)
head (ddsCollapsed)
# Collapse technical replicates
ddsCollapsed <- collapseReplicates(dds, groupby = metadata$sample_id)
head(ddsCollapsed)

# Differential expression analysis
ddsCollapsed <- DESeq(ddsCollapsed)

# Define comparisons
comparisons <- c("PPMS", "RRMS", "SPMS")
sig_gene_lists <- list()

# Loop through comparisons
for (comp in comparisons) {
  res <- results(ddsCollapsed, contrast = c("condition", comp, "control"))
  write.csv(as.data.frame(res), file = paste0("lncRNA_", comp, "_vs_Control_results.csv"))
  
  # Filter significant genes
  sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 2)
  sig_gene_lists[[comp]] <- rownames(sig_genes)
  write.csv(as.data.frame(sig_genes), file = paste0("lncRNA_", comp, "_vs_Control_significant_genes.csv"))
  
  # Generate and save Volcano plot
  png(filename = paste0("Volcano_lncRNA_", comp, "_vs_Control.png"), width=1200, height=1000, res=150)
  print(EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste("lncRNA", comp, "vs Control"),
                        pCutoff = 0.05, FCcutoff = 2))
  dev.off()
  
  # Generate and save Heatmap
  topGenes <- head(order(res$padj), 50)
  png(filename = paste0("Heatmap_top50_lncRNA_", comp, "_vs_Control.png"), width=1200, height=1200, res=150)
  pheatmap(assay(ddsCollapsed)[topGenes, ], scale="row", annotation_col=as.data.frame(colData(ddsCollapsed)))
  dev.off()
  
  # Generate and save top 50 genes table image
  table_top50 <- as.data.frame(res[topGenes,])
  png(filename = paste0("Table_top50_lncRNA_", comp, "_vs_Control.png"), width=1400, height=2000, res=150)
  grid.table(table_top50)
  dev.off()
}

# PCA plot
png("PCA_lncRNA.png", width=1200, height=1000, res=150)
plotPCA(vst(ddsCollapsed, blind=FALSE), intgroup="condition")
dev.off()

# Venn Diagram for overlapping significant lncRNAs
venn.diagram(sig_gene_lists, filename="Venn_overlapping_lncRNA.png", fill=c("blue", "red", "green"),
             category.names=comparisons, output=TRUE)
