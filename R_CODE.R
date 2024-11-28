library(DESeq2)
library(vsn)
library(hexbin)
library(devEMF)
library(limma)

#FUNCTIONS
plotPCAWithSampleNames = function(x, targets=targets, intgroup=colnames(targets)[1], ntop=500)
{
  library(RColorBrewer)
  library(genefilter)
  library(lattice)
  
  # pca
  #rv = rowVars(assay(x))
  rv = rowVars(x)
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  #pca = prcomp(t(assay(x)[select,]))
  pca = prcomp(t(x[select,]))
  
  # proportion of variance
  variance = pca$sdev^2 / sum(pca$sdev^2)
  variance = round(variance, 3) * 100
  
  # sample names
  names = colnames(x)
  #names = as.character(x$sample)
  
  # factor of groups
  #fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  fac = factor(apply(as.data.frame(targets[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  
  # colors
  if( nlevels(fac) >= 10 )
    colors = rainbow(nlevels(fac))
  else if( nlevels(fac) >= 3 )
    colors = brewer.pal(nlevels(fac), "Set1")
  else
    colors = c( "dodgerblue3", "firebrick3" )
  
  # plot
  xyplot(
    PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=1.5,
    aspect = "fill",
    col = colors,
    xlab = list(paste("PC1 (", variance[1], "%)", sep=""), cex=0.8),
    ylab = list(paste("PC2 (", variance[2], "%)", sep=""), cex=0.8),
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...);
      ltext(x=x, y=y, labels=names, pos=1, offset=0.8, cex=0.7)
    },
    main = draw.key(
      key = list(
        rect = list(col = colors),
        text = list(levels(fac)),
        rep = FALSE
      )
    )
  )
}

gene_count= read.csv("gene_count_matrix.csv", row.names=1)
transcript_count= read.csv("transcript_count_matrix.csv", row.names=1)
design_exp= read.csv("assess_design.csv", row.names=1)


dds_gene= DESeqDataSetFromMatrix(countData = gene_count,
                                   colData = design_exp,
                                   design = ~ batch + group)

dds_transcript <- DESeqDataSetFromMatrix(countData = transcript_count,
                                         colData = design_exp,
                                         design = ~ batch + group)


# Perform pre-filtering for gene counts
keep_gene=rowSums(counts(dds_gene) >= 10) >= 3
dds_gene=dds_gene[keep_gene,]

# Perform pre-filtering for transcript counts
keep_transcript=rowSums(counts(dds_transcript) >= 10) >= 3
dds_transcript=dds_transcript[keep_transcript,]

# Run DE analysis for gene counts
dds_gene=DESeq(dds_gene)

# Run DE analysis for transcript counts
dds_transcript=DESeq(dds_transcript)

# Make dispersion plot for gene counts
plotDispEsts(dds_gene)

# Make dispersion plot for transcript counts
plotDispEsts(dds_transcript)

# Perform rlog transformation for gene counts
rld_gene=rlog(dds_gene)

# Perform rlog transformation for transcript counts
rld_transcript=rlog(dds_transcript)

# Create PCA plot for gene counts
plotPCAWithSampleNames(assay(rld_gene),targets=design_exp, intgroup="group")

plotPCAWithSampleNames(assay(rld_transcript),targets=design_exp, intgroup="group")

# Calculate normalized counts
norm_counts_gene=counts(dds_gene, normalized = TRUE)
norm_counts_transcript=counts(dds_transcript, normalized = TRUE)

# Log-normalized counts
log_norm_counts_gene=log2(norm_counts_gene + 1)  # Adding 1 to avoid log(0)
log_norm_counts_transcript=log2(norm_counts_transcript + 1)

# rlog-transformed counts
rlog_counts_gene=assay(rlog(dds_gene))
rlog_counts_transcript=assay(rlog(dds_transcript))


meanSdPlot(log_norm_counts_gene, ylim = c(0, 3), main = "Log-normalized Gene Counts")
meanSdPlot(rlog_counts_gene, ylim = c(0, 3), main = "rlog-transformed Gene Counts")

# Extract DE results for B vs A, and C vs. A contrasts using a standard function that returns unshrunken LFCs
# For standard NULL hypothesis (LFC=0)
results_BvsA <- results(dds_gene, contrast=c("group", "B", "A"), alpha=0.05, lfcThreshold = 0)
results_CvsA <- results(dds_gene, contrast=c("group", "C", "A"), alpha=0.05, lfcThreshold = 0)
summary(results_BvsA,alpha=0.05)
summary(results_CvsA,alpha=0.05)

# For NULL hypothesis of LFC<1
results_BvsA_LFC1 <- results(dds_gene, contrast=c("group", "B", "A"), alpha=0.05, lfcThreshold = 1)
results_CvsA_LFC1 <- results(dds_gene, contrast=c("group", "C", "A"), alpha=0.05, lfcThreshold = 1)
summary(results_BvsA_LFC1,alpha=0.05)
summary(results_CvsA_LFC1,alpha=0.05)

#MA plots for unshrunk LFCs
DESeq2::plotMA(results_BvsA, ylim= c(-4,4), main="B vs A (LFC=0)")
DESeq2::plotMA(results_BvsA_LFC1, ylim= c(-4,4), main="B vs A (LFC<1)")
DESeq2::plotMA(results_CvsA, ylim= c(-4,4), main="C vs A (LFC=0)")
DESeq2::plotMA(results_CvsA_LFC1, ylim= c(-4,4), main="C vs A (LFC<1)")


# Extract DE results with shrunken LFCs for standard NULL hypothesis (LFC=0)
res_shrink_BvsA<- lfcShrink(dds_gene, contrast=c("group", "B", "A"), type="ashr", lfcThreshold=0)
res_shrink_CvsA <- lfcShrink(dds_gene, contrast=c("group", "C", "A"), type="ashr", lfcThreshold=0)
summary(res_shrink_BvsA,alpha=0.05)
summary(res_shrink_CvsA,alpha=0.05)


# Extract DE results with shrunken LFCs for NULL hypothesis of LFC<1
res_shrink_BvsA_LFC1 <- lfcShrink(dds_gene, contrast=c("group", "B", "A"), type="ashr", lfcThreshold=1)
res_shrink_CvsA_LFC1 <- lfcShrink(dds_gene, contrast=c("group", "C", "A"), type="ashr", lfcThreshold=1)
summary(res_shrink_BvsA_LFC1,alpha=0.05)
summary(res_shrink_CvsA_LFC1,alpha=0.05)

#Getting significant genes from DE results with shrunken LFC and NULL hypothesis of LFC < 1
# Filter significant genes for B vs A contrast
sig_genes_BvsA <- subset(res_shrink_BvsA_LFC1, padj < 0.05)
sig_genes_BvsA <- sig_genes_BvsA[order(sig_genes_BvsA$log2FoldChange), ]
write.csv(sig_genes_BvsA, file = "Significant.Genes.BvsA.csv", row.names = TRUE)

# Filter significant genes for C vs A contrast
sig_genes_CvsA <- subset(res_shrink_CvsA_LFC1, padj < 0.05)
sig_genes_CvsA <- sig_genes_CvsA[order(sig_genes_CvsA$log2FoldChange), ]
write.csv(sig_genes_CvsA, file = "Significant.Genes.CvsA.csv", row.names = TRUE)

#MA plots for shrunken LFCs
DESeq2::plotMA(res_shrink_BvsA, ylim= c(-4,4), alpha=0.05, main="B vs A shrunk (LFC=0)")
DESeq2::plotMA(res_shrink_BvsA_LFC1, ylim= c(-4,4), alpha= 0.05, main="B vs A  shrunk (LFC<1)")
DESeq2::plotMA(res_shrink_CvsA, ylim= c(-4,4), alpha= 0.05, main="C vs A shrunk (LFC=0)")
DESeq2::plotMA(res_shrink_CvsA_LFC1, ylim= c(-4,4), alpha=0.05,  main="C vs A shrunk (LFC<1)")

#Batch correction
dds_no_batch= DESeqDataSetFromMatrix(countData=gene_count,colData=design_exp, design= ~ group)
keep_no_batch=rowSums(counts(dds_no_batch) >= 10) >= 3

dds_no_batch= dds_no_batch[keep_no_batch,]
dds_no_batch= DESeq(dds_no_batch)

rld_no_batch= rlog(dds_no_batch, blind=TRUE)
plotPCAWithSampleNames(assay(rld_no_batch),targets=design_exp, intgroup='group')

mydesign= model.matrix(design(dds_no_batch),colData(dds_no_batch))
b.corrected= limma::removeBatchEffect(assay(rld_no_batch), batch=colData(dds_gene)$batch, design=mydesign)
plotPCAWithSampleNames(b.corrected,targets=design_exp, intgroup='group')
write.csv(b.corrected,file="BatchCorrected.Rlog.csv",quote=FALSE)




