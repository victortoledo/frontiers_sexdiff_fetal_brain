#------------------------------------------------------------------------
#title: "DSEq2 for expression values"
#author: "Victor Toledo victor.toledo.btos@gmail.com"
#date: "10/08/2021"
#------------------------------------------------------------------------


################### Load all the packages and set the folder
setwd("~/Documents/projects/ms-victor/sex_differences_fetal_brain/DESeq2/_m/")
options(stringsAsFactors = FALSE)
library(DESeq2)
library(calibrate)
library(gplots)
library(RColorBrewer)
######################### Read files
countdata = read.table('~/Documents/projects/ms-victor/sex_differences_fetal_brain/cemitool/_h/corrected_exp.txt', header=TRUE, row.names=1, check.names = FALSE)
countdata = as.matrix(countdata)
head(countdata)

######################### Experimental conditions and number of replicates
condition = factor(c(rep("female",50), rep("male",69)))
coldata = data.frame(row.names=colnames(countdata), condition)

######################### Run DESeq2
dds = DESeqDataSetFromMatrix(countData=round(countdata), colData=coldata, design=~condition)
dds = DESeq(dds)

######################### Plot dispersion
par(mar=c(1,1,1,1))
pdf("dispersion_plot.pdf")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

######################### Regularized log transform for heatmaps/clustering plots
rld = rlogTransformation(dds)
head(assay(rld))

pdf("hist.pdf")
hist(assay(rld))
dev.off()

########### Heatmap
mycols = brewer.pal(8, "Dark2")[1:length(unique(condition))]
sampleDists = as.matrix(dist(t(assay(rld))))
pdf("heatmap.pdf")
heatmap(as.matrix(sampleDists), key=F, trace="none",
        col=colorpanel(100, "black", "white"),
        ColSideColors=mycols[condition], RowSideColor=mycols[condition],
        margin=c(10,10), main="Sample Distance Matrix")
dev.off()

######################### Get the differential expression results
res = results(dds)
table(res$padj<0.05)
res = res[order(res$padj),]
resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] = "Gene"
head(resdata)
write.csv(resdata, file="~/Documents/projects/ms-victor/sex_differences_fetal_brain/DESeq2/_m/deseq2_results.csv")

######################### MA plot

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
pdf("MA_plot.pdf")
maplot(resdata, main="MA Plot")
dev.off()

######################### Volcano plot

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
pdf("volcano_plot.pdf")
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8)
dev.off()
