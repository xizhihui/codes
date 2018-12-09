library(limma)
library(edgeR)
library(DESeq2)

load('lncRNA_partial.txt')

groups = ifelse(grepl('11A', colnames(counts)), 'normal', 'tumor')
design = model.matrix(~factor(groups))
groups = data.frame(groups)
colnames(groups) = c('condition')
rownames(groups) = colnames(counts)

## limma 
{
    lm.voom = voom(counts, normalize.method = 'quantile', plot=T)
    lm.fit = lmFit(lm.voom, design = design)
    lmd.ebs = eBayes(lm.fit)
    top20.lm = topTable(lmd.ebs, coef = 2, adjust='BH', lfc=1, p.value=0.05, number=20)
}

## edgeR 
{
    ed.DE = DGEList(counts, group=groups$condition)
    ed.DE = calcNormFactors(ed.DE)
    ed.d = estimateGLMCommonDisp(ed.DE)
    ed.d = estimateGLMTrendedDisp(ed.d)
    ed.d = estimateGLMTagwiseDisp(ed.d)
    ed.fit = glmQLFit(ed.d, design = design)
    ed.test = glmQLFTest(ed.fit, coef=2)
    top20.ed = topTags(ed.test, sort.by='logFC', p.value = 0.05, n=20)
    top20.ed = as.data.frame(top20.ed)
}

## DESeq2
{
    ds.dds = DESeqDataSetFromMatrix(counts, colData=groups, design = ~condition)
    ds.dds = DESeq(ds.dds)
    ds.res = results(ds.dds, alpha = 0.05, lfcThreshold = 1)
    ds.res.data = as.data.frame(ds.res)
    top20.ds = ds.res.data[order(-abs(ds.res.data$log2FoldChange)),][1:20,]
}
