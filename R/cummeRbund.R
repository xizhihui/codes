library(cummeRbund)

savepng = function(file, obj) {
	png(file)
	print(obj)
	dev.off()
}

path = "path/to/your/cufflinks/output/"

cuff = readCufflinks(dir=path)

######### informations and assesstion about data ################
runinfo(cuff)
replicates(cuff)
gene.features = annotation(genes(cuff))
gene.fpkm = fpkm(genes(cuff))
gene.repfpkm = repFpkm(genes(cuff))
gene.count = count(genes(cuff))
isoforms.fpkm = fpkm(isoforms(cuff))
gene.diff = diffData(genes(cuff))
# samples(), featureNames(), fpkmMatrix(), repFpkmMatrix()
# countMatrix()

######### quality assessment of data ################
# dispersion plot
d = dispersionPlot(genes(cuff))
d
# boxplot, dendroplot (JS distance), FPKM distributions, outlier
pBoxRep = csBoxplot(genes(cuff), replicates=T)
pBoxRep
pBox = csBoxplot(genes(cuff))
pBoxs
pDendro = csDendro(genes(cuff), replicates=T)
pDendro
# fpkmSCV plot at different level, cross-replicate variability between conditions
genes.scv = fpkmSCVPlot(genes(cuff))
isoforms.scv = fpkmSCVPlot(isoforms(cuff))
# csDensity plot, distributions of FPKM scores accross samples
csdens = csDensity(genes(cuff))
csdens
csdensRep = csDensity(genes(cuff), replicates=T)
csdensRep
# pairwise scatterplots matrix, 不比较矩阵，只旋转2个样本的话，使用csScatter
# identify global changes and trends in gene expression between pairs of conditions
csscatter = csScatterMatrix(genes(cuff))
csscatterSingle = csScatter(genes(cuff), 'sample1', 'sample2', smooth=T)
csscatter; csscatterSingle;
# MA plot, determine any systematic bias that may be present between conditions
m = MAplot(genes(cuff), 'sample1', 'sample2')
mcount = MAplot(genes(cuff), 'sample1', 'sample2', useCount=T)
m; mcount;
# Volcano plot
v = csVolcanoMatrix(genes(cuff))
#v <- csVolcano(genes(cuff),"hESC","Fibroblasts")
v

########### differential expression gene set ###############
#### 获取目标基因集，比如差异基因
# Sig gene id
sigGeneIds = getSig(cuff,alpha=0.05,level="genes")
hESCvsFibroblast.sigGeneIds<-getSig(cuff,"hESC","Fibroblasts",alpha=0.05,level="genes")
# Sig gene set, sigGeneIds你可以自己指定
sigGenes = getGenes(cuff, sigGeneIds)
# fpkm(sigGenes); fpkm(isoforms(sigGenes));
# 获取某些样本的
sigGenes.plury = getGenes(cuff, sigGeneIds, sampleIdList=c("hESC","Fibroblasts"))

#### 可视化目标基因集
# heatmap of FPKM
h = csHeatmap(sigGenes, cluster="both")
h.rep = csHeatmap(sigGenes, cluster="both", replicates=T)
h; h.rep;
isoh = csHeatmap(isoforms(sigGenes), cluster='both', labRow=F)
tsoh = csHeatmap(TSS(sigGenes), cluster='both', labRow=F)
# barplot of FPKM
b = expressionBarplot(sigGenes)
b
# comparisons between two conditions
s = csScatter(sigGenes, "hESC","Fibroblasts", smooth=T)
s
# volcano
v = csVolcano(sigGenes, "hESC","Fibroblasts")
v


########### differential expression individual gne ###############
# 当sigGeneId 只含一个基因时。
mygene = getGenes(sigGenes, 'PINK1')
# 可以像上面一样，画单个基因的表达量图
g1 = expressionPlot(mygene)
g2 = expressionBarplot(mygene, replicates=T)
g3 = expressionBarplot(isoforms(mygene))
g4 = csPie(mygene, level='isoforms')
# 基因结构可视化
genetrack = makeGeneRegionTrack(mygene)
plotTracks(genetrack)

########### data exploration ###############
# 有多个处理条件时，获取各个condition下的显著基因数及交集，类似venn图
mySigMat = sigMatrix(cuff, level='genes', alpha=0.05)
# 获取其中2个处理的显著基因Id集
hESC_vs_iPS.sigIsoformIds<-getSig(cuff,x='hESC',y='iPS',alpha=0.05,level='isoforms')
mySigTable = getSigTable(cuff, alpha=0.01, level='genes')
# 相似性分析
myDistHeat = csDistHeat(genes(cuff))
# MDS分析，PCA分析, NMF
genes.PCA = PCAplot(genes(cuff), "PC1", "PC2")
genes.MDS = MDSplot(genes(cuff))
genes.nmf = csNMF(genes(cuff))
# genes.PCA.rep = PCAplot(genes(cuff), "PC1", "PC2", replicates=T)
genes.PCA; genes.MDS; genes.nmf;
# PAM聚类分析
ic = csCluster(mygene, k=4)
icp = csClusterPlot(ic)
icp
# condition-specific，当值为1时，该feature只在值为1的那个feature表达
mygene.spec = csSpecificity(mygene)
aimgene = rownames(mygene.spec[mygene.spec == 1])
aimgeneSimilar = findSimilar(cuff, aimgene[0], n=20)
aimgeneSimilar.expr = expressionPlot(aimgeneSimilar, logMode=T, showErrorbars=F)
