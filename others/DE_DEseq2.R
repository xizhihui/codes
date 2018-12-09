###########################################################
# countData: 以样品为列名,基因为行名的表达矩阵
# colData的形式: 以样品作为行名,列是样品对应的分组类型
#
# cts是count matrix
##########################################################

# expr_raw = read.table('lncRNA_matrix.txt')
# # 去除在10%的样本中counts数低于1的
# expr = expr_raw[apply(expr_raw, 1, function(x) sum(x >1) > length(x)*0.1), ]

# groups = grepl(pattern = '11A', x = colnames(expr))
# groups = ifelse(groups, 'tumor', 'normal')

# coldata = data.frame(groups)
# rownames(coldata) = colnames(expr)
# colnames(coldata) = c('condition')

DE_DEseq = function(expr, coldata, condition, padj=0.1, output='.', prefix='DESeq2', gene=NULL) {
	require(DESeq2)
	# 导入矩阵
	dds <- DESeqDataSetFromMatrix(countData = expr,
	                              colData = coldata,
	                              design = condition)
	# 生成DESeq数据集
	dds <- DESeq(dds)
	# 剔除低counts基因
	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep,]
	# 设定条件因子,抛弃无元素因子
	dds$condition <- droplevels(dds$condition)

	dds <- DESeq(dds)
	res <- results(dds)
	resOrdered <- res[order(res$pvalue),]
	write.table(resOrdered[res$padj < padj, ], file=file.path(output, paste0(prefix, '_expr_matrix_padj.txt')))
	write.table(resOrdered, file=file.path(output, paste0(prefix, '_expr_matrix_all.txt'))
	
	# visualization
	# 数据校正后的箱图
	log10_dds = log10(assay(dds))[['cooks']]
	colnames(log10_dds) = 1:length(colnames(log10_dds))
	boxplot(log10_dds, range=0, las=2, xlab='')

	# volcano, MA-plot
	resLFC = lfcShrink(dds)
	pdf(file=file.path(output, past0(prefix, '_MA_plot.pdf')))
	plotMA(resLFC, ylim=c(-2,2))
	dev.off()

	# pca
	dds_vst = vst(dds)
	pdf(file=file.path(output, paste0(prefix, '_PCA_plot.pdf')))
	plotPCA((dds_vst, intgroup=c('condition')))
	dev.off()

	# Distribution of counts
	pdf(file=file.path(output, paste0(prefix, '_counts_distribution_with_minpadj.pdf')))
	plotCounts(dds, gene=ifelse(gene, gene, which.min(resOrdered$padj)), intgroup='condition')

	# 表达矩阵热图
	require(pheatmap)
	select = order(rowMeans(counts(dds, normalized=T)), decreasing=T)[1:20]
	df <- as.data.frame(colData(dds)[,c('condition')])
	pheatmap(assay(dds_vst)[select,], 
			cluster_rows=FALSE, 
			show_rownames=FALSE, 
			cluster_cols=FALSE, 
			annotation_col=df)

	# 样品距离热图
	sampleDists = dist(t(assay(dds_vst)))
	require(RColorBrewer)
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep='-')
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)
	pheatmap(sampleDistMatrix,
			clustering_distance_rows=sampleDists,
			clustering_distance_cols=sampleDists,
			col=colors)
}