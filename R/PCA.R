# @Author: XiZhihui
# @Date:   2018-10-15 09:42:14
# @Last Modified by:   XiZhihui
# @Last Modified time: 2018-10-15 10:00:38

PCA = function (dat, prefix='PCA', output='.', plot3d=TRUE) {
	dat = t(as.matrix(dat))
	dat.class = rownames(dat)
	dat.pca = prcomp(dat, scale.=TRUE)
	pca.sum = summary(dat.pca)

	getOutput = function(name) {
		out_path = file.path(output, paste(prefix, name, sep='_'))
		out_path
	}
	# 输出特征向量
	write.table(dat.pca$rotation,
		file=getOutput('feature_vector.xls'),
		quote=F, sep='\t')

	# 输出新表
	write.table(predict(dat.pca), 
		file=getOutput('new_table.xls'),
		quote=F, sep='\t')

	# 输出解释变异比重
	write.table(pca.sum$importance,
		file=getOutput('new_table.xls'),
		quote=F, sep='\t')

	# 输出柱状图
	pdf(file=getOutput('barplot.pdf'))
	barplot(pca.sum$importance[2,1:10]*100,
		xlab="PC",
		ylab="percent",
		col="skyblue")
	dev.off()

	# 碎石图
	pdf(file=getOutput('gravel_map.pdf'))
	plot(pca.sum$importance[2,1:10]*100,
		type='o',
		xlab="PC",
		ylab="percent",
		col="skyblue")

	# 2d
	suppressPackageStartupMessages(require(ggplot2))
	pdf(file=getOutput('2d.pdf'))
	p = ggplot(data=dat.pca, aes(PCA1, PCA2)) + 
		geom_point(aes(color = group)) +
		geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=1, linetype=2) +
		annotate("text", x=PCA.mean$PCA1, y=PCA.mean$PCA2,label=PCA.mean$group) +
		theme_bw()+
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	dev.off()

	# 3d
	suppressPackageStartupMessages({
		require(pca3d)
		require(rgl)
	})
	pdf(file=getOutput('3d.pdf'))
	pca3d(data.pca, components = 1:3, group = c(rep("con",5),rep("A",5),rep("B",3)),show.centroids=TRUE,new=TRUE)
	dev.off()
}