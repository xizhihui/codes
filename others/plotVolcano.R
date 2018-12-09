# @Author: XiZhihui
# @Date:   2018-10-12 20:01:44
# @Last Modified by:   XiZhihui
# @Last Modified time: 2018-10-12 20:23:23

# usage in R:
#	source('path/to/plotvolcano.R')
#	volcano(yourdataframe) 

# df: a data.frame with columns named 'log2FoldChange' and 'padj'
# 颜色区分：padj < 0.05 & log2FC > 1时是红， log2FC=1时是蓝，其余为灰色
plotvolcano = function(df, padj=0.05, log2FC=1) {
	df = as.data.frame(df)
	df$color = ifelse(
			df$padj < 0.05 & abs(df$log2FoldChange) >= 1,
			ifelse(df$log2FoldChange > 1, 'red', 'blue'),
			'gray')
	color = c(red='red', gray='gray', blue='blue')

	require(ggplot2)

	p = ggplot(df, aes(log2FoldChange, -log10(padj), color = color) +
		geom_point() +
	    theme_bw() +
	    scale_color_manual(values = color) +
	    labs(x="log2 (fold change)",y="-log10 (p-adj)") +
	    geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
	    geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
	    theme(legend.position = "none",
	          panel.grid=element_blank(),
	          axis.title = element_text(size = 16),
	          axis.text = element_text(size = 14))
	p
}