rowAUCs = function(expr, group) {
	options(stringsAsFactors = F)
	suppressPackageStartupMessages(require(pROC))
	expr = as.data.frame(expr)
	auc_per_gene = apply(expr, 1, function(x) auc(group, x))
	auc_per_gene = data.frame(auc=auc_per_gene)
	rownames(auc_per_gene) = colnames(expr)
	return(auc_per_gene)
}