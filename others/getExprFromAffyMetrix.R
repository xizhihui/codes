# @Author: XiZhihui
# @Date:   2018-10-12 21:58:14
# @Last Modified by:   XiZhihui
# @Last Modified time: 2018-10-12 22:15:19
# @Description: 本函数从AffyMetrix的.CEL源文件获得表达数据

getExprFromAffyMetrix = function(filepath, method='rma', destfile=NULL) {
	suppressPackageStartupMessages(library(affy))

	filenames = list.files(filepath, pattern='*.CEL')
	the.data = ReadAffy(filenames=file.path(filepath, filenames))
	old.names = sampleNames(the.data)
	sampleNames(the.data) = gsub('.CEL', '', filenames)
	pm.data = pm(the.data)
	mm.data = mm(the.data)
	pdata = pData(the.data)
	probes = geneNames(the.data)

	if (method == 'rma') {
		eSet = rma(the.data)
		eSet = 2 ^ (exprs(eSet))
	}
	if (method == 'mas') {
		eSet = mas5(this.data)
		eSet = exprs(eSet)
	}
	if (method == 'gcrma') {
		suppressPackageStartupMessages(library(gcrma))
		eSet = gcrma(the.data)
	}
	# mas5calls: P-present, A-absent, M-marginal
	calls = mas5calls(the.data)
	calls = exprs(calls)
	probes = apply(calls, 1, function(x) any(x == 'P'))

	expr = eSet[probes, ]
	if (destfile) save(expr, file=destfile)
	return(expr)
}