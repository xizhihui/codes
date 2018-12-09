#######################
#  本代码根据从gdc下载的clinical_manifest文件将下载的clinical.xml文件进行整合,
#  整合结果就是一个比较全的临床信息矩阵
#  使用方法: getClinFromTCGA(manifest, data_path)
#  		manifest: TCGA下载的manifest文件
#  		data_path：使用gdc-client下载的数据存储路径
#  输出文件: clinical_matrix.txt
#######################

getClinFromTCGA = function(manifestFile, dataPath, destfile=NULL) {
	options(stringsAsFactors=F)

	# 由于临床数据是xml文件,加载相应的包
	
	suppressStartupMessages({
		require((XML)
		require(tidyr)
	})

	# manifest_path = commandArgs(T)[1]
	# path = commandArgs(T)[2]

	clinical_manifest = read.delim(manifest_path)

	ids = clinical_manifest[,1:2]
	clini_paths = paste(path, ids[,1], '/', ids[,2], sep='')

	n = length(clini_paths)
	df = data.frame()
	for (i in 1:n) {
	    directory = clini_paths[i]
	    fileID = ids$id[i]
	    # 读取xml文件,并生成data.frame
	    xmlfile = try(xmlParse(file=directory, encoding="UTF-8"), silent=T)
	    if('try-error' %in% class(xmlfile)) {
	        next
	    }
	    xmldf = as.data.frame(xmlToDataFrame(xmlfile))
	        
	    rownames(xmldf) = c('one', 'two')
	    xmldf = as.data.frame(t(xmldf))
	    result = unite(xmldf, 'infos', one, two)
	    result = as.data.frame(apply(result, 2, function(x) sub('(_NA)|(NA_)', '', x)))
	    colnames(result) = fileID
	    if (length(df) == 0) {
	        df = result
	    } else {
	        df = cbind(df, result)
	    }
	}

	colnames(df) = df['bcr_patient_barcode',]
	df = t(df)
	write.table(df,
	            file=ifelse(destfile, destfile, 'clinical_matrix.txt'),
	            col.names=T,
	            row.names=T)
	return(df)
}