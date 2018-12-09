#######################
#  本代码根据从gdc下载的RNAseq counts/FPKM_manifest文件将下载的数据文件进行整合,
#  整合结果就是某个癌症的表达矩阵文件
#  使用方法: getCountsFromTCGA(manifest, ID_transfer, data_path)
#       manifest: TCGA下载的manifest文件
#       ID_transfer: 下载的ID对照文件
#       data_path：使用gdc-client下载的数据存储路径
#  输出文件: counts_matrix.txt
#######################

# path.manifest = commandArgs(T)[1]
# path.ID = commandArgs(T)[2]
# path.data = commandArgs(T)[3]

getCountsFromTCGA = function(manifest, IDTransfer, dataPath, destfile=NULL) {
    options(stringsAsFactors=F)
    # 读取文件
    expr_manifest = read.delim(file=manifest)
    ID_transfer = read.delim(file=IDTransfer)
    path = dataPath

    # 获取所有FPKM文件的路径
    ids = expr_manifest[, c(1,2)]
    expr_paths = paste(path, '/', ids[,1], '/', ids[,2], sep='')

    # 读取所有的表达数据并合并
    df = data.frame()
    n = length(expr_paths)
    for (i in 1:n) {
        file_id = ids[i,1]
        mydata = read.table(gzfile(expr_paths[i]))
        names(mydata) = c('ENSG_ID', file_id)
        if (length(df) == 0) {
            df = mydata
        } else {
            df = merge(df, mydata, by='ENSG_ID')
        }
    }

    # 把file_id改成类似TCGA-FF-8061-01A的形式
    match_id = match(ids[,1], ID_transfer$id)
    match_case_id = ID_transfer$cases.0.samples.0.submitter_id[match_id]

    expr_matrix = df[,-1]
    rownames(expr_matrix) = df$ENSG_ID
    colnames(expr_matrix) = match_case_id

    # 存储备用
    write.table(file=ifelse(destfile, destfile, 'counts_matrix.txt'),
                expr_matrix,
                row.names = T,
                col.names = T)
    return(expr_matrix)
}