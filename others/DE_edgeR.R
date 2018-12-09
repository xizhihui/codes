########################################################
# 使用edgeR进行差异分析
# 使用方式: 
#         source('DE_edgeR.R')
#         edgeR_DGE(exprSet, group, type)
# exprSet: 表达矩阵
# group: 样品分组
# type: 差异分析方法, classical, lrt, qlt
# + classcial
# + glm: likelihood ratio test/ quasi-likelihood F-test
# + + quasi-likelihood(qlf): 推荐用于差异分析，因为他对错误率限制较好。
# + + likelihood(lrt)：对与单细胞RNA-测序和没有重复的数据较好
##################### example ###########################
# suppressPackageStartupMessages(library(edgeR))
# setwd('~/practice/180716_edgeR/')

# rawdata = read.table('rawdata.txt')
# head(rawdata)
# rawdata = rawdata[-(1:5),]

# groups = grepl('01A', colnames(rawdata))
# groups = ifelse(groups, 'tumor', 'normal')
# table(groups)

# edgeR_DGE(exprSet = rawdata, group=groups, type='classical')
########################################################

edgeR_DGE <- function(exprSet, group, type, cpm=c(100,4)) {
    require(edgeR)
    # 生成DEGList
    DGE = DGEList(counts=exprSet, group=group)
    DGE.old = DGE
    
    # 生成design
    design = model.matrix(~group)
    cat('design生成完成\n')

    # cpm过滤
    keep = rowSums(edgeR::cpm(DGE) > cpm[1]) >= cpm[2]
    DGE = DGE[keep,]
    cat('过滤完成\n')

    # 进行校正
    DGE = calcNormFactors(DGE)
    cat('计算校正因子\n')

    # 检测离群值和关系
    png('plotMDS.png')
    plotMDS(DGE, method='bcv', col=as.numeric(DGE$samples$group))
    legendCol = unique(as.numeric(DGE$samples$group))
    legendGroup = unique(as.character(DGE$samples$group))
    legend("bottomright", legendGroup, col=legendCol, pch=20)
    dev.off()
    cat('MDS 出图完毕\n')

    if (type == 'classical') {
        # 计算离散度dispersion
        d = estimateCommonDisp(DGE)
        d = estimateTagwiseDisp(d)
        test = exactTest(d)
        cat('使用classical进行差异分析完毕\n')
    } else {
        # 计算离散度dispersion
        d = estimateGLMCommonDisp(DGE)
        d = estimateGLMTrendedDisp(d)
        d = estimateGLMTagwiseDisp(d)
        
        if (type == 'qlf') {
            fit = glmQLFit(d, design)
            test = glmQLFTest(fit, coef=2)
            cat('使用qlf差异分析完毕')
        } else if (type == 'lrt') {
            fit = glmFit(d, design)
            test = glmLRT(fit, coef=2)
            cat('使用lrt差异分析完毕')
        }
    }
    
    png('plotBCV.png')
    plotBCV(d)
    dev.off()
    
    png('plotSmear.png')
    de = decideTestsDGE(test, adjust.method="BH", p.value = 0.05)
    tags = rownames(d)[as.logical(de)]
    plotSmear(test, de.tags=tags)
    abline(h=c(-4,4), col='blue')
    dev.off()
    
    cat('BCV,Smear 出图完毕\n正则保存数据.')

    finalDGE = topTags(test, n=nrow(exprSet))
    finalDGE = as.data.frame(finalDGE)
    write.table(file='DGE_edgeR.txt', finalDGE)
    save(test, fit, file='DE_edgeR.RData')
    cat('数据保存完毕,运行完成')
}