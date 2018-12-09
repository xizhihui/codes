# knitr设置，出报告用的
## ----style, echo=FALSE, results='hide', message=FALSE--------------------
{
    library(knitr)
    # 错误信息和注释不输出，提示信息输出
    opts_chunk$set(error=FALSE, prompt=TRUE, comment=NA)
    # 设置图片宽高等
    opts_chunk$set(fig.width=7, fig.height=7, out.width="3.5in", fig.align="center", fig.path="")
    opts_chunk$set(dpi=300, dev="png", dev.args=list(pointsize=15))
    options(width=83)
}
    

## ----targets-------------------------------------------------------------
# 读取样品信息
targets <- read.delim("targets.txt", stringsAsFactors=FALSE)
targets

## ----group---------------------------------------------------------------
# 进行分组
Group <- factor(targets$Population)
Group

## ----readcounts----------------------------------------------------------
Sample <- targets$Sample
fn <- paste0(Sample,".bismark.cov.gz")
data <- list()
for(i in 1:length(Sample)) {
    data[[i]] <- read.delim(file=fn[i], header=FALSE)[,-(3:4)]
    names(data[[i]]) <- c("Chr", "Position", "Meth", "Un")
}

## ----combine-------------------------------------------------------------
position <- sapply(data, function(x) paste(x[,1], x[,2], sep="-") )
position_all <- unique(unlist(position))
counts <- matrix(0L, nrow=length(position_all), ncol=2*length(Sample))
for(i in 1:length(Sample)) {
    m <- match(position[[i]], position_all)
    counts[m, c(2*i-1,2*i)] <- as.matrix(data[[i]][, 3:4])
}

## ----counts--------------------------------------------------------------
rownames(counts) <- position_all
Sample2 <- rep(Sample, each=2)
Sample2 <- factor(Sample2)
Meth <- rep(c("Me","Un"), length(Sample))
Meth <- factor(Meth, levels=c("Un","Me"))
colnames(counts) <- paste(Sample2, Meth, sep="-")
head(counts)

## ----DGEList, message=FALSE----------------------------------------------
library(edgeR)
options(digits=3)
Chr <- gsub("-.*$", "", position_all)
Position <- gsub("^.*-", "", position_all)
Genes <- data.frame(Chr=Chr, Position=Position)
y <- DGEList(counts, genes=Genes, group=rep(Group,each=2))

## ----counts_total--------------------------------------------------------
counts_total <- t(rowsum(t(counts), Sample2)) # 存疑，rowsum(x,y)是什么
head(counts_total)

## ----keep----------------------------------------------------------------
keep <- rowSums(counts_total >= 10) == 6
table(keep)

## ----filter--------------------------------------------------------------
y <- y[keep, , keep.lib.sizes=FALSE]

## ----norm----------------------------------------------------------------
TotalReadCount <- colMeans(matrix(y$samples$lib.size, nrow=2, ncol=6))
y$samples$lib.size <- rep(TotalReadCount, each=2)
y$samples

## ----Mvalues-------------------------------------------------------------
Beta <- y$counts[, Meth=="Me"] / counts_total[keep, ]
logCPM <- cpm(y, log=TRUE, prior.count=2)   # counts-per-million
M <- logCPM[, Meth=="Me"] - logCPM[, Meth=="Un"]
colnames(Beta) <- colnames(M) <- Sample

## ----mdsplot, fig.width=12, fig.height=6, out.width="6in", fig.cap="The MDS plots of the methylation levels of the data set. Methylation levels are measured in beta values (left) and M-values (right). Samples are separated by the cell population in the first dimension in both MDS plots."----
par(mfrow=c(1,2))
plotMDS(Beta, col=rep(1:3, each=2), main="Beta-values")
plotMDS(M, col=rep(1:3, each=2), main="M-values")

## ----design--------------------------------------------------------------
design <- model.matrix(~ Sample2 + Meth)
colnames(design) <- gsub("Sample2","",colnames(design))
colnames(design) <- gsub("Meth","",colnames(design))
colnames(design)[1] <- "Int"
design <- cbind(design, 
    Me2=c(0,0,0,0,1,0,1,0,0,0,0,0),
    Me3=c(0,0,0,0,0,0,0,0,1,0,1,0))
design

## ----estimateDisp--------------------------------------------------------
y <- estimateDisp(y, design=design, trend="none")
y$common.dispersion
summary(y$prior.df)

## ----fit-----------------------------------------------------------------
fit <- glmFit(y, design)

## ----BvsL----------------------------------------------------------------
contr <- makeContrasts(BvsL=0.5*(Me2+Me3), levels=design)

## ----lrt-----------------------------------------------------------------
lrt <- glmLRT(fit, contrast=contr)

## ----topTags-------------------------------------------------------------
topTags(lrt)

## ----decideTests---------------------------------------------------------
summary(decideTests(lrt))

## ----plotMDfit, fig.cap="MD plot showing the log-fold change of the methylation level and average abundance of each CpG site. Significantly up and down methylated CpG's are highlighted in red and blue, respectively."----
plotMD(lrt)

## ----promoters, message=FALSE--------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
genes_Mm <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
pr <- promoters(genes_Mm, upstream=2000, downstream=1000)
pr

## ----sites---------------------------------------------------------------
Position <- as.numeric(Position)
sites <- GRanges(seqnames=Chr, ranges=IRanges(start=Position, end=Position))

## ----overlap-------------------------------------------------------------
olap <- findOverlaps(query=pr, subject=sites)
olap

## ----counts2-------------------------------------------------------------
counts2 <- counts[subjectHits(olap), ]
counts2 <- rowsum(counts2, queryHits(olap))

## ----annotation, message=FALSE-------------------------------------------
ind <- as.numeric(rownames(counts2))
rownames(counts2) <- pr$gene_id[ind]
library(org.Mm.eg.db)
anno <- select(org.Mm.eg.db, keys=pr$gene_id, columns="SYMBOL", 
               keytype="ENTREZID")
anno <- data.frame(Symbol=anno$SYMBOL[ind])
y2 <- DGEList(counts2, genes=anno, group=rep(Group,each=2))

## ----counts2_total-------------------------------------------------------
counts2_total <- t(rowsum(t(counts2), Sample2))

## ----filtering2----------------------------------------------------------
keep2 <- rowSums(counts2_total >= 10) == 6
table(keep2)
y2 <- y2[keep2,,keep.lib.sizes=FALSE]

## ----norm2---------------------------------------------------------------
TotalReadCount2 <- colMeans(matrix(y2$samples$lib.size, nrow=2, ncol=6))
y2$samples$lib.size <- rep(TotalReadCount2, each=2)
y2$samples

## ----mdsplot2, fig.width=12, fig.height=6, out.width="6in", fig.cap="The MDS plots of the methylation levels at gene promoters. Methylation levels are measured in beta values (left) and M-values (right). Samples are separated by the cell population in the first dimension in both MDS plots."----
Beta2 <- y2$counts[, Meth=="Me"] / counts2_total[keep2, ]
logCPM2 <- cpm(y2, log=TRUE, prior.count=2)
M2 <- logCPM2[, Meth=="Me"] - logCPM2[, Meth=="Un"]
colnames(Beta2) <- colnames(M2) <- Sample
par(mfrow=c(1,2))
plotMDS(Beta2, col=rep(1:3, each=2), main="Beta-values")
plotMDS(M2, col=rep(1:3, each=2), main="M-values")

## ----estimateDisp2-------------------------------------------------------
y2 <- estimateDisp(y2, design=design, trend="none", robust=TRUE)
y2$common.dispersion
summary(y2$prior.df)

## ----plotBCV2, width="3.8in", fig.cap="Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene. The plot shows the square-root estimates of the common and tagwise NB dispersions."----
plotBCV(y2)

## ----fit2----------------------------------------------------------------
fit2 <- glmFit(y2, design)

## ----lrt2----------------------------------------------------------------
lrt2 <- glmLRT(fit2, contrast=contr)

## ----topTags2------------------------------------------------------------
topTags(lrt2)

## ----decideTests2--------------------------------------------------------
summary(decideTests(lrt2))

## ----plotMDfit2, fig.cap="MD plot showing the log-fold change of the methylation level and average abundance of each CpG site. Significantly up and down methylated CpG's are highlighted in red and blue, respectively."----
plotMD(lrt2)

## ----Readin1-------------------------------------------------------------
load("rna.RData")
y_rna
rna_DE <- read.csv("BvsL-fc3.csv", row.names="GeneID")
head(rna_DE)

## ----lfc-----------------------------------------------------------------
tp <- topTags(lrt2, n=Inf, p=0.05)$table
m <- match(row.names(tp), row.names(rna_DE))
lfc <- tp[,c("Symbol","logFC")]
names(lfc)[2] <- "ME"
lfc$RNA <- rna_DE$logFC[m]
lfc <- lfc[!is.na(lfc$RNA), ]
head(lfc)

## ----Cor-----------------------------------------------------------------
cor(lfc$ME, lfc$RNA)

## ----Regression, fig.cap="Scatter plot of the log-fold changes of methylation levels in gene promoters (x-axis) vs the log fold-changes of gene expression (y-axis). The plot shows results for the genes of which the promoters are significantly differentially methylated between basal and luminal. The red line shows the least squares line with zero intercept. A strong negative correlation is observed."----
plot(lfc$ME, lfc$RNA, main="Basal vs Luminal", xlab="log-FC Methylation", 
    ylab="log-FC Gene Expression", pch=16, cex=0.8, col="gray30")
u <- lm(lfc$RNA ~ 0 + lfc$ME)
coef(summary(u))
abline(u, col="red", lwd=2)
abline(h=0, v=0, col="gray10", lty=2, lwd=2)

## ----Fry-----------------------------------------------------------------
ME <- data.frame(GeneID=row.names(lfc), weights=lfc$ME)
fry(y_rna, index=ME, design=y_rna$design, contrast=c(0,0,1,0,0,-1))

## ----Barcode, fig.width=8, fig.height=6.4, out.width="5.4in", out.height="4.3in", fig.cap="Barcode plot showing strong negative correlation between gene expression and DNA methylation in gene promoters."----
m <- match(row.names(rna_DE), row.names(tp))
gw <- tp$logFC[m]
gw[is.na(gw)] <- 0
barcodeplot(rna_DE$logFC, gene.weights=gw, labels=c("Luminal","Basal"), 
            main="Basal vs Luminal")
legend("topright", col=c("red","blue"), lty=1, lwd=2,
       legend=c("Up-methylation in Basal", "Up-methylation in Luminal") )

## ----info----------------------------------------------------------------
sessionInfo()
