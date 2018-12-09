
########### 这些都是必须列 ###############
samples <- dplyr::transmute(mocks,
				SampleID=paste(sample, unit, sep="."),
				Replicate=unit,
				Factor=sample,
				Tissue=ifelse(is.null(tissue), rep(NA, length(sample)), tissue),
				Condition=ifelse(is.null(condition), rep(NA, length(sample)), condition),
				ControlID=paste(mock_sample, mock_unit, sep="-"),
				bamReads=add_bams(sample, unit),
				bamControl=add_bams(mock_sample, mock_unit),
				Peaks=add_Peaks(type, sample, unit))

# total experiment
experiment = ChIPQC(samples)
experiment
ChIPQCreport(experiment)

# single sample
sample = ChIPQCsample("chip.bam")
ChIPQCreport(sample)

#################### encode example ####################
# 1. read metadata
samples = read.csv(file.path(system.file("extdata", package="ChIPQC"),
					"example_QCexperiment.csv"))
# SampleID Tissue Factor Replicate bamReads	Peaks
# CTCF_1 A549 CTCF 1 reads/SRR568129.bam peaks/SRR568129_chr22_peaks.bed
# CTCF_2 A549 CTCF 2 reads/SRR568130.bam peaks/SRR568130_chr22_peaks.bed
# cMYC_1 A549 cMYC 1 reads/SRR568131.bam peaks/SRR568131_chr22_peaks.bed
# cMYC_2 A549 cMYC 2 reads/SRR568132.bam peaks/SRR568132_chr22_peaks.bed

# 2. ChIPQCexperiment object
#	BiocParallel is used inside
exampleExp = ChIPQC(samples, annotaiton="hg19")
singleExp = ChIPQCsample(sample, 
						annotaiton="hg19",	# optional, genome assembly
						chromosomes=18,		# optional, a vector?
						mapQCth=15,			# default 15, MAPQ value
						blacklist=blacklist,	# a file (bed?) or Granges object
						)
# check the quality control metrics
# QCmetrics(exampleExp)

# 3. get a summary QC reprot 
ChIPQCreport(exampleExp, facetBy=c("Tissue","Factor"))	


######################### plotting QCmetrics for experimental sample groups ####################
##### 这个是自己画图, 以下都是返回的 ggplot 对象，意味着什么你懂吧
# coverage histogram
plotCoverageHist(tamoxifen,facetBy=c("Tissue","Factor"))
# cross-coverage plot
plotCC(tamoxifen,facetBy=c("Tissue","Factor"))
# Plotting Relative Enrichment of reads in Genomic Intervals
plotRegi(tamoxifen,facetBy=c("Tissue","Condition"))
# Plots of composite peak profile
plotPeakProfile(tamoxifen,facetBy=c("Tissue","Condition"))
# Barplots of the relative number of reads that overlap peaks vs. background reads, 
# as well and the proportion of reads overlapping blacklisted regions
plotRap(tamoxifen,facetBy=c("Tissue","Condition"))
plotFribl(tamoxifen,facetBy=c("Tissue","Condition"))
# sample clustering
plotCorHeatmap(tamoxifen,attributes=c("Tissue","Factor","Condition","Replicate"))
# pca plot
plotPrincomp(tamoxifen,attributes=c("Tissue","Condition"))

# facetBy/colourBy/lineBy param can change the sample grouped togother
# facetBy accepts a character vector corresponding to metadata column names 
# and will group plots according to those supplied.
# they can be provided by addMetaData
 extraMetadata <- data.frame(Sample = rownames(QCmetrics(exampleExp)),
 							FRiBL = QCmetrics(exampleExp)[,("RiBL%")],
							SSD =QCmetrics(exampleExp)[,("SSD")]
							)
plotCC(exampleExp,facetBy="Sample",lineBy="Tissue",addMetaData=extraMetadata,colourBy="SSD")


######################### single sample analysis #########################
 CTCF1 = ChIPQCsample("reads/SRR568129.bam", peaks="peaks/SRR568129_chr22_peaks.bed")
 # CTCF1 = QCsample(exampleExp,1)	# 从前面的 exampleExp 里面取一个 sample 来做演示
 # CTCF1 就是一个包含了 peak 在基因结构上的 counts ratio 对象 和一个 Granges 对象（bed文件来的）
 QCmetrics(CTCF1)
 plotCC(CTCF1)

 ##### QCsample 可以把已有的 ChIPQCsample 对象，同新的待分析的样本结合起来
 sampleList = QCsample(tamoxifen)
 tamoxifen = ChIPQC(sampleSheet, samples=sampleList)