# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE47nnn/GSE47197/matrix/GSE47197_series_matrix.txt.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1009/matrix/GSE1009_series_matrix.txt.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1009/suppl/GSE1009_RAW.tar
# http://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&mode=csv&series=47197

getGSELink <- function (studyID='GSE1009', down=F, destdir='./') {
	supp_link = paste0('ftp://ftp.ncbi.nlm.nih.gov/geo/series/',
						substr(studyID,1,nchar(studyID)-3),
						'nnn/', studyID,
						'/suppl/', studyID,
						'_RAW.tar')

	meta_link = paste0('http://www.ncbi.nlm.nih.gov/geo/browse/?view=sample&mode=csv&series=',
						substr(studyID, 4, nchar(studyID)))

	matrix_link = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
						substr(studyID, 1, nchar(studyID) - 3),
						"nnn/", studyID, 
						"/matrix/", studyID,
						"_series_matrix.txt.gz")

	print(paste0('RAW  data:    ', supp_link))
	print(paste0('meta data:    ', meta_link))
	print(paste0('matrix   :    ', matrix_link))

	if (down) {
		download.file(supp_link, destfile=paste0(destdir, studyID, '_RAW.tar'))
		download.file(supp_link, destfile=paste0(destdir, studyID, '_meta.csv'))
		download.file(supp_link, destfile=paste0(destdir, studyID, '_matrix.txt.gz'))
	}
}