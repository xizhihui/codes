# 如果不是Ensembl gtf, 是gencode gtf的话,要把gene_biotype改成gene_type
perl -alne '/gene_id "(.*?)";.*?gene_biotype "(.*?)";/; print $1."\t".$2;' ~/database/reference/hg38/Homo_sapiens.GRCh38.92.gtf | sort | uniq | grep -E 'non_coding|3prime_overlapping_ncRNA|antisense|lincRNA|retained_intron|sense_intronic|sense_overlapping|macro_lncRNA|bidirectional_lncRNA' > lncRNA.txt