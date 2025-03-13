#！/bin/bash
wd=$your_working_directory

### prepare a list of ~100 bam files for each batch
bamFileList=$wd/bam_file_list.txt

### all analysis were based on GRCh38 reference genome
reference=$wd/GRCh38_full_analysis_set_plus_decoy_hla.fa

### curated STR catalogs were obtained from the GitHub repository of HipSTR and GangSTR
### STR catalogs were merged and filtered to keep non-overlapping site with motif ranging from 2-6 bp and referene allele < 150 bp

hipstrCatalog=$wd/hipstr_catalog.bed

gangstrCatalog=$wd/gangstr_catalog.bed


# HipSTR genotyping
HipSTR \
	--fasta $reference \
	--regions $hipstrCatalog \
	--bams $bams \
	--max-reads 2000000 \
	--max-str-len 500 \
	--silent \
	--str-vcf $wd/hipstr_output.vcf.gz

bcftools index -tf $wd/hipstr_output.vcf.gz

# GangSTR genotyping
GangSTR \
	--ref $reference \
	--regions $gangstrCatalog \
	--bam $bams \
	--max-proc-read 2000000 \
	--quiet \
	--out $wc/gangstr_output

bgzip $wc/gangstr_output.vcf

bcftools index -tf $wc/gangstr_output.vcf.gz