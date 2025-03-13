#!/bin/bash

wd=$your_working_directory

### remove sites within segmental duplication regions in GRCh38
segdup=${wd}/hg38.segdups.bed.gz

### prepare the lists of VCF files outputted by HipSTR and GangSTR
rawHipstrVcfList=${wd}/raw_hipstr_vcf_list.txt

rawGangstrVcfList=${wd}/raw_gangstr_vcf_list.txt

### merge individual output VCF files of HipSTR
mergeSTR --vcfs-list ${rawHipstrVcfList} --vcftype hipstr \
--out ${wd}/merged_hipstr_variants

bgzip -f ${wd}/merged_hipstr_variants.vcf

bcftools index -tf ${wd}/merged_hipstr_variants.vcf.gz

### merge individual output VCF files of GangSTR
mergeSTR --vcfs-list ${rawGangstrVcfList} --vcftype gangstr \
	--out ${wd}/merged_hipstr_variants

bgzip -f ${wd}/merged_hipstr_variants.vcf
bcftools index -tf ${wd}/merged_hipstr_variants.vcf.gz

### dump low-quality sites and genotypes of merged HipSTR variants
dumpSTR --vcf ${wd}/merged_hipstr_variants.vcf.gz --vcftype hipstr \
	--zip --drop-filtered \
	--min-locus-callrate 0.5 \
	--min-locus-hwep 1e-20 \
	--filter-regions $segdup \
	--out ${wd}/dumped_hipstr_variants

bgzip ${wd}/dumped_hipstr_variants.vcf.gz

### dump low-quality sites and genotypes of merged GangSTR variants

dumpSTR --vcf ${wd}/merged_gangstr_variants.vcf.gz --vcftype gangstr \
	--zip --drop-filtered \
	--min-locus-callrate 0.5 \
	--min-locus-hwep 1e-20 \
	--filter-regions $segdup \
	--gangstr-filter-spanbound-only \
	--gangstr-filter-badCI \
	--out ${wd}/dumped_gangstr_variants

bgzip ${wd}/dumped_gangstr_variants.vcf.gz