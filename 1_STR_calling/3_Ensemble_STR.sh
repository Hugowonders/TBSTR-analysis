#!/bin/bash

wd=$your_working_directory

reference=$wd/GRCh38_full_analysis_set_plus_decoy_hla.fa

### prepare the HipSTR correction file proviced in the GitHub repository

### correct inconsistent alleles in the dumped HipSTR VCF file
python3 $wd/Hipstr_correction.py ${wd}/dumped_gangstr_variants.vcf.gz \
	${wd}/dumped_gangstr_variants.corrected.vcf

bgzip -f ${wd}/dumped_gangstr_variants.corrected.vcf
bcftools index -tf ${wd}/dumped_gangstr_variants.corrected.vcf.gz

### consensus ensemble of HipSTR and GanSTR variants
EnsembleTR \
	--vcfs ${wd}/dumped_gangstr_variants.corrected.vcf.gz,${wd}/dumped_gangstr_variants.vcf.gz \
	--ref reference \
	--out ${wd}/raw_ensemble_variants.vcf

### Compress and index the output file
bgzip -f ${wd}/raw_ensemble_variants.vcf
bcftools index -tf ${wd}/raw_ensemble_variants.vcf.gz