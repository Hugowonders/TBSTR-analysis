#!/bin/bash

wd=$your_working_directory

bcftools filter ${wd}/raw_ensemble_variants.vcf.gz \
	-e 'F_PASS(GT=\"mis\")>0.5' | \
	bcftools filter \
	-i 'FORMAT/SCORE>=0.9' \
	-S '.' \
	--write-index=tbi \
	-Oz \
	-o ${wd}/filtered_ensemble_variants.vcf.gz