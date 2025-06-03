#!/bin/bash

## Call level filtering after EnsembleTR
wd=$your_working_directory

bcftools filter ${wd}/raw_ensemble_variants.vcf.gz \

### Remove loci with missing rate greater than 0.5
	-e 'F_PASS(GT=\"mis\")>0.5' | \
	bcftools filter \

### Remove calls with concensus genotype score less than 0.9
	-i 'FORMAT/SCORE>=0.9' \

### Set missing genotypes as "."
	-S '.' \

### Compress and index the output file
	--write-index=tbi \
	-Oz \
	-o ${wd}/filtered_ensemble_variants.vcf.gz