#!/bin/bash
# Annotating variant consequences of genome-wide STRs
# Variant Effect Predictor v111
wd=$your_working_directory

# Input VCF file
vcf=$wd/your_vcf_file.vcf.gz

# Chache directory for VEP
vepCacheDir=$wd/vep_chache_dir

# Reference genome
fasta=$wd/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.gz

# Annotation files for VEP plugins
# Refer to https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
utrFile=$wd/uORF_5UTR_GRCh38_PUBLIC.txt
spliceVault=$wd/SpliceVault_data_GRCh38.tsv.gz
mutfuncDB=$wd/mutfunc_data.db
maveDB=$wd/MaveDB_variants.tsv.gz

# Run VEP in offline mode
vep -i $vcf --fork 3 -o $wd/vep_out.vcf --assembly GRCh38 \
--cache --cache_version 111 --dir $vepCacheDir \
--fasta $fasta --force_overwrite --offline \
--hgvs --symbol --canonical --numbers --allele_number \
--biotype --protein --domains --uniprot --tsl \
--plugin SpliceVault,file=$spliceVault \
--plugin TSSDistance \
--plugin UTRAnnotator,file=$utrFile \
--plugin mutfunc,db=$mutfuncDB \
--plugin MaveDB,file=$maveDB,single_aminoacid_changes=0 \
--plugin Downstream

# Compress the output file
gzip -f $wd/vep_out.vcf