#!/bin/bash
# @ job_name = vf_gtmatrix
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/volcanofinder/log/gtmatrix_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/volcanofinder/log/gtmatrix_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 1:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-5


# Tue  8 Feb 2022 11:54:18 CET
# calc sfs of the autosomes
# first step generate GT matrix from vcf

#-----------------------------------------------------------------------------------------------------------------------
# 2) generate sfs
#- need unfolded sfs w ancestral allele call
#- ie generate sfs from polarised vcf

# calc sfs -> normalise so all sites sum to 1 -> then drop the first line (the 0 category)
#-----------------------------------------------------------------------------------------------------------------------

# 1. vcf -> genotype matrix ( genotypes coded as 0, 1, 2, -1 (missing data))

module load BCFTOOLS/1.12
cd /scratch/devel/hpawar/admix/volcanofinder/input/2feb22

ID=$SLURM_ARRAY_TASK_ID

# tag for VCF file (vcf file with format "vcffile".vcf)
vcffile="4_"$ID"_Egor.biallelic.filt";      
# Get GT field with bcftools
bcftools query -f '[%GT\t]\n' ${vcffile}.vcf.gz > ${vcffile}.GT
# Replace 0/0 by 0 # Replace 0/1 by 1 # Replace 1/1 by 2 # Replace ./. by -1 # Replace . by -1
#sed -i "s/0\/0/0/g;s/0\/1/1/g;s/1\/0/1/g;s/1\/1/2/g;s/\.\/./-1/g;s/\./-1/g" ${vcffile}.GT - not working b/c | instead of /

sed -i "s/0|0/0/g;s/0|1/1/g;s/1|0/1/g;s/1|1/2/g;s/\.|./-1/g;s/\./-1/g" ${vcffile}.GT

#-----------------------------------------------------------------------------------------------------------------------
# 2. then combine .GT files per chr -> genome-wide GT
# 3. calc sfs in R
