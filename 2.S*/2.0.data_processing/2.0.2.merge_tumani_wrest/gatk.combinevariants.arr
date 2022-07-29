#!/bin/bash
# @ job_name = combinevariants
# @ initialdir = /scratch/devel/hpawar/admix/sstar/vcf_24jun22
# @ output = /scratch/devel/hpawar/admix/sstar/vcf_24jun22/log/combinevariants_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/vcf_24jun22/log/combinevariants_%j_%a.err
# @ cpus_per_task = 4
# @ wall_clock_limit = 20:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID


#mkdir -p /scratch/devel/hpawar/admix/sstar/vcf_24jun22

#-----------------------------------------------------------------------------------------------------------------------
# Tue 21 Jun 2022 16:37:46 CEST
# Marina (MAE)
#This is the path to Tumanis folder where there are the BAMs and the gVCFs (by chromosome) that I have generated.
#/scratch/devel/malvarest/WGS_Gorillas/Tumani

#The VCF folder is the: 05_UnifiedGenotyper

#Hope this works for you. If there is something missing please let me know.
#-----------------------------------------------------------------------------------------------------------------------

# MK: I have not applied any filter to the original genotypes, but merged them with the other individuals. 
#So you may start by merging the "new" Tumani with chr9 of the "old sample set", and ideally the very same way:

#-----------------------------------------------------------------------------------------------------------------------

# 1: merge
module load java GATK/3.6 tabix

java -Xmx4g -Djava.io.tmpdir=/scratch_tmp -jar /apps/GATK/3.6/GenomeAnalysisTK.jar  \
-T CombineVariants -R /home/devel/marcmont/scratch/snpCalling_hg19/chimp/assembly/BWA/hg19.fa \
--variant /scratch/devel/mkuhlwilm/gvcfs/greatapeN_9.vcf.gz \
--variant /scratch/devel/malvarest/WGS_Gorillas/Tumani/05_UnifiedGenotyper/Gbg-Tumani_chr9.g.vcf.gz | bgzip -f > /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz

tabix /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz

