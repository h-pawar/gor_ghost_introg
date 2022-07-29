#!/bin/bash
# @ job_name = filt_imbalancedhets
# @ initialdir = /scratch/devel/hpawar/admix/sstar/vcf_24jun22
# @ output = /scratch/devel/hpawar/admix/sstar/vcf_24jun22/log/filt_imbhets_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/vcf_24jun22/log/filt_imbhets_%j_%a.err
# @ cpus_per_task = 4
# @ wall_clock_limit = 24:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID

# 2: filter imbalanced hets
# follows from merging the "new" Tumani with chr9 of the "old sample set" - which used the script - /scratch/devel/hpawar/admix/sstar/scripts/gatk.combinevariants.arr

module load java/1.8.0u31 xz gcc/6.3.0 bcftools/1.9 SAMTOOLS/1.0 BEDTools/2.26.0 R/3.2.0 tabix
R CMD BATCH --vanilla --slave  /scratch/devel/hpawar/admix/sstar/scripts/filt.hets.MK.R
exit

#-----------------------------------------------------------------------------------------------------------------------
# how MK was running - (but now removign the ty arg - b/c only have gor & chr 9) in filt_hets.arr
#module load java/1.8.0u31 tabix xz gcc/6.3.0 bcftools/1.9 SAMTOOLS/1.0 BEDTools/2.26.0 R/3.2.0 tabix
## filter allele imbalance on hets
#ID=$SLURM_ARRAY_TASK_ID
#fil=$(sed -n ${ID}p /scratch/devel/mkuhlwilm/arch/species.lst)
#spec=$(echo $fil | cut -f 6 -d "/" | cut -f 1 -d ".")
#echo $spec 
#for chrom in {1..22}; do echo $chrom; tabix -f /scratch/devel/mkuhlwilm/arch/subsets/1_"$spec"_"$chrom".vcf.gz;done
#tabix -f /scratch/devel/mkuhlwilm/arch/subsets/1_"$spec"_X.vcf.gz
#R CMD BATCH --vanilla --slave "--args ${spec}" /home/devel/mkuhlwilm/programs/filter_hets.R /home/devel/mkuhlwilm/logs/hefi_${spec}
#exit
#-----------------------------------------------------------------------------------------------------------------------

