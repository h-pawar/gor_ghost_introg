#!/bin/bash
# @ job_name = proc_tumani
# @ initialdir = /scratch/devel/hpawar/admix/sstar/vcf_24jun22
# @ output = /scratch/devel/hpawar/admix/sstar/vcf_24jun22/log/proc_t_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/vcf_24jun22/log/proc_t_%j_%a.err
# @ cpus_per_task = 4
# @ wall_clock_limit = 8:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 9

ID=$SLURM_ARRAY_TASK_ID

#-----------------------------------------------------------------------------------------------------------------------
# follows from 
#/scratch/devel/hpawar/admix/sstar/scripts/gatk.combinevariants.arr # merge newly processed with rest of the samples
# /scratch/devel/hpawar/admix/sstar/scripts/filt.hets.MK.R, called by /scratch/devel/hpawar/admix/sstar/scripts/filt.hets.MK.arr # filter imbalanced hets
#-----------------------------------------------------------------------------------------------------------------------

module load java/1.8.0u31 tabix xz gcc/6.3.0 bcftools/1.9 SAMTOOLS/1.0 BEDTools/2.26.0 R/3.2.0

fil="/scratch/devel/hpawar/admix/sstar/scripts/gorilla_noT.lst" # modified with the proper name of the new Tumani

# amend header - to only gorillas incl new tumani
#bcftools view -S /scratch/devel/hpawar/admix/sstar/scripts/gorilla_noT.lst -h /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz > /scratch/devel/hpawar/admix/sstar/vcf_24jun22/hdr.txt 

zcat /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz | sed -e 's/;MLEAC=.*;MLEAF=.*;/;/g' | /apps/BCFTOOLS/1.9/bin/bcftools view -S $fil -a -U -V indels - | /apps/BCFTOOLS/1.9/bin/bcftools view -m 2 - | /apps/BCFTOOLS/1.9/bin/bcftools filter --exclude 'GT~"\."' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '(MIN(FMT/DP) >= 6) && (MAX(FMT/DP) <= 100) && (MIN(FMT/MQ) > 20)'  | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '((FMT/MQ0)/(FMT/DP)) < 0.1' |  intersectBed -wa -v -header -a stdin -b /scratch/project/devel/mkuhlwilm/RM.bed | intersectBed -wa -v -header -a - -b <(cat /scratch/devel/hpawar/admix/sstar/vcf_24jun22/Na_gor_chr9.bed )  | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep "#|mapability=1" | /apps/BCFTOOLS/1.9/bin/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n' - | sed 's/\//|/g' | cat /scratch/devel/hpawar/admix/sstar/vcf_24jun22/hdr.txt - | bgzip -f > /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz 



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# this step is unnecessary
# Remove the old tumani (so there is only the newly processed version in the input for s*) & subsets to only gorilla vcf
#bcftools view -s^ "Gorilla_beringei_graueri-Tumani" /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz -Oz -o /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.1.vcf.gz
#bcftools view -S /scratch/devel/hpawar/admix/sstar/scripts/gor.samples.chr9.30jun22 /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz -Oz -o /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.1.vcf.gz
#tabix /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.1.vcf.gz
# this step has run

# this step gave errors
# MK: I think the main issue comes from not trimming alternative alleles - 
#there are SNPs where you have an alternative allele in the ALT column, 
#but no individual with any read, that can cause problems. Thats why I use 


# 3: finally, adapt and modify the filtering from the follwing:
#fil="/scratch/devel/hpawar/admix/sstar/scripts/gorilla_noT.lst" # modified with the proper name of the new Tumani
#spec=$(echo $fil | cut -f 6 -d "/" | cut -f 1 -d ".")
#echo $spec
#zcat /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.1.vcf.gz | sed -e 's/;MLEAC=.*;MLEAF=.*;/;/g' | /apps/BCFTOOLS/1.9/bin/bcftools view -S $fil -a -U -V indels - | /apps/BCFTOOLS/1.9/bin/bcftools view -m 2 - | /apps/BCFTOOLS/1.9/bin/bcftools filter --exclude 'GT~"\."' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '(MIN(FMT/DP) >= 6) && (MAX(FMT/DP) <= 100) && (MIN(FMT/MQ) > 20)' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '((FMT/MQ0)/(FMT/DP)) < 0.1' | intersectBed -wa -v -header -a stdin -b /scratch/project/devel/mkuhlwilm/RM.bed | intersectBed -wa -v -header -a - -b <(cat /scratch/devel/hpawar/admix/sstar/vcf_24jun22/Na_gor_chr9.bed ) | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep "#|mapability=1" | /apps/BCFTOOLS/1.9/bin/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n' - | sed 's/\//|/g' | cat /scratch/devel/mkuhlwilm/arch/subsets/gehdr.txt - | bgzip -f > /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz 

#This filters empty genotypes, low/high coverage, repeatmask, imbalanced hets, and mapability, to spit out the input for S*.
