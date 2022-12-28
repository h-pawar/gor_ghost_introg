#!/bin/bash
# @ job_name = vf_allelefreq
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/volcanofinder/log/allelefreq_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/volcanofinder/log/allelefreq_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 6:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-5

# script amended - Mon  7 Feb 2022 18:51:48 CET - to incorporate /scratch/devel/hpawar/admix/volcanofinder/scripts/test.filter.anc.alleles.arr
# steps to add additional filtering step - to rm sites where have multiple ancestral alleles
# most jobs finished in 2h, chr 1-5 did not

# try first for chr21, if all fine run for 1-20
#-----------------------------------------------------------------------------------------------------------------------

ID=$SLURM_ARRAY_TASK_ID

#Fri  4 Feb 2022 11:26:00 CET
# 1) generate the allele freq file (required input for volcanofinder)
    # - needs derived allele counts on biallelic sites 

module load BCFTOOLS/1.12
cd /scratch/devel/hpawar/admix/volcanofinder/input/2feb22
# filtered & annotated with ancestral allele (human ref vs macaque - now gor vs human ref)
vcfdir="/scratch/devel/mkuhlwilm/arch/subsets/full/4_gorilla_"
vcf=$vcfdir$ID".vcf.gz"

#-----------------------------------------------------------------------------------------------------------------------

# 0) subset to only easterns
subsvcfdir="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22"
subsvcf=$subsvcfdir"/4_"$ID"_Egor.vcf.gz"

## = commenting out for failed reps 1-5,b/c have already run these steps successfully

##bcftools view -s  Gorilla_beringei_beringei-Bwiruka,Gorilla_beringei_beringei-Imfura,Gorilla_beringei_beringei-Kaboko,Gorilla_beringei_beringei-Kahungye,Gorilla_beringei_beringei-Katungi,Gorilla_beringei_beringei-Maisha,Gorilla_beringei_beringei-Nyamunwa,Gorilla_beringei_beringei-Semehe,Gorilla_beringei_beringei-Tuck,Gorilla_beringei_beringei-Turimaso,Gorilla_beringei_beringei-Umurimo,Gorilla_beringei_beringei-Zirikana,Gorilla_beringei_graueri-9732_Mkubwa,Gorilla_beringei_graueri-A929_Kaisi,Gorilla_beringei_graueri-A967_Victoria,Gorilla_beringei_graueri-Dunia,Gorilla_beringei_graueri-Itebero,Gorilla_beringei_graueri-Mukokya,Gorilla_beringei_graueri-Ntabwoba,Gorilla_beringei_graueri-Pinga,Gorilla_beringei_graueri-Tumani $vcf -Oz -o $subsvcf

#-----------------------------------------------------------------------------------------------------------------------

# 1) filter by biallelic sites
##module unload gcc/6.3.0
##module load gcc/latest zlib/1.2.8 VCFTOOLS/0.1.15 TABIX/0.2.6

out="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_"$ID"_Egor.biallelic"

##vcftools --gzvcf ${subsvcf}  --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out $out

pref="4_"$ID"_Egor.biallelic"

##bgzip -c "$pref".recode.vcf > "$pref".vcf.gz # compress
##tabix -p vcf "$pref".vcf.gz # create index file
##rm "$pref".recode.vcf # remove intermediate files

#-----------------------------------------------------------------------------------------------------------------------

# 1.1) filter to remove loci with multiple ancestral alleles 
    # b/c don't know the ancestral state of these loci (insertions/deletions)

##module unload gcc/latest zlib/1.2.8 VCFTOOLS/0.1.15 TABIX/0.2.6
##module load BCFTOOLS/1.12


# extract rows where AA=x (rather than AA=GAT eg)
bcftools view -H 4_"$ID"_Egor.biallelic.vcf.gz |  cut -f1-2,8 | awk -F ";" '{print $1}' | awk ' length($3) == 4 ' |  awk -v OFS='\t' '{print $1,$2}' > 4_"$ID"_Egor.biallelic.keep.regions.txt


# keep sites where there is one ancestral allele given
bcftools view -R 4_"$ID"_Egor.biallelic.keep.regions.txt 4_"$ID"_Egor.biallelic.vcf.gz -Oz -o 4_"$ID"_Egor.biallelic.filt.vcf.gz


#-----------------------------------------------------------------------------------------------------------------------


# 2) generate allele freq file

filtvcf="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_"$ID"_Egor.biallelic.filt.vcf.gz"

module unload gcc/6.3.0
module load gcc/latest zlib/1.2.8 VCFTOOLS/0.1.15 TABIX/0.2.6

vcftools --counts2 --derived --gzvcf $filtvcf --stdout | awk 'NR<=1 {next} {print $2"\t"$6"\t"$4"\t0"}' > 4_"$ID"_SF2.input
# add header line
echo -e "position\tx\tn\tfolded" | cat - 4_"$ID"_SF2.input > tmp.txt  && mv tmp.txt 4_"$ID"_SF2.input

#-----------------------------------------------------------------------------------------------------------------------

# filter to remove 0 sites (ref homozygotes)
awk '{ if ($2 != 0) { print } }' /scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_"$ID"_SF2.input > /scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_"$ID"_SF2_polym.input
