## step 1: get segregating sites
module load java/1.8.0u31 TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools
# for each chromosome  $chrom
rm /scratch/devel/mkuhlwilm/arch/private/segsiteT_"$chrom".vcf.gz
module load TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools
bcftools view -m 2 -O z -V indels /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz > /scratch/devel/mkuhlwilm/arch/private/segsiteT_"$chrom".vcf.gz
tabix -f -p vcf /scratch/devel/mkuhlwilm/arch/private/segsiteT_"$chrom".vcf.gz

## step 2: filter for repeats
module load java/1.8.0u31 TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools
intersectBed -wa -v -header -a /scratch/devel/mkuhlwilm/arch/private/segsiteT_"$chrom".vcf.gz -b /project/devel/mkuhlwilm/RM.bed | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep "#|mapability=1" | bgzip > /scratch/devel/mkuhlwilm/arch/private/segsite_filt2_"$chrom".vcf.gz
tabix -f -p vcf /scratch/devel/mkuhlwilm/arch/private/segsite_filt2_"$chrom".vcf.gz

## step 3: identify imbalanced heterozygous sites
module load java/1.8.0u31 tabix xz gcc/6.3.0 bcftools/1.9 SAMTOOLS/1.0 BEDTools/2.26.0 R/3.2.0 tabix
ID=$SLURM_ARRAY_TASK_ID
fil=$(sed -n ${ID}p /scratch/devel/mkuhlwilm/arch/species.lst)
spec="gorilla" 
R CMD BATCH --vanilla --slave "--args ${spec}" /home/devel/mkuhlwilm/programs/filter_hets.R

## step 4: filtering for S*: only sites where all individuals have high quality genotypes
# list of gorilla individuals
fil="/scratch/devel/mkuhlwilm/arch/gorilla_noT.lst"
spec="gorilla"
/apps/BCFTOOLS/1.9/bin/bcftools view -S $fil -h /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz > /scratch/devel/mkuhlwilm/arch/subsets/gehdr.txt
# for each chromosome $chrom
zcat /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz | sed -e 's/;MLEAC=.*;MLEAF=.*;/;/g' | /apps/BCFTOOLS/1.9/bin/bcftools view -S $fil -a -U -V indels - | /apps/BCFTOOLS/1.9/bin/bcftools view -m 2 - | /apps/BCFTOOLS/1.9/bin/bcftools filter --exclude 'GT~"\."' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '(MIN(FMT/DP) >= 6) && (MAX(FMT/DP) <= 100) && (MIN(FMT/MQ) > 20)' | /apps/BCFTOOLS/1.9/bin/bcftools filter --include '((FMT/MQ0)/(FMT/DP)) < 0.1' | intersectBed -wa -v -header -a stdin -b /scratch/project/devel/mkuhlwilm/RM.bed | intersectBed -wa -v -header -a - -b <(zcat /scratch/devel/mkuhlwilm/arch/filter/Na_"$spec"_"$chrom".bed.gz ) | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep "#|mapability=1" | /apps/BCFTOOLS/1.9/bin/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT[\t%GT]\n' - | sed 's/\//|/g' | cat /scratch/devel/mkuhlwilm/arch/subsets/gehdr.txt - | bgzip -f > /scratch/devel/mkuhlwilm/arch/subsets/3_"$spec"_"$chrom".vcf.gz

## step 5: filtering for hmmix/Skov method
## get weights files using an R script
module load java/1.8.0u31 TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools
R CMD BATCH --vanilla --slave "--args ${chrom}" /home/devel/mkuhlwilm/programs/get_clbl.R

## get mutation rate files 
module load java/1.8.0u31 TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools
/apps/BCFTOOLS/1.9/bin/bcftools view -S /scratch/devel/mkuhlwilm/arch/"$spec".lst /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chr".vcf.gz -G | intersectBed -wa -v -header -a stdin -b /project/devel/mkuhlwilm/RM.bed | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep "#|mapability=1" | bedtools merge -i stdin > /scratch/devel/mkuhlwilm/arch/mura/chr_"$chr"mura_"$spec".bed
/apps/BCFTOOLS/1.9/bin/bcftools view -m2 -M2 -S /scratch/devel/mkuhlwilm/arch/"$spec".lst /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chr".vcf.gz | intersectBed -wa -header -a stdin -b /scratch/devel/mkuhlwilm/arch/mura/chr_"$chr"mura_"$spec".bed | vcftools --vcf - --counts --stdout --remove-indels > /dev/shm/mydata/"$chr"_"$spec".freq
awk -v c="$chr" '$1==c' /scratch/devel/mkuhlwilm/arch/N2_"$spec"_weights_float.txt > /dev/shm/mydata/"$chr"_"$spec"_weight.txt
python /home/devel/mkuhlwilm/programs/skov/Introgression-detection-master4/Estimate_mutationrate.py /dev/shm/mydata/"$chr"_"$spec".freq 100000 1000 /dev/shm/mydata/"$chr"_"$spec"_weight.txt /scratch/devel/mkuhlwilm/arch/mura/chr_"$chr"mura_"$spec".mut

spec=$(echo $spec | cut -f 6 -d "/" | cut -f 1 -d ".")
sed -e 's/0\.0/0\.1/g' /scratch/devel/mkuhlwilm/arch/mura/chr_1mura_"$spec".mut > /scratch/devel/mkuhlwilm/arch/mura/mura_"$spec".mut
for chrom in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
  do echo $chrom
  sed -e 's/0\.0/0\.1/g' /scratch/devel/mkuhlwilm/arch/mura/chr_"$chrom"mura_"$spec".mut >> /scratch/devel/mkuhlwilm/arch/mura/mura_"$spec".mut
  done

## get observations file using an R script
module load TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools
# per chromosome: filter data individual-wise
R CMD BATCH --vanilla --slave "--args ${chrom}" /home/devel/mkuhlwilm/programs/getgeno1.R

# per chromosome: turn into observation file format
R CMD BATCH --vanilla --slave "--args ${chrom}" /home/devel/mkuhlwilm/programs/getgeno2.R

# per group: check output and write files
R CMD BATCH --vanilla --slave "--args ${chrom}" /home/devel/mkuhlwilm/programs/getgeno3.R

## step 6: convert to eigenstrat format
R CMD BATCH --vanilla --slave "--args ${chrom}" /home/devel/mkuhlwilm/programs/geno_eigen.R
