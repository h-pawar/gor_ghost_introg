#!/bin/bash
# @ job_name = sstar
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/sstar/log/chr9_newT_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/chr9_newT_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 2
# @ wall_clock_limit = 24:00:00
# @ class = normal
# @ requeue = 1
# @ array = 9

#-----------------------------------------------------------------------------------------------------------------------
# Fri  8 Jul 2022 10:17:00 CEST
## Run S* on newly generated chr 9 (after GBG individual 9 Tumani has been reprocessed)
# load venv - where have installed necessary python modules for sstar - bitarry, numpy, pandas, numexpr and tables
module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17
source venv2/bin/activate

ID=$SLURM_ARRAY_TASK_ID
# only run this script for array = 9 (chr 9)

chrom=$ID

#-----------------------------------------------------------------------------------------------------------------------
# GBB = EM (Gorilla beringei beringei, eastern clade - mountain gorilla)
# GBG = EL (Gorilla beringei graueri, eastern clade - eastern lowland)
# GGD = WC (Gorilla gorilla diehli, western clade - cross river)
# GGG = WL (Gorilla gorilla gorilla, western clade - western lowland)

# Calc S* for the subspecies comparisons on newly generated chr 9 (with newly processed gbg-tumani in the input data)

#-----------------------------------------------------------------------------------------------------------------------

# 1) WL outgroup, EL ingroup (GGG outg, GBG ing)
echo "GBG"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz  /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/vcf_24jun22/gorilla.subsp.newT.slst  -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GGG -target-pops GBG -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBG_newT_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 2) WL outgroup, EM ingroup (GGG outg, GBB ing)
echo "GBB"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz  /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/vcf_24jun22/gorilla.subsp.newT.slst  -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GGG -target-pops GBB -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBB_newT_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 3) EL outgroup, WL ingroup (GBG outg, GGG ing)
echo "GGG"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz  /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/vcf_24jun22/gorilla.subsp.newT.slst  -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GBG -target-pops GGG -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GGG_newT_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 4) EL outgroup, WC ingroup (GBG outg, GGD ing)
echo "GGD"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz  /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/vcf_24jun22/gorilla.subsp.newT.slst  -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GBG -target-pops GGD -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GGD_newT_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg
#-----------------------------------------------------------------------------------------------------------------------

exit