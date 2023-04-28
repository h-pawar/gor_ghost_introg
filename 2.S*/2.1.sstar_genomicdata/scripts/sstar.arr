#!/bin/bash
# @ job_name = sstar
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/sstar/log/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 2
# @ wall_clock_limit = 24:00:00
# @ class = normal
# @ requeue = 1
# @ array = 23

## this is running Sstar on the genome

# load venv - where have installed necessary python modules for sstar - bitarry, numpy, pandas, numexpr and tables

module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 
source venv2/bin/activate

ID=$SLURM_ARRAY_TASK_ID

chrom=$ID

#-----------------------------------------------------------------------------------------------------------------------
## the program needs some parameters
# -vcfz = gzipped VCF file to analyze
# -indf = sample_pop_mapping_file.txt : cols: sample	pop	super_pop	gender
# -ref-pops -target-pops = populations to use as in- and outgroup; could be GG, GB (among that also GBB, GBG)
# -ancbsg = ancestralized genome file in weird format
# -archaic-vcf = "archaic" genome; I just feed in an ancestralized vcf file, but ignore later output based on that
# -r callable regions in weird format

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# GBB = EM (Gorilla beringei beringei, eastern clade - mountain gorilla)
# GBG = EL (Gorilla beringei graueri, eastern clade - eastern lowland)
# GGD = WC (Gorilla gorilla diehli, western clade - cross river)
# GGG = WL (Gorilla gorilla gorilla, western clade - western lowland)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# B) Calc S* for the subspecies comparisons on the autosomes - run for array 1-22

# 1) WL outgroup, EL ingroup (GGG outg, GBG ing)
#echo "GBG"
#python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_"$chrom".vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GGG -target-pops GBG -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBG_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 2) WL outgroup, EM ingroup (GGG outg, GBB ing)
#echo "GBB"
#python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_"$chrom".vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GGG -target-pops GBB -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBB_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg


# 3) EL outgroup, WL ingroup (GBG outg, GGG ing)
#echo "GGG"
#python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_"$chrom".vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GBG -target-pops GGG -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GGG_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 4) EL outgroup, WC ingroup (GBG outg, GGD ing)
#echo "GGD"
#python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_"$chrom".vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GBG -target-pops GGD -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GGD_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_"$chrom"b.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

#exit

#-----------------------------------------------------------------------------------------------------------------------

# Calc S* for the subspecies comparisons on the X chr (chr 23) -  run this section

# 1) WL outgroup, EL ingroup (GGG outg, GBG ing)
echo "GBG"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_X.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GGG -target-pops GBG -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBG_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_Xb.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 2) WL outgroup, EM ingroup (GGG outg, GBB ing)
echo "GBB"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_X.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GGG -target-pops GBB -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBB_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_Xb.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg


# 3) EL outgroup, WL ingroup (GBG outg, GGG ing)
echo "GGG"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_X.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GBG -target-pops GGG -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GGG_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_Xb.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

# 4) EL outgroup, WC ingroup (GBG outg, GGD ing)
echo "GGD"
python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py -s-star -vcfz /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_X.vcf.gz -indf /scratch/devel/hpawar/admix/sstar/scripts/gorilla.subsp.slst -winlen 40000 -winstep 30000 -winchrom $chrom -ref-pops GBG -target-pops GGD -o /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GGD_"$chrom".star -ancbsg /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.bsg -archaic-vcf /scratch/devel/mkuhlwilm/pseudoarc/ga_Xb.vcf.gz -r /scratch/devel/mkuhlwilm/pseudoarc/gor_callable.bed.bbg

exit

#-----------------------------------------------------------------------------------------------------------------------



