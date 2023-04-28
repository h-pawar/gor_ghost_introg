#!/bin/bash
# @ job_name = abc_emp
# @ initialdir = /scratch/devel/hpawar/admix/abc/emp.data/
# @ output = /scratch/devel/hpawar/admix/abc/emp.data/log/stats_chr%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/emp.data/log/stats_chr%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 4:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-23

# array = 1-22

ID=$SLURM_ARRAY_TASK_ID

# calc summary statistics for the empirical data in the same way (in 40kb windows) as for the simulated data in test.abc.model.v4.R
#  heterozygosity, number of segregating sites, pairwise Fst, pi, Tajima's D)

#-----------------------------------------------------------------------------------------------------------------------

module load gcc/6.3.0 R/3.4.2 xz/5.2.2  BCFTOOLS/1.6
if [[ "$ID" = 23 ]]; then ID="X"; fi

#-----------------------------------------------------------------------------------------------------------------------
# Tue 24 Aug 2021 10:52:20 CEST
# need to recalculate summary stats for empirical data b/c order of individuals in the vcf != order of individuals simulated
  # & output sums of het, seg sites (instead of means)   
Rscript --vanilla /scratch/devel/hpawar/admix/abc/simul/scripts/stats.emp.data.clean3.R $ID
#-----------------------------------------------------------------------------------------------------------------------

