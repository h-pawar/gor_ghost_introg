#!/bin/bash
# @ job_name = het_emp
# @ initialdir = /scratch/devel/hpawar/admix/abc/emp.data/
# @ output = /scratch/devel/hpawar/admix/abc/emp.data/log/seg_chr%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/emp.data/log/seg_chr%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 4:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-22

# array = 1-22

ID=$SLURM_ARRAY_TASK_ID

# calc summary statistics for the empirical data in the same way (in 40kb windows) as for the simulated data in test.abc.model.v4.R
#  heterozygosity, number of segregating sites, pairwise Fst, pi, Tajima's D)

#-----------------------------------------------------------------------------------------------------------------------

module load gcc/6.3.0 R/3.4.2 xz/5.2.2  BCFTOOLS/1.6

#------------------------------------------------------------------------------------------------------------------------
#Tue 31 Aug 2021 12:36:40 CEST - change way of calculating seg sites
# first run
#Rscript --vanilla /scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.31aug.R $ID
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
#Mon 27 Sep 2021 16:21:04 CEST
# calc for empirical data:
                # 2) number of population-wise fixed sites and the number of population-wise segregating sites; 
                # 3) fixed sites per individual
Rscript --vanilla /scratch/devel/hpawar/admix/abc/simul/scripts/segs.emp.27sep.R $ID
#------------------------------------------------------------------------------------------------------------------------

 
