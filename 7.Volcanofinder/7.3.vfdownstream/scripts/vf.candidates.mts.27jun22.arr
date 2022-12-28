#!/bin/bash
# @ job_name = vf_mts
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/volcanofinder/log/mts_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/volcanofinder/log/mts_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 6:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID

module load R/4.0.1
#Rscript --vanilla  /scratch/devel/hpawar/admix/volcanofinder/scripts/vf.candidates.mts.27jun22.R

#-----------------------------------------------------------------------------------------------------------------------
#Thu 21 Jul 2022 16:51:47 CEST
Rscript --vanilla  /scratch/devel/hpawar/admix/volcanofinder/scripts/vf.overlap.s*.skov.mts.21jul22.R
#-----------------------------------------------------------------------------------------------------------------------

