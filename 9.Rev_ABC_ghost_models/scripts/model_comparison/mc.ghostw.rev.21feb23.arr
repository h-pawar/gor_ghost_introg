#!/bin/bash
# @ job_name = mc_ghostw
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/ghostw_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/ghostw_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 4:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-200

# Tue 21 Feb 2023 09:51:06 CET
# generate model comparison simulations for revised ghost -> w_anc model - sampling all parameters from priors

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/mc.ghostw.rev.24feb23.R  /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/mc_ghostwrev_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'mc.ghostw.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------

 
