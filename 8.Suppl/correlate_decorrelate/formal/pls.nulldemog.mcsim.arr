#!/bin/bash
# @ job_name = mc_nullpls
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 2:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 21-200

# 1-200

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#Mon 21 Feb 2022 12:43:09 CET
# model comparison simulations for ABC PLS 'reduced dimensionality ABC'
R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/decorrelatestats/pls.nulldemog.mcsim.R  /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/test.mc_nullpls_simul_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'mc.ghoste_sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------

