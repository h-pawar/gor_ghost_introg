#!/bin/bash
# @ job_name = mc_nullplusghost
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/null_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/null_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 2:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 8-200

# Tue 12 Apr 2022 16:51:42 CEST
# generate model comparison simulations for null pl ghost
# using weighted median posteriors from final abc parameter inference for null model
# scale up to 200 (target 10,000 reps total)

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------

# 250 windows per iter - (per rep = 50 iter)
# model comparison simulations

# null demog model + non-interacting ghost
# filter out sites fixed in all gorilla individuals before calculating the summary stats (same filter prev applied to empirical data)

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/mc.nullplghost.12apr22.R  /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/mc_nullplusghost_12apr22_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'mc.null.plusghost.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------

     
