#!/bin/bash
# @ job_name = mc_ghoste
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/ghoste_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/ghoste_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 4:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 11-200

# scale up to 200 (target 10,000 reps total)

# Tue 12 Apr 2022 16:51:42 CEST
# generate model comparison simulations for ghoste model B

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/mc.ghoste.12apr22.R  /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/mc_ghoste_12apr22_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'mc.ghoste.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------

     
