#!/bin/bash
# @ job_name = ghost_abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/ghost/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/ghost/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-700

# infer  ghoste model parameter inference 

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata


#-----------------------------------------------------------------------------------------------------------------------

# Tue  5 Apr 2022 15:34:56 CEST - amended with final null posteriors from paramteter inference

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.1mar22.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/ghostabcsimul_1mar22_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'ghost.abc.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------

