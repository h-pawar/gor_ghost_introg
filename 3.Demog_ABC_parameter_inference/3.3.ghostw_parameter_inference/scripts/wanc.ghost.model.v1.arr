#!/bin/bash
# @ job_name = wghost_abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/ghost/w%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/ghost/w%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 429-543, 639-700

# array = 2-700 # send first for a test job

# generate ghostw ABC parameter inference simulations

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix

mkdir /dev/shm/mydata
cd /dev/shm/mydata

#-----------------------------------------------------------------------------------------------------------------------
#Fri  8 Apr 2022 11:46:28 CEST
# final weighted median posteriors from final null abc parameter inference
R CMD BATCH --vanilla --slave "--args ${ID}"  /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.1apr22.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/wghostabcsimul_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'wanc.ghost.abc.sim*' -execdir rm -- '{}' \;

exit 
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
