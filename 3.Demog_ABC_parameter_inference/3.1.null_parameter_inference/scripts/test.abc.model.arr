#!/bin/bash
# @ job_name = abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 2-700

ID=$SLURM_ARRAY_TASK_ID


#-----------------------------------------------------------------------------------------------------------------------

# Generate ms simulations for abc + output relevant summary statistics
#-----------------------------------------------------------------------------------------------------------------------

module load gcc/6.3.0 R/3.4.2 tabix 

mkdir /dev/shm/mydata
cd /dev/shm/mydata

#-----------------------------------------------------------------------------------------------------------------------
# Sat 11 Sep 2021 11:35:35 CEST

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R  /scratch/devel/hpawar/admix/abc/simul/log/abcsimul_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;


exit
#-----------------------------------------------------------------------------------------------------------------------

