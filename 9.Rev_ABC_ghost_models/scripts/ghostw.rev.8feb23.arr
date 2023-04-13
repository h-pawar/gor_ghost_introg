#!/bin/bash
# @ job_name = ghostw_abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

# array = 2-700 # send first for a test job

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------
# Wed  8 Feb 2023 10:59:23 CET
# re-perform parameter inference for ghost models, sampling all parameters from priors

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghostw.rev.8feb23.R /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/ghostwrev_simul_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'ghostw.rev.abc.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------
