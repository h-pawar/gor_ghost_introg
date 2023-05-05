#!/bin/bash
# @ job_name = crossval
# @ initialdir = /scratch/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/crossval_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/crossval_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 4:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

#Wed  1 Mar 2023 09:41:50 CET

#ID=$SLURM_ARRAY_TASK_ID

module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------

R CMD BATCH --vanilla /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/crossvalidation.1mar23.R  /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/crossvalidation.log

#-----------------------------------------------------------------------------------------------------------------------

 
