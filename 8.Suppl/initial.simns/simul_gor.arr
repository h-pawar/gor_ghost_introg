#!/bin/bash
# @ job_name = simul_gor
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/sstar/log/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 2
# @ wall_clock_limit = 4:00:00
# @ class = normal
# @ requeue = 1
# @ array = 134,25,113,49

ID=$SLURM_ARRAY_TASK_ID


module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix 
source venv2/bin/activate

#Rscript --vanilla /scratch/devel/hpawar/admix/sstar/scripts/simul_gor.R $ID

Rscript --vanilla /scratch/devel/hpawar/admix/sstar/scripts/troubleshoot.simul.R $ID

