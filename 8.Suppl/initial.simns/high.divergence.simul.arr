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
# @ array = 2-158

ID=$SLURM_ARRAY_TASK_ID


module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix 
source venv2/bin/activate

# to generate ms simulations (& apply s*) where requiring stepwise segregating sites (val2 in the R script - informed by the array number)
#Rscript --vanilla /scratch/devel/hpawar/admix/sstar/scripts/troubleshoot.simul.R $ID 

# to generate ms simulations (& apply s*) using a higher divergence time for the split b/n E & W gorillas (& require stepwise seg sites in order to generate glms)
Rscript --vanilla /scratch/devel/hpawar/admix/sstar/scripts/high.divergence.simul.R $ID
