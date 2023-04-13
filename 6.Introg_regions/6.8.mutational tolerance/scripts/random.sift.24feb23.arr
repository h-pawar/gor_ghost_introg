#!/bin/bash
# @ job_name = random_gerp
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log/random_sift_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log/random_sift_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 24:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID

 module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
 
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/revisions_8feb23/random.sift.24feb23.R


