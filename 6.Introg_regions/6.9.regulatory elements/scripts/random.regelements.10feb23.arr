#!/bin/bash
# @ job_name = random_regelements
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log/random_regelements_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log/random_regelements_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 24:00:00
# @ class = normal
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID
# array=1

module load gcc/6.3.0 openssl/1.0.2q R/4.0.1


#mkdir -p /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/revisions_8feb23
# mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log


# Fri 10 Feb 2023 16:52:20 CET
# perform 100 iterations of random intersects - generating random regions of sufficient callable sites
# & assess proportion of bp in regulatory regions
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/revisions_8feb23/random.regelements.10feb23.R
