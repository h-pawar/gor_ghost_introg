#!/bin/bash
# @ job_name = pca_sstarskov
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/pca_sstarskov_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/pca_sstarskov_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 6:00:00
# @ class = normal
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID
# array=1

module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14

# Tue 26 Jul 2022 11:36:34 CEST
# generate PCAs for intersect regions (of 99% s* - strict_skov )
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/pca_sstarskov_26jul22.R