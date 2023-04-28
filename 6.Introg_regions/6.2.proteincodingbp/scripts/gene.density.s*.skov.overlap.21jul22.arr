#!/bin/bash
# @ job_name = gene_density_skov
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/gene_density_s*_skov_overlap_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/gene_density_s*_skov_overlap_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 4:00:00
# @ class = normal
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID
# array=1

module load gcc/6.3.0 openssl/1.0.2q R/4.0.1

# Wed  4 May 2022 19:20:52 CEST
# calculate gene density introgressed regions (s*-skov overlap)
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/gene.density.s*.skov.overlap.21jul22.R
