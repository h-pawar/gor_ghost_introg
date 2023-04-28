#!/bin/bash
# @ job_name = gene_density_s*_skov_overlap_22oct22
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/gene_density_s*_skov_overlap_22oct22_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/gene_density_s*_skov_overlap_22oct22_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 24:00:00
# @ class = normal
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID
# array=1

module load gcc/6.3.0 openssl/1.0.2q R/4.0.1

# Sat 22 Oct 2022 12:17:00 CEST
# perform 100 iterations of random intersects - generating random regions of sufficient callable sites & calc gene density
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/random.gene.density.s*skov.overlap.21oct22.R
