#!/bin/bash
# @ job_name = random_setype
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log/random_setype_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/log/random_setype_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 12:00:00
# @ class = normal
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID

module load gcc/6.3.0 openssl/1.0.2q R/4.0.1

#Mon 20 Mar 2023 14:24:43 CET
# perform 100 iterations of random intersects - generating random regions of sufficient callable sites
# calculate types of strong enhancers in the random regions
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/revisions_8feb23/random.mg_se_locations_20mar23.R
