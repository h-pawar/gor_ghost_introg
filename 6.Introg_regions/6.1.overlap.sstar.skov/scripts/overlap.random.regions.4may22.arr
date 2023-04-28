#!/bin/bash
# @ job_name = random_overlap
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/random_overlap_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/random_overlap_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 8:00:00
# @ class = normal
# @ requeue = 1
# @ array = 1

#ID=$SLURM_ARRAY_TASK_ID
# array=1

module load gcc/6.3.0 openssl/1.0.2q R/4.0.1


# perform 100 iterations of random intersects

# Wed 20 Jul 2022 17:22:17 CEST
Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/overlap.random.regions.20jul22.R
#-----------------------------------------------------------------------------------------------------------------------

