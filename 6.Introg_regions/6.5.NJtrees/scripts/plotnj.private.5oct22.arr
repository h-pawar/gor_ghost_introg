#!/bin/bash
# @ job_name = plot_nj
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/plotnj_privateintersect_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/plotnj_privateintersect_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 10:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-21

ID=$SLURM_ARRAY_TASK_ID

module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14

# Thu 29 Sep 2022 16:02:10 CEST
# generate basic plots of nj trees & bootstraps

m=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $1}')
SPE=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $2}')
spe=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $3}')

Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/plotnj.private.5oct22.R $m $SPE $spe
