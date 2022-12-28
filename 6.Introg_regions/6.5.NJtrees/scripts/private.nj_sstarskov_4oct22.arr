#!/bin/bash
# @ job_name = private.nj_sstarskov
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/private.nj_sstarskov_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/private.nj_sstarskov_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 6:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-21

ID=$SLURM_ARRAY_TASK_ID

module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14

# Tue  4 Oct 2022 10:01:25 CEST
# generate nj trees of private regions (intersect sstar-skov regions which are found only in one id of the pop)


m=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $1}')
SPE=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $2}')
spe=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $3}')

Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/private.nj_sstarskov_4oct22.R $m $SPE $spe
