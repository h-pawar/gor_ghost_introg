#!/bin/bash
# @ job_name = nj_sstarskov
# @ initialdir = /scratch/devel/hpawar/
# @ output = /scratch/devel/hpawar/admix/sstar/log/nj_sstarskov_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/sstar/log/nj_sstarskov_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 8:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 21

ID=$SLURM_ARRAY_TASK_ID

module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14

# Wed 29 Jun 2022 10:00:31 CEST
# generate nj trees for overlap regions (s*-skov)

m=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $1}')
SPE=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $2}')
spe=$(sed -n ${ID}p  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj.comp.list | awk -F ":" '{print $3}')

Rscript --vanilla  /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/nj_sstarskov_26jul22.R $m $SPE $spe
