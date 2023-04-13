#!/bin/bash
# @ job_name = ghoste_abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/17feb23%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/17feb23%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-290


#-----------------------------------------------------------------------------------------------------------------------
# E: 21623 meaningful simulations generated (t_archintrog < t8)
# W: 13941 meaningful simulations generated (t_archintrog < t8)


#(21623/35700)
#[1] 0.6056863

#(1- 0.6056863)*35700
#[1] 14077
#> ((1- 0.6056863)*35700)/51
#[1] 276.0196

# -> target 290 simns for ghoste
#-----------------------------------------------------------------------------------------------------------------------


ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------
# Fri 17 Feb 2023 11:55:48 CET
# ABC ghost parameter inference - sample all parameters from priors
#1) add condition to retain only simulations where archaic introgression time < species split time → 
# generate more simulations to get to ~35.7k (to make revised parameter inference complete)

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/revisions_feb23/ghoste.rev.17feb23.R /scratch/devel/hpawar/admix/abc/simul/test/ghost/revisions_8feb23/log/ghosterev_simul_17feb23_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'ghoste.rev.abc.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------
