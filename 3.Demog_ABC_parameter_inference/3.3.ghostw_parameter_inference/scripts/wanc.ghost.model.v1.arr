#!/bin/bash
# @ job_name = wghost_abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/ghost/w%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/ghost/w%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 429-543, 639-700

# array = 2-700 # send first for a test job

#mkdir -p /scratch/devel/hpawar/admix/abc/simul/log/ghost

ID=$SLURM_ARRAY_TASK_ID
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata
#-----------------------------------------------------------------------------------------------------------------------
# Sat 11 Sep 2021 11:35:35 CEST
# new simulations
# 1) change way of calculating segregating sites 
  # to match empirical data need to remove sites which are fixed across all gorilla populations (ie 1/1 across all pops)
  # b/c empirical data being mapped to human ref introduces a large number of such sites
# 2) do not simulate where t8<t7 (& t5 or t6 > t7) - ie need (t5/t6<t7, t7<t8)
# 3) fix the time of admixture b/n the extant lineages (atm high levels of confusion b/c timing of admixture is overlapping the divergence times b/n the sp)

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.v1.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/wghostabcsimul_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'wghost.abc.sim*' -execdir rm -- '{}' \;

#exit

#-----------------------------------------------------------------------------------------------------------------------
#Fri  8 Apr 2022 11:46:28 CEST
# final weighted median posteriors from final null abc parameter inference
R CMD BATCH --vanilla --slave "--args ${ID}"  /scratch/devel/hpawar/admix/abc/simul/scripts/wanc.ghost.model.1apr22.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/wghostabcsimul_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'wanc.ghost.abc.sim*' -execdir rm -- '{}' \;

exit 
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
