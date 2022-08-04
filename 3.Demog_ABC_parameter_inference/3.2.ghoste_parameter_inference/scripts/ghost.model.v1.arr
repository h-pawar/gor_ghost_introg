#!/bin/bash
# @ job_name = ghost_abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/ghost/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/ghost/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-700

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

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.v1.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/ghostabcsimul_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'ghost.abc.sim*' -execdir rm -- '{}' \;


#exit

#-----------------------------------------------------------------------------------------------------------------------
#Thu  9 Dec 2021 11:33:47 CET
# re-infer model inference with ghost admixture
# i.e. infer all (or more) parameters instead of fixing them. 
# new ghost -> e_anc simulations
#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.v2.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/ghostabcsimulv2_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'ghost.abc.v2.sim*' -execdir rm -- '{}' \;


#exit
#-----------------------------------------------------------------------------------------------------------------------
#Tue  1 Mar 2022 10:13:28 CET
# infer  ghoste model parameter inference - building from posteriors from this new null model of 28feb22 (with logit transf in the abc) 
# has not output tajimas d for EM - not sure why Thu 17 Mar 2022 00:34:06 CET
# Tue  5 Apr 2022 15:34:56 CEST - amended with final null posteriors from paramteter inference
R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.1mar22.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/ghostabcsimul_1mar22_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'ghost.abc.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#Fri 11 Mar 2022 10:39:12 CET
# think have amended the error from *1mar22.R
# ghoste model parameter inference - building from posteriors from this new null model of 28feb22 (with logit transf in the abc) 

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/ghost.model.11mar22.R  /scratch/devel/hpawar/admix/abc/simul/log/ghost/ghostabcsimul_11mar22_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'ghost.abc.sim*' -execdir rm -- '{}' \;


#exit
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
