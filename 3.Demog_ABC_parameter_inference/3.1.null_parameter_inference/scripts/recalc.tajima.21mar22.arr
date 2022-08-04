#!/bin/bash
# @ job_name = abcnull_tajima
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/taj%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/taj%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 6:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 2,59,82,89,113,126,152,156,166,189,242,247,258,285,295,324,354,361,370,378,380,419,427,428,453,457,466,512,514,515,536,608,610,631,636,664,679


ID=$SLURM_ARRAY_TASK_ID

# array = 1-700
# 2h sufficient when using mclapply across 6 cores
#-----------------------------------------------------------------------------------------------------------------------

# Recalculate tajimas d for the null demog simulations
# regenerate these simulations, & only calculate the values for tajimas d - outputting the mean & sd per iter, & the object with the value per window.
#-----------------------------------------------------------------------------------------------------------------------
# this section worked (reps 1-700), but some reps did not give output => run the following section 
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.21mar22.R  /scratch/devel/hpawar/admix/abc/simul/log/tajima_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;

#exit

#-----------------------------------------------------------------------------------------------------------------------
# Sun 27 Mar 2022 18:45:48 CEST
# regenerate ms simns & recalc tajimas d for failed iter of the following reps
#2,3,89,126,152,160,166,189,242,258,285,294,324,343,354,380,386,419,427,428,433,453,457,466,512,515,536,595,608,610,664,679

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/recalc.tajima.failedreps.R  /scratch/devel/hpawar/admix/abc/simul/log/tajima_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;

exit
#-----------------------------------------------------------------------------------------------------------------------

