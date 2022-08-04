#!/bin/bash
# @ job_name = abc_simul
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 6
# @ wall_clock_limit = 9:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 2-700

ID=$SLURM_ARRAY_TASK_ID

# MK - If half the data takes 40min, you can launch jobs with 51 iterations on 4 CPUs, which should take ~8.5h, so you can set 9h for the time. 
# array = 1-2000
#-----------------------------------------------------------------------------------------------------------------------

# Generate ms simulations for abc + output relevant summary statistics : heterozygosity, number of segregating sites, pairwise Fst
# using script /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.R 

#-----------------------------------------------------------------------------------------------------------------------

# whether need to load python here? only for sstar or also for ms?
#module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix 
#source venv2/bin/activate
#-----------------------------------------------------------------------------------------------------------------------

#module load gcc/6.3.0 R/4.0.1 tabix 
module load gcc/6.3.0 R/3.4.2 tabix # for v3

mkdir /dev/shm/mydata
cd /dev/shm/mydata

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.R  /scratch/devel/hpawar/admix/abc/simul/log/abcsimul_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;

#-----------------------------------------------------------------------------------------------------------------------

# try to output 3 extra stats - /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v2.R
# heterozygosity, number of segregating sites, pairwise Fst, Tajima's D, nucleotide diversity (pi), r2 (LD between pairs of SNPs).

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v2.R  /scratch/devel/hpawar/admix/abc/simul/log/abcsimul_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;
#-----------------------------------------------------------------------------------------------------------------------

# output 2 extra stats - /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v3.R
# heterozygosity, number of segregating sites, pairwise Fst, Tajima's D, nucleotide diversity (pi) - last 2 using pegas & the GT tables, rather htan popgenome & the LD-based stat r2 which were taking too long to calc

#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v3.R  /scratch/devel/hpawar/admix/abc/simul/log/abcsimul_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;

#-----------------------------------------------------------------------------------------------------------------------

# Mon 31 May 2021 16:41:02 CEST - used to generate first batch of 2000 simulations in dir /scratch/devel/hpawar/admix/abc/simul/test/out/
# amended script - to fix errors v3 was throwing
#R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v4.R  /scratch/devel/hpawar/admix/abc/simul/log/abcsimul_${ID}.log
#find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;


#exit

#-----------------------------------------------------------------------------------------------------------------------
# Sat 11 Sep 2021 11:35:35 CEST
# new simulations
# 1) change way of calculating segregating sites 
  # to match empirical data need to remove sites which are fixed across all gorilla populations (ie 1/1 across all pops)
  # b/c empirical data being mapped to human ref introduces a large number of such sites
# 2) do not simulate where t8<t7 (& t5 or t6 > t7) - ie need (t5/t6<t7, t7<t8)
# 3) fix the time of admixture b/n the extant lineages (atm high levels of confusion b/c timing of admixture is overlapping the divergence times b/n the sp)

R CMD BATCH --vanilla --slave "--args ${ID}" /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R  /scratch/devel/hpawar/admix/abc/simul/log/abcsimul_${ID}.log
find /dev/shm/mydata/ -type f -mtime +1 -name 'abc.sim*' -execdir rm -- '{}' \;


exit
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# MK:
# the last step is there to remove remaining simulations
# concerning the time, you must figure out how much time it approximately takes to run this script one time

#You will run the script 1000 independent times, each creating 102 random combinations of values, 
#each with a random number assigned. 
#Before doing so, I recommend running 1 job (or start with 1 full iteration) to see how long it takes, and adjusting the time in the submission script. 
#It should be approximately the time of one iteration*102/6.

# ie send test job with iter=1
#-----------------------------------------------------------------------------------------------------------------------
# MK
#2) Actually, the number of CPUs could be smaller than the number of threads, since each CPU may be able to run several threads (hyperthreading). 
#I am not entirely sure how well that works (Im not tech-bro enough for this...). 
#Probably it will not speed up all parts, but possibly some of it. I suggest you try the following:  cpus_per_task = 4 in the array job; 
#in the script for mclapply use mc.cores=6. 
#Then, try to run the script with iter=6 and measure how long that takes. 
#If it takes, say, 30min, then you could load jobs of 9h for 
#102 iterations (102iterations/6iterations*30min/60min=8.5h).

# already calling this in R script
# mc.cores=6 # but change the number of cpu cores requested to 4 (instead of requesting 6) ** 
#-----------------------------------------------------------------------------------------------------------------------



# /scratch/devel/hpawar/admix/abc/simul/scripts # scripts dir
# /scratch/devel/hpawar/admix/abc/simul/test # output for test job
# /scratch/devel/hpawar/admix/abc/simul/log # log dir
