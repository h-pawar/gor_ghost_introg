#!/bin/bash
# @ job_name = pls
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/pls_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/simul/log/modelcomp/pls_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 12:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1

#Tue 15 Feb 2022 14:16:05 CET
# perform box-cox transformation, then pls of the summary statistics
	# was taking a long time to run the pls step interactively -> sending as a job

#ID=$SLURM_ARRAY_TASK_ID

#-----------------------------------------------------------------------------------------------------------------------

module load R/4.0.1
#Rscript --vanilla   /scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/decorrelatestats/uncorrelate.stats.26jan22.R

#-----------------------------------------------------------------------------------------------------------------------
#Wed 16 Feb 2022 14:53:19 CET
# run the pls step with validation
# have already run the pls step without validation, by running interactively :  /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/uncorrelate.stats.26jan22.R
Rscript --vanilla   /scratch/devel/hpawar/admix/abc/simul/scripts/modelcomp/decorrelatestats/pls.16feb22.R
#-----------------------------------------------------------------------------------------------------------------------

