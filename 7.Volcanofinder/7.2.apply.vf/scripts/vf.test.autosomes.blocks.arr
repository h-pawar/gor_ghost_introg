#!/bin/bash
# @ job_name = volcanofinder_blocks
# @ initialdir = /home/devel/hpawar
# @ output = /scratch/devel/hpawar/admix/volcanofinder/log/test/vf_blocks_23apr22_%j_%a.out
# @ error = /scratch/devel/hpawar/admix/volcanofinder/log/test/vf_blocks_23apr22_%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 168:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 481,483


ID=$SLURM_ARRAY_TASK_ID

#-----------------------------------------------------------------------------------------------------------------------
# Sat 23 Apr 2022 11:55:20 CEST
# apply volcanofinder to chr split into blocks

# cat /scratch/devel/hpawar/admix/volcanofinder/scripts/blocks_chr_23apr22.txt | wc -l
# ie 555 array jobs
# some reps completed at 24h, many did not
#-----------------------------------------------------------------------------------------------------------------------

#  col1 = block number (i 1..x) col2 = chr number (1..22), col3 = number of blocks x
BLOCK=$(sed -n ${ID}p /scratch/devel/hpawar/admix/volcanofinder/scripts/blocks_chr_23apr22.txt | awk -F ":" '{print $2}')
CHR=$(sed -n ${ID}p /scratch/devel/hpawar/admix/volcanofinder/scripts/blocks_chr_23apr22.txt | awk -F ":" '{print $3}')
NBLOCK=$(sed -n ${ID}p /scratch/devel/hpawar/admix/volcanofinder/scripts/blocks_chr_23apr22.txt | awk -F ":" '{print $4}')


FreqFile="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_"$CHR"_SF2_polym.input"
SpectFile="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/new_4_autosomes_sfs.input.txt"
OutFile="/scratch/devel/hpawar/admix/volcanofinder/output/4_"$CHR"_test"

/home/devel/hpawar/volcanofinder_v1.0/VolcanoFinder -big 1000 $FreqFile $SpectFile -1 1 1 $OutFile $BLOCK $NBLOCK
