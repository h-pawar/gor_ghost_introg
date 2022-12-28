#!/bin/bash
# @ job_name = window_info
# @ initialdir = /scratch/devel/hpawar/admix/abc/emp.data/
# @ output = /scratch/devel/hpawar/admix/abc/emp.data/log/info_chr%j_%a.out
# @ error = /scratch/devel/hpawar/admix/abc/emp.data/log/info_chr%j_%a.err
# @ total_tasks = 1
# @ cpus_per_task = 4
# @ wall_clock_limit = 12:00:00
# @ class = lowprio
# @ requeue = 1
# @ array = 1-22

# array = 1-22

ID=$SLURM_ARRAY_TASK_ID

# annotate 40kb windows from the empirical data with how informative they are (amount of data/snps within the window)
	#  after filtering for repeats and mapability
		# then next step is filter by these windows - & only retain summary statistics in sufficiently informative windows
		#  heterozygosity, number of segregating sites, pairwise Fst, pi, Tajima's D)

#-----------------------------------------------------------------------------------------------------------------------

module load TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools java/latest

chrom=$ID
if [[ "$chrom" = 23 ]]; then chrom="X"; fi

Rscript --vanilla /scratch/devel/hpawar/admix/abc/simul/scripts/get_clbl.1jul.R $chrom
