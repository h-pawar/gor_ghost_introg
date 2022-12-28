# Tue 12 Jul 2022 14:27:49 CEST
# S* calculated on chr 9 with newly processed tumani sample - script - /scratch/devel/hpawar/admix/sstar/scripts/sstar.newT.chr9.arr
# now apply glm to data (with chr 9 plus newly processed tumani), reinfer fragments & replot distribution of fragments genomewide
	# following procedure of when reinferred for chr9 without tumani - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/chr9noT.applyglm.8jun22.R
	# & following /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/check.postabc.glm.25apr22.R

#-----------------------------------------------------------------------------------------------------------------------
## check output of S* on newly processed data looks as expected (compared to previous)
#[hpawar@login2 ~]$ head /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBG_newT_9.star
#chrom	winstart	winend	n_snps	n_ind_snps	n_region_ind_snps	ind_id	pop	s_star	num_s_star_snps	s_star_snps	hap_1_window_pval	hap_2_window_pval	hap_1_window_match_pct	hap_2_window_match_pct	hap_1_window_pval_local	hap_2_window_pval_local	hap_1_window_match_pct_local	hap_2_window_match_pct_local	hap_1_window_pval_table	hap_2_window_pval_table	hap_1_window_match_pct_table	hap_2_window_match_pct_table	hap_1_window_match_N_table	hap_2_window_match_N_table	hap_1_window_match_f_table	hap_2_window_match_f_table	hap_1_window_match_len_table	hap_2_window_match_len_table	hap_1_window_match_mapped_table	hap_2_window_match_mapped_table	hap_1_window_match_mh_table	hap_2_window_match_mh_table	hap_1_window_match_sfs_table	hap_2_window_match_sfs_table	hap_1_window_match_orig_sfs_table	hap_2_window_match_orig_sfs_table	hap_1_window_match_arc_table	hap_2_window_match_arc_table	hap_1_window_match_tot_sites_table	hap_2_window_match_tot_sites_table	hap_1_s_start	hap_1_s_end	hap_2_s_start	hap_2_s_end	s_start	s_end	ref_nocals	n_s_star_snps_hap1	n_s_star_snps_hap2	s_star_haps	callable_bases
#chr9	180000	220000	49	2	48	Gbg-Tumani	GBG	0	0	.	1.0	1.0	0.0	0.0	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0	0	0	0	NA	NA	36904
#chr9	180000	220000	49	3	49	Gorilla_beringei_graueri-9732_Mkubwa	GBG	10492	2	204920,210411.0	1.0	0.0	0.0	1.0	1.0	0.0	0.0	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0	0	204920	210412	204920	210412	0	0	2	2,2	36904

#[hpawar@cnb1 ~]$ head /scratch/devel/hpawar/admix/sstar/out.subsp.comp/GBG_22.star
#chrom	winstart	winend	n_snps	n_ind_snps	n_region_ind_snps	ind_id	pop	s_star	num_s_star_snps	s_star_snps	hap_1_window_pval	hap_2_window_pval	hap_1_window_match_pct	hap_2_window_match_pct	hap_1_window_pval_local	hap_2_window_pval_local	hap_1_window_match_pct_local	hap_2_window_match_pct_local	hap_1_window_pval_table	hap_2_window_pval_table	hap_1_window_match_pct_table	hap_2_window_match_pct_table	hap_1_window_match_N_table	hap_2_window_match_N_table	hap_1_window_match_f_table	hap_2_window_match_f_table	hap_1_window_match_len_table	hap_2_window_match_len_table	hap_1_window_match_mapped_table	hap_2_window_match_mapped_table	hap_1_window_match_mh_table	hap_2_window_match_mh_table	hap_1_window_match_sfs_table	hap_2_window_match_sfs_table	hap_1_window_match_orig_sfs_table	hap_2_window_match_orig_sfs_table	hap_1_window_match_arc_table	hap_2_window_match_arc_table	hap_1_window_match_tot_sites_table	hap_2_window_match_tot_sites_table	hap_1_s_start	hap_1_s_end	hap_2_s_start	hap_2_s_end	s_start	s_end	ref_nocals	n_s_star_snps_hap1	n_s_star_snps_hap2	s_star_haps	callable_bases
#chr22	16200000	16240000	3	0	3	Gorilla_beringei_graueri-9732_Mkubwa	GBG	0	0	1.0	1.0	None	None	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0	0	NA	NA	39	0	0	.	32413
#chr22	16200000	16240000	3	0	3	Gorilla_beringei_graueri-A929_Kaisi	GBG	0	0	1.0	1.0	None	None	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0	0	NA	NA	39	0	0	.	32413

#-----------------------------------------------------------------------------------------------------------------------

# paths to chr 9 sstar fragments inferred when calculating with newly processed tumani in the vcf
#ls -lhtr /scratch/devel/hpawar/admix/sstar/out.subsp.comp/
#-rw-r--r-- 1 hpawar devel 8.2M Jul  8 13:26 GBG_newT_9.star
#-rw-r--r-- 1 hpawar devel  12M Jul  8 13:42 GBB_newT_9.star
#-rw-r--r-- 1 hpawar devel  34M Jul  8 14:19 GGG_newT_9.star
#-rw-r--r-- 1 hpawar devel 1.2M Jul  8 14:26 GGD_newT_9.star
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix
#-----------------------------------------------------------------------------------------------------------------------
load(file="/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",verbose=T)   # also need the predictions for the function applied to empirical data below
	# this is output from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/check.postabc.glm.25apr22.R

options(scipen=100)
library(mgcv)
vals=seq(15,800,5)

# create a data frame of the same segregating sites:
newdatA2<-data.frame(sS=vals)


predictions<-list()
# for the 4 scenarios, for the 3rd CI - 99.5% CI
for (i in 1:4){
  # predict S* values at a given CI:
out.pred<-predict(mods[[i]][[3]],newdata=newdatA2,type="response") 
predictions[[i]]<-out.pred
}

#-----------------------------------------------------------------------------------------------------------------------



# read in all comparisons for all chromosomes, and process them
chroms=1:23
# will need to specify either GB or GG
#cn1<-list(c("GB","GG"))

  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))

#-----------------------------------------------------------------------------------------------------------------------
# run interactively - to check order of chromosomes being read in - all fine

#starout<-list()
#for (chrom in (chroms)) {
#starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 
#}

# head(starout[[9]])
#  s_star n_snps num_s_star_snps callable_bases hap_1_window_pval
#1  10492     47               2          36904                 1
#2      0     47               0          36904                 1
#3  10492     47               2          36904                 1
#4      0     47               0          36904                 1
#5      0     47               0          36904                 1
#6      0     47               0          36904                 1
#  hap_2_window_pval                                 ind_id chrom winstart
#1                 1   Gorilla_beringei_graueri-9732_Mkubwa  chr9   180000
#2                 1    Gorilla_beringei_graueri-A929_Kaisi  chr9   180000
#3                 1 Gorilla_beringei_graueri-A967_Victoria  chr9   180000

#nrow(starout[[9]])
#[1] 9684

#ls -lhtr /scratch/devel/hpawar/admix/sstar/out.subsp.comp/

#starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

#GGD_newT_9.star

#> head(starout[[9]])
#  s_star n_snps num_s_star_snps callable_bases hap_1_window_pval
#1  10492     47               2          36904                 1
#2      0     47               0          36904                 1
#3  10492     47               2          36904                 1
#4      0     47               0          36904                 1
#5      0     47               0          36904                 1
#6      0     47               0          36904                 1
#  hap_2_window_pval                                 ind_id chrom winstart
#1                 1   Gorilla_beringei_graueri-9732_Mkubwa  chr9   180000
#2                 1    Gorilla_beringei_graueri-A929_Kaisi  chr9   180000
#3                 1 Gorilla_beringei_graueri-A967_Victoria  chr9   180000
#4                 1         Gorilla_beringei_graueri-Dunia  chr9   180000
#5                 1       Gorilla_beringei_graueri-Itebero  chr9   180000
#6                 1       Gorilla_beringei_graueri-Mukokya  chr9   180000
#> nrow(starout[[9]])
#[1] 29016
	# now much more fragments inferred
#-----------------------------------------------------------------------------------------------------------------------


applyglm_function2<-function(nput,ci) {

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 
}
# read in data for s* applied to chr 9 for target individuals minus tumani (whose chr 9 had sparse data issues)
starout[[9]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_newT_9.star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 

staroutgbb<-do.call(rbind,starout)

#get values per individual - for each comparison
staroutgbb[,8]<-paste(staroutgbb[,8],staroutgbb[,9])
starperind.gbb<-list()
# col 7 = ind_id
indiv.gbb<-unique(staroutgbb[,7])

# here, I only want to use windows where there is enough data (33,333 bp out of 40,000, with at least 30 segregating sites)
for (ind in (1:length(indiv.gbb))) {
starperind.gbb[[ind]]<- staroutgbb[which(staroutgbb[,7]==indiv.gbb[ind] & staroutgbb[,4]>=33333 &staroutgbb[,2]>=30),c(1,8,2)] 
}

# is this correct, b/c haven't added identifiers of which ids each df of the list belongs to..  **
allstars<-do.call(rbind,starperind.gbb)

#the jitter is applied to the numbers of segregating sites (column 3 & 4), to avoid too many exact steps in the model. # ie for col3  - the n_snps
# jitter to make data less stepwise
allstars[,3]<-as.numeric(jitter(allstars[,3]))

## either choose a subset of windows (for thinning, maybe for plotting etc), or everything
#    subse=sample(1:nrow(allstars),50000)
subse=1:nrow(allstars)

# create a data frame of the same segregating sites:
newdatA=data.frame(sS=allstars[,3])

# get the predicted S* for each window given the segregating sites (using the glm for the 99.5% CI)
# for the 1 scenario relevant for this ingroup & for the 3rd CI - 99.5% CI -> for GBB mods[[2]]
  # predict S* values at a given CI:
out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response") 

#return(out.pred)

# write out these predictions_gb
#save(out.pred,file="/scratch/devel/hpawar/admix/sstar/out.subsp.comp/predictions_gbb_fromglm")    
#save(out.pred,file=paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/predictions_",cn1[[1]][[nput]],"_",ci,"_fromglm",sep=""))    


# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/stepwisesegs/out.subsp.comp

# final window ABC
save(out.pred,file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/predictions_",cn1[[1]][[nput]],"_",ci,"_chr9newT_from_glm_25apr22",sep=""))    


#}



ldif_mod1<-out.pred[]-allstars[,1]

# which windows are outside the expectation for the XX% CI for both comparisons?  
fset3_mod1<-which( allstars[,1]>out.pred[])

#save(fset3_mod1,file=paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/fset_",cn1[[1]][[nput]],"_",ci,"_fromglm",sep=""))    

save(fset3_mod1,file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/fset_",cn1[[1]][[nput]],"_",ci,"_chr9newT_from_glm_25apr22",sep=""))    


test<- allstars[fset3_mod1,]
#save(test,file=paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/sstaroutliers_",cn1[[1]][[nput]],"_",ci,"_fromglm",sep=""))    

save(test,file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/sstaroutliers_",cn1[[1]][[nput]],"_",ci,"_chr9newT_from_glm_25apr22",sep=""))    


# outlier windows / total windows
length(fset3_mod1)/nrow(ldif_mod1)
# [1] 0.2219823

# ranges for the n_snps & S* vals
xr=range(allstars[,3]);yr=range(allstars[,1])

# take the mean & sd of the S* vals, & the n_snps for the empirical data

return(rbind(length(fset3_mod1)/nrow(ldif_mod1), 
cbind(mean(allstars[,1]), sd(allstars[,1])),
cbind(mean(allstars[,3]), sd(allstars[,3])),
xr,yr
  ))

}

# CI 95%
applyglm_function2(1,1)
applyglm_function2(2,1)
applyglm_function2(3,1)
applyglm_function2(4,1)

# CI 99%
applyglm_function2(1,2)
applyglm_function2(2,2)
applyglm_function2(3,2)
applyglm_function2(4,2)

# CI 99.5%
applyglm_function2(1,3)
applyglm_function2(2,3)
applyglm_function2(3,3)
applyglm_function2(4,3)
#-----------------------------------------------------------------------------------------------------------------------
# output =
 # outlier windows / total windows
 # mean, sd of S* vals 
 # mean, sd of n_snps
 # ranges for n_snps
 # ranges for S* vals
# (for each of the 4 scenarios)
#-----------------------------------------------------------------------------------------------------------------------
# CI 95%
> # CI 95%

> applyglm_function2(1,1)
                            [,1]                       [,2]
                      0.02960561                 0.02960561
       -778698481709949.75000000 84744495606943616.00000000
                    103.29564902                52.33689692
xr                   29.80000351              1391.07975604
yr -9223372036854775808.00000000            372874.00000000
> applyglm_function2(2,1)
                            [,1]                       [,2]
                      0.04016825                 0.04016825
       -891641476640233.25000000 90681611581894480.00000000
                    105.47134030                53.24566474
xr                   29.80003626              1403.06617945
yr -9223372036854775808.00000000            362993.00000000
> applyglm_function2(3,1)
                            [,1]                       [,2]
                      0.02501515                 0.02501515
       -103826464148122.95312500 30945432233388544.00000000
                    103.29600850                52.33699712
xr                   29.80003681              1391.16441226
yr -9223372036854775808.00000000           1380657.00000000
> applyglm_function2(4,1)
               [,1]              [,2]
        0.000605273       0.000605273
    91358.609670126   49579.397937034
       60.723441992      32.914622299
xr     29.800802485    1141.055797948
yr -10000.000000000 1305388.000000000
> # CI 99%
> applyglm_function2(1,2)
                            [,1]                       [,2]
                      0.01643928                 0.01643928
       -778698481709949.75000000 84744495606943616.00000000
                    103.29572912                52.33691499
xr                   29.80012279              1391.19935805
yr -9223372036854775808.00000000            372874.00000000
> applyglm_function2(2,2)
                            [,1]                       [,2]
                      0.02362684                 0.02362684
       -891641476640233.25000000 90681611581894480.00000000
                    105.47158217                53.24578732
xr                   29.80007973              1403.18516068
yr -9223372036854775808.00000000            362993.00000000
> applyglm_function2(3,2)
                            [,1]                       [,2]
                      0.01285396                 0.01285396
       -103826464148122.95312500 30945432233388544.00000000
                    103.29567485                52.33677417
xr                   29.80010493              1391.18006993
yr -9223372036854775808.00000000           1380657.00000000
> applyglm_function2(4,2)
               [,1]              [,2]
        0.000124615       0.000124615
    91358.609670126   49579.397937034
       60.723836479      32.914451066
xr     29.800379573    1141.022898366
yr -10000.000000000 1305388.000000000
> 
> # CI 99.5%
> applyglm_function2(1,3)
                            [,1]                       [,2]
                      0.01113447                 0.01113447
       -778698481709949.75000000 84744495606943616.00000000
                    103.29562457                52.33695828
xr                   29.80019603              1391.17560884
yr -9223372036854775808.00000000            372874.00000000
> applyglm_function2(2,3)
                            [,1]                       [,2]
                      0.01608432                 0.01608432
       -891641476640233.25000000 90681611581894480.00000000
                    105.47150888                53.24579405
xr                   29.80000503              1403.19670726
yr -9223372036854775808.00000000            362993.00000000
> applyglm_function2(3,3)
                             [,1]                        [,2]
                      0.008023816                 0.008023816
       -103826464148122.953125000 30945432233388544.000000000
                    103.295814712                52.336930128
xr                   29.800073463              1391.195418219
yr -9223372036854775808.000000000           1380657.000000000
> applyglm_function2(4,3)
                 [,1]                [,2]
        0.00007120859       0.00007120859
    91358.60967012621   49579.39793703419
       60.72445663555      32.91514897781
xr     29.80010169400    1141.16574183805
yr -10000.00000000000 1305388.00000000000


# will need to update spreadsheet **
#-----------------------------------------------------------------------------------------------------------------------

# original vals (with old processing of tumani chr 9)

# Summary : 
# window-based ABC
  # CIs:.      95%    99%  99.5%
# gbg (el i) - 2.98% 1.65% 1.12%
# gbb (em i) - 4.05% 2.38% 1.62%
# ggg (wl i) - 2.52% 1.29% 0.80%
# ggd (wc i) - 0.060% 0.012% 0.007%


#-----------------------------------------------------------------------------------------------------------------------

# new vals with newly processed tumani chr 9 data

# Summary : 
# window-based ABC
  # CIs:.      95%    99%  99.5%
# gbg (el i) - 2.96% 1.64% 1.11%
# gbb (em i) - 4.02% 2.36% 1.60%
# ggg (wl i) - 2.50% 1.29% 0.80%
# ggd (wc i) - 0.060% 0.012% 0.007%

