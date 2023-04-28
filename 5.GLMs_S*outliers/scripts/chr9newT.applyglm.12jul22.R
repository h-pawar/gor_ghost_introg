# Tue 12 Jul 2022 14:27:49 CEST
# S* calculated on final dataset 
# now apply glm to data, reinfer fragments & replot distribution of fragments genomewide

#-----------------------------------------------------------------------------------------------------------------------
# module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix
#-----------------------------------------------------------------------------------------------------------------------
load(file="/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/glm",verbose=T)
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


  # update to per subspecies comparisons
# 1) WL outgroup, EL ingroup (GGG outg, GBG ing) = "GBG"
# 2) WL outgroup, EM ingroup (GGG outg, GBB ing) = "GBB"
# 3) EL outgroup, WL ingroup (GBG outg, GGG ing) = "GGG"
# 4) EL outgroup, WC ingroup (GBG outg, GGD ing) = "GGD"


cn1<-list(c("GBG","GBB","GGG","GGD"))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


applyglm_function2<-function(nput,ci) {

starout<-list()
for (chrom in (chroms)) {
starout[[chrom]]<-read.table(paste("/scratch/devel/hpawar/admix/sstar/out.subsp.comp/",cn1[[1]][[nput]],"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] 
}
# read in data for s* applied to empirical data
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

allstars<-do.call(rbind,starperind.gbb)

# jitter is applied to the numbers of segregating sites (column 3 & 4), to make data less stepwise
allstars[,3]<-as.numeric(jitter(allstars[,3]))

## either choose a subset of windows (for thinning, maybe for plotting etc), or everything
subse=1:nrow(allstars)

# create a data frame of the same segregating sites:
newdatA=data.frame(sS=allstars[,3])

# get the predicted S* for each window given the segregating sites (using the glm for the 99.5% CI)

  # predict S* values at a given CI:
out.pred<-predict(mods[[nput]][[ci]],newdata=newdatA,type="response") 

# final window ABC
save(out.pred,file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/predictions_",cn1[[1]][[nput]],"_",ci,"_chr9newT_from_glm_25apr22",sep=""))    


ldif_mod1<-out.pred[]-allstars[,1]

# which windows are outside the expectation for the XX% CI for both comparisons?  
fset3_mod1<-which( allstars[,1]>out.pred[])

save(fset3_mod1,file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/fset_",cn1[[1]][[nput]],"_",ci,"_chr9newT_from_glm_25apr22",sep=""))    

test<- allstars[fset3_mod1,]

save(test,file=paste("/scratch/devel/hpawar/admix/sstar/simul_postabc/stepwisesegs/final_8apr22/check_22apr22/sstaroutliers_",cn1[[1]][[nput]],"_",ci,"_chr9newT_from_glm_25apr22",sep=""))    


# outlier windows / total windows
length(fset3_mod1)/nrow(ldif_mod1)


# ranges for the n_snps & S* vals
xr=range(allstars[,3]);yr=range(allstars[,1])

# take the mean & sd of the S* vals, & the n_snps for the empirical data

return(rbind(length(fset3_mod1)/nrow(ldif_mod1), 
cbind(mean(allstars[,1]), sd(allstars[,1])),
cbind(mean(allstars[,3]), sd(allstars[,3])),
xr,yr
  ))

}
#-----------------------------------------------------------------------------------------------------------------------
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


