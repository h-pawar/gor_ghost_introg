# Wed 29 Sep 2021 09:01:56 CEST
# have regenerated the stats - need to reprocess

#-----------------------------------------------------------------------------------------------------------------------

#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1) identify which windows had infinite fsts -> need to deduct their proportion*40000 from the total

files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", pattern="_informativewindows"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-stats_inform1
}

#replace NaN values with 0
#MK: If they are NaN (not a number), they are due to division by zero (I guess actually often 0/0) and might be corrected to 0 instead. 
find.nas_fun1<-function(x) {
for (j in 1:length(genomewide[[x]])){
genomewide[[x]][[j]][[1]]<-as.numeric(genomewide[[x]][[j]][[1]])
}
out_1<-list()
for (j in 1:length(genomewide[[x]])){
out_1[[j]]<-lapply(genomewide[[x]][[j]], function(y) replace(y, !is.finite(y), 0))
}
return(out_1)
}



finalwindows=list()
for (i in 1:length(genomewide)){
finalwindows[[i]]<-find.nas_fun1(i)
}

#-----------------------------------------------------------------------------------------------------------------------

format_stats_fun<-function(x,y) {
otest=list()
for (i in 1:length(finalwindows[[x]])){
otest[[i]]<-t(as.data.frame(finalwindows[[x]][[i]][[y]]))
}
outest<-data.frame(matrix(unlist(otest), nrow=length(otest), byrow=TRUE),stringsAsFactors=FALSE)
return(outest)
}

stats_asdf_fun<-function(x){
stats_o=list()
for (i in 1:22){
stats_o[[i]]<-format_stats_fun(i,x)
}
ptest<-do.call(rbind, stats_o)
return(ptest)
}

fst_df<-stats_asdf_fun(4)

inf_fsts=list()
for (i in 1:ncol(fst_df)){
inf_fsts[[i]]<-which(fst_df[,i] > 1)
}

fst_df1<-fst_df[-c(sort(unlist(inf_fsts))),]

#-----------------------------------------------------------------------------------------------------------------------
# Tue 28 Sep 2021 11:46:21 CEST
# 2) read in files w newly calculated stats
#- number of population-wise fixed sites and the number of population-wise segregating sites; 
# - the fixed sites per individual
# generated using
    #/scratch/devel/hpawar/admix/abc/simul/scripts/segs.emp.27sep.R  # called by 
    #/scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr
#then ran interactively - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.segs.by.clbl.27sep.R
#-----------------------------------------------------------------------------------------------------------------------

segsfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", pattern="_informativewindows"), sep = "")

segsgenomewide=list()
for (i in 1:length(segsfiles)){
load(file=segsfiles[[i]],verbose=T)
segsgenomewide[[i]]<-stats_inform1
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# A) number of population-wise fixed sites

format_fixedsitesperid_fun<-function(x,y,z) {
htest=list()
for (i in 1:length(segsgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(segsgenomewide[[x]][[i]][[2]][z,][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


fixedsitesperid_asdf_fun<-function(y,z){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_fixedsitesperid_fun(i,y,z)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}



fixedsitesperid_wl<-fixedsitesperid_asdf_fun(1,1)
# remove stats from windows with infinite fsts
fixedsitesperid_wl1<-fixedsitesperid_wl[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_wc<-fixedsitesperid_asdf_fun(2,1)
fixedsitesperid_wc1<-fixedsitesperid_wc[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_el<-fixedsitesperid_asdf_fun(3,1)
fixedsitesperid_el1<-fixedsitesperid_el[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_em<-fixedsitesperid_asdf_fun(4,1)
fixedsitesperid_em1<-fixedsitesperid_em[-c(sort(unlist(inf_fsts))),]


# normalise by data coverage 
rbind(sum(fixedsitesperid_wl1/767730079)*1000,
sum(fixedsitesperid_wc1/767730079)*1000,
sum(fixedsitesperid_el1/767730079)*1000,
sum(fixedsitesperid_em1/767730079)*1000)

#-----------------------------------------------------------------------------------------------------------------------
# B) number of population-wise segregating sites

segsitesperid_wl<-fixedsitesperid_asdf_fun(1,2)
segsitesperid_wl1<-segsitesperid_wl[-c(sort(unlist(inf_fsts))),]
segsitesperid_wc<-fixedsitesperid_asdf_fun(2,2)
segsitesperid_wc1<-segsitesperid_wc[-c(sort(unlist(inf_fsts))),]
segsitesperid_el<-fixedsitesperid_asdf_fun(3,2)
segsitesperid_el1<-segsitesperid_el[-c(sort(unlist(inf_fsts))),]
segsitesperid_em<-fixedsitesperid_asdf_fun(4,2)
segsitesperid_em1<-segsitesperid_em[-c(sort(unlist(inf_fsts))),]

rbind(sum(segsitesperid_wl1/767730079)*1000,
sum(segsitesperid_wc1/767730079)*1000,
sum(segsitesperid_el1/767730079)*1000,
sum(segsitesperid_em1/767730079)*1000)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2.2) fixed sites per individual

format_fixedperid_fun<-function(x,y) {
htest=list()
for (i in 1:length(segsgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(segsgenomewide[[x]][[i]][[3]][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


fixedperid_asdf_fun<-function(y){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_fixedperid_fun(i,y)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}

fixedperid_wl<-fixedperid_asdf_fun(1)
fixedperid_wc<-fixedperid_asdf_fun(2)
fixedperid_el<-fixedperid_asdf_fun(3)
fixedperid_em<-fixedperid_asdf_fun(4)

fixedperid_wl1<-fixedperid_wl[-c(sort(unlist(inf_fsts))),]
fixedperid_wc1<-fixedperid_wc[-c(sort(unlist(inf_fsts))),]
fixedperid_el1<-fixedperid_el[-c(sort(unlist(inf_fsts))),]
fixedperid_em1<-fixedperid_em[-c(sort(unlist(inf_fsts))),]
#-----------------------------------------------------------------------------------------------------------------------

# 3) sum of all fixed sites divided by the sum of all sites with data
# *1000 for fixed sites per kb

cbind(mean((colSums(fixedperid_wl1)/767730079))*1000,sd((colSums(fixedperid_wl1)/767730079))*1000)
cbind(mean((colSums(fixedperid_el1)/767730079))*1000,sd((colSums(fixedperid_el1)/767730079))*1000)
cbind(mean((colSums(fixedperid_em1)/767730079))*1000,sd((colSums(fixedperid_em1)/767730079))*1000)
sum(fixedperid_wc1/767730079)*1000

#-----------------------------------------------------------------------------------------------------------------------
