#!/usr/bin/env
#module load R/4.0.1

#Fri 21 Jan 2022 17:33:15 GMT
# which of the summary stats are correlated with each other?

#-----------------------------------------------------------------------------------------------------------------------

library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)

# read in for empirical data
files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", pattern="_informativewindows"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-stats_inform1
}


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

# summary stats for all retained windows (genome-wide) - for fst, pi, tajimasd
fst_df<-stats_asdf_fun(4)
pi_df<-stats_asdf_fun(5)
tajima_df<-stats_asdf_fun(6)

inf_fsts=list()
for (i in 1:ncol(fst_df)){
inf_fsts[[i]]<-which(fst_df[,i] > 1)
}


#-----------------------------------------------------------------------------------------------------------------------

fst_df1<-fst_df[-c(sort(unlist(inf_fsts))),]

pi_df1<-pi_df[-c(sort(unlist(inf_fsts))),]
tajima_df1<-tajima_df[-c(sort(unlist(inf_fsts))),]

#-----------------------------------------------------------------------------------------------------------------------

# 2) read in heterozygosity files 
hetfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/", pattern="_informativewindows"), sep = "")


hetgenomewide=list()
for (i in 1:length(hetfiles)){
load(file=hetfiles[[i]],verbose=T)
hetgenomewide[[i]]<-stats_inform1
}


# 2.1) calc heterozygosity per individual

format_hetperid_fun<-function(x,y) {
htest=list()
for (i in 1:length(hetgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(hetgenomewide[[x]][[i]][[2]][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


hetperid_asdf_fun<-function(y){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_hetperid_fun(i,y)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}

hetperid_wl<-hetperid_asdf_fun(1)
hetperid_wc<-hetperid_asdf_fun(2)
hetperid_el<-hetperid_asdf_fun(3)
hetperid_em<-hetperid_asdf_fun(4)

# removing windows remvoed for the other summary stats b/c infinite fsts
hetperid_wl1<-hetperid_wl[-c(sort(unlist(inf_fsts))),]
hetperid_wc1<-hetperid_wc[-c(sort(unlist(inf_fsts))),]
hetperid_el1<-hetperid_el[-c(sort(unlist(inf_fsts))),]
hetperid_em1<-hetperid_em[-c(sort(unlist(inf_fsts))),]


# 3) sum of all heterozygous sites divided by the sum of all sites with data
# *1000 for het per kb

het_mu1<-c(
mean((colSums(hetperid_wl1)/767730079))*1000,
sum(hetperid_wc1/767730079)*1000,
mean((colSums(hetperid_el1)/767730079))*1000,
mean((colSums(hetperid_em1)/767730079))*1000
    )

het_sd1<-c(
sd((colSums(hetperid_wl1)/767730079))*1000,
sd((colSums(hetperid_el1)/767730079))*1000,
sd((colSums(hetperid_em1)/767730079))*1000
    )

#-----------------------------------------------------------------------------------------------------------------------

#- number of population-wise fixed sites and the number of population-wise segregating sites; & - the fixed sites per individual
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



pop_fixedsites<-rbind(sum(fixedsitesperid_wl1/767730079)*1000,
sum(fixedsitesperid_wc1/767730079)*1000,
sum(fixedsitesperid_el1/767730079)*1000,
sum(fixedsitesperid_em1/767730079)*1000)


#-----------------------------------------------------------------------------------------------------------------------
# B) number of population-wise segregating sites
# z=2 for B

segsitesperid_wl<-fixedsitesperid_asdf_fun(1,2)
segsitesperid_wl1<-segsitesperid_wl[-c(sort(unlist(inf_fsts))),]
segsitesperid_wc<-fixedsitesperid_asdf_fun(2,2)
segsitesperid_wc1<-segsitesperid_wc[-c(sort(unlist(inf_fsts))),]
segsitesperid_el<-fixedsitesperid_asdf_fun(3,2)
segsitesperid_el1<-segsitesperid_el[-c(sort(unlist(inf_fsts))),]
segsitesperid_em<-fixedsitesperid_asdf_fun(4,2)
segsitesperid_em1<-segsitesperid_em[-c(sort(unlist(inf_fsts))),]

pop_segsites<-rbind(sum(segsitesperid_wl1/767730079)*1000,
sum(segsitesperid_wc1/767730079)*1000,
sum(segsitesperid_el1/767730079)*1000,
sum(segsitesperid_em1/767730079)*1000)


#-----------------------------------------------------------------------------------------------------------------------
# 2.2) calc fixed sites per individual

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



all_fixedperid<-rbind(
(cbind(mean((colSums(fixedperid_wl1)/767730079))*1000,sd((colSums(fixedperid_wl1)/767730079))*1000)),
(cbind(sum(fixedperid_wc1/767730079)*1000,0)),
(cbind(mean((colSums(fixedperid_el1)/767730079))*1000,sd((colSums(fixedperid_el1)/767730079))*1000)),
(cbind(mean((colSums(fixedperid_em1)/767730079))*1000,sd((colSums(fixedperid_em1)/767730079))*1000)))


#-----------------------------------------------------------------------------------------------------------------------

# order of parameters

#het_mu1 # mu het - per kbp & normalised by data coverage of empirical data
#het_sd1 # sd het
#pop_fixedsites # population fixed sites  
#pop_segsites # population segregating sites
#all_fixedperid # fixed sites per individual
#fst.out, # fst
#fout.pi, # pi
#fout.taj # tajima's d


# add in het_sd_WC as 0 here - so that same number of rows as simulated data & can filter after   
target=c(het_mu1,
     het_sd1[1],0,het_sd1[2:3],
pop_fixedsites,
pop_segsites,
all_fixedperid,
sapply(fst_df1, mean),
sapply(pi_df1, mean),
sapply(pi_df1, sd),
sapply(tajima_df1, mean),
sapply(tajima_df1, sd))


names(target)<-c(
"het_mu_WL", "het_mu_WC", "het_mu_EL", "het_mu_EM",
"het_sd_WL", "het_sd_WC","het_sd_EL", "het_sd_EM",
"fixedsites_WL", "fixedsites_WC", "fixedsites_EL", "fixedsites_EM",
"segsites_WL", "segsites_WC", "segsites_EL", "segsites_EM",
"fixedperid_mu_WL", "fixedperid_mu_WC", "fixedperid_mu_EL", "fixedperid_mu_EM",
"fixedperid_sd_WL", "fixedperid_sd_WC", "fixedperid_sd_EL", "fixedperid_sd_EM",
"fst_WL.WC","fst_WL.EL","fst_WL.EM","fst_WC.EL","fst_WC.EM","fst_EL.EM",
"pi_mu_WL", "pi_mu_WC", "pi_mu_EL", "pi_mu_EM",
"pi_sd_WL", "pi_sd_WC", "pi_sd_EL", "pi_sd_EM",
 "tajima_mu_WL", "tajima_mu_EL", "tajima_mu_EM",
  "tajima_sd_WL", "tajima_sd_EL", "tajima_sd_EM")


#-----------------------------------------------------------------------------------------------------------------------

colnames(fst_df1)<-c("fst_WL.WC","fst_WL.EL","fst_WL.EM","fst_WC.EL","fst_WC.EM","fst_EL.EM")
colnames(pi_df1)<-c("pi_WL", "pi_WC", "pi_EL", "pi_EM")
colnames(tajima_df1)<-c("tajima_WL", "tajima_EL", "tajima_EM")

pi_df2<-pi_df1[,-c(2)] # remove pi_WC


test<-cbind(fst_df1, pi_df2, tajima_df1)
cormat <- round(cor(test),2)

#-----------------------------------------------------------------------------------------------------------------------

# normalise the heterozygosities

hetperid_wl2<-((rowSums(hetperid_wl1)/767730079))*1000
hetperid_wc2<-(((hetperid_wc1)/767730079))*1000
hetperid_el2<-((rowSums(hetperid_el1)/767730079))*1000
hetperid_em2<-((rowSums(hetperid_em1)/767730079))*1000



fixedsitesperid_wl2<-(((fixedsitesperid_wl1)/767730079))*1000
fixedsitesperid_wc2<-(((fixedsitesperid_wc1)/767730079))*1000
fixedsitesperid_el2<-(((fixedsitesperid_el1)/767730079))*1000
fixedsitesperid_em2<-(((fixedsitesperid_em1)/767730079))*1000


segsitesperid_wl2<-(((segsitesperid_wl1)/767730079))*1000
segsitesperid_wc2<-(((segsitesperid_wc1)/767730079))*1000
segsitesperid_el2<-(((segsitesperid_el1)/767730079))*1000
segsitesperid_em2<-(((segsitesperid_em1)/767730079))*1000


fixedperid_wl2<-((rowSums(fixedperid_wl1)/767730079))*1000
fixedperid_wc2<-(((fixedperid_wc1)/767730079))*1000
fixedperid_el2<-((rowSums(fixedperid_el1)/767730079))*1000
fixedperid_em2<-((rowSums(fixedperid_em1)/767730079))*1000


test1<-cbind(hetperid_wl2,hetperid_wc2,hetperid_el2,hetperid_em2,
fixedsitesperid_wl2,fixedsitesperid_wc2,fixedsitesperid_el2,fixedsitesperid_em2,
segsitesperid_wl2,segsitesperid_wc2,segsitesperid_el2,segsitesperid_em2,
fixedperid_wl2,fixedperid_wc2,fixedperid_el2,fixedperid_em2
    )

colnames(test1)<-c("het_wl","het_wc","het_el","het_em","fixedsites_wl","fixedsites_wc","fixedsites_el","fixedsites_em", "segsites_wl","segsites_wc","segsites_el","segsites_em","fixedperid_wl","fixedperid_wc","fixedperid_el","fixedperid_em")

test2<-cbind(test1,test)

cormat2 <- round(cor(test2),2)

#-----------------------------------------------------------------------------------------------------------------------

# plot correlations between the summary statistics 

# change these to be consistent - ie all lower case
rownames(cormat2)[17:28]<-tolower(rownames(cormat2)[17:28])
colnames(cormat2)[17:28]<-tolower(colnames(cormat2)[17:28])


library(reshape2)
library(ggplot2)

melted_cormat2 <- melt(cormat2)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(cormat2)
upper_tri

# Melt the correlation matrix

melted_cormat2 <- melt(upper_tri, na.rm = TRUE)
# Heatmap

ggheatmap2 <- ggplot(data = melted_cormat2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_classic() + 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 10, hjust = 1))+
 coord_fixed()

pdf("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/corrl.emp.summarystats.unordered.24jan22.replot.pdf") 
print(ggheatmap2)
dev.off()


