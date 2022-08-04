#!/usr/bin/env
# Tue 29 Jun 2021 15:12:06 CEST
# calc summary statistics for the empirical data in the same way (in 40kb windows) as for the simulated data in test.abc.model.v4.R
#  heterozygosity, number of segregating sites, pairwise Fst, pi, Tajima's D)
# module load gcc/6.3.0 R/3.4.2 xz/5.2.2  BCFTOOLS/1.6

#-----------------------------------------------------------------------------------------------------------------------
# Tue 17 Aug 2021 17:06:13 CEST
#bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_21.vcf.gz
  # order of individuals in the vcf != order of individuals simulated   
    # need to regenerate stats for the empirical data**

# order in vcf - subspecies - order of pop in simns 
#1:12 - b.b e_moun (4)
#13:21 - b.g e_lowl (3)
#22:22 - g.d w_cros (2)
#23:49 - g.g w_lowl (1)

# => I should output stats in order of : 23:49, 22:22 , 13:21, 1:12 (in order to have the same order of the output as the simulations)

#-----------------------------------------------------------------------------------------------------------------------

# Mon 23 Aug 2021 16:56:58 CEST - 
#I should regenerate the summary stats for the empirical data
# [- outputting sds of pi & tajimas d -> in the final stage]
#- outputting sum of segregating sites per population instead of the mean
#- outputting the sum of heterozygosity per window

#MK:
#yes, it is very important to compare the same thing, 
#so better re-calculate it the same way you did for the simulations. 
#Also, I think you should not include "id" in the ABC, because this may introduce unnecessary noise.
#Lets see if it looks better when including the missing sds and correcting these steps. 
#Depending on what you see, you may also consider using a higher tol value (0.01 for example).

#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
require(data.table)
options(stringsAsFactors=F)
library(ape)
library(pegas)

#chr=22 # for testing, otherwise send from bash wrapper script
chr=args[1]

#emp.data<-system(paste("bcftools query -f '%POS [%GT ]\n' /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chr,".vcf.gz",sep=""),intern=T)
# # read in alleles directly from bcftools : %TGT = Translated genotype (e.g. C/A)
emp.data<-system(paste("bcftools query -f '%POS [%GT %TGT ]\n' /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chr,".vcf.gz",sep=""),intern=T)
emp.data<-do.call(rbind,strsplit(emp.data,split=" ")) ## ie list merged vertically into df
pos<-emp.data[,1] 
pos<-as.numeric(pos)
window<-seq(min(pos), max(pos), 40000) # splits chromosome into 40kb non-overlapping regions
emp.data1<-cbind(pos, emp.data[,2:ncol(emp.data)])

# split into 2 dfs one where GTs as 0|1, the other with tgts as C|A
even_indexes<-seq(2,ncol(emp.data1),2)
odd_indexes<-seq(1,ncol(emp.data1),2)
emp.data.gts<-emp.data1[,c(1,even_indexes)]
emp.data.tgts<-emp.data1[,odd_indexes]

window<-seq(min(pos), max(pos), 40000) # splits chromosome into 40kb non-overlapping regions

#intvl<-c(window[i],window[i+1]) # loop through the intervals to extract GTs for each 40kb region

# for the first intvl
#intvl<-c(window[1],window[1+1])
#intvl
#[1] 16211755 16251755

# now extract GTs for this 40kb region
#win.intvl<-emp.data1[intvl[1] <= pos & pos < intvl[2] ,]
#win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]

#config<-c(54,2,18,24)
config<-c(24,18,2,54)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
sum(config)
len=40000

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# function to calculate summary statistics per window:

stats.in.win_function<-function(i) {
intvl<-c(window[[i]],window[[i+1]]) # interval to extract GTs from (40kb region)
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 1)
# count heterozygous GTs per individual & divide by window length

# not just 0s & 1s here, also alternate alleles - 
het_fun<-function(id){ 
	(length(which(win.intvl[,id]=="0|1")) + 
	length(which(win.intvl[,id]=="1|2")) +
	length(which(win.intvl[,id]=="2|1")) +
	length(which(win.intvl[,id]=="0|2")) +
	length(which(win.intvl[,id]=="3|2")))/len }

out.winhet<-lapply(2:ncol(win.intvl), het_fun) # apply for all 49 individuals, from col 2-50
# output this & then take mean, sd over all the segments 

# subsets of the heterozygosities - split per pop (WL, WC, EL, EM) - as nested list (of 4)
# => I should output stats in order of : 23:49, 22:22 , 13:21, 1:12
                  # WL , WC, EL, EM

#  output the sum of heterozygosity per window
# & take the mean & sds only at the end when I have all the informative windows
# MK: that is right, it should be per population, not per window. 

out.winhet.perpop<-cbind(sum(unlist(out.winhet)[23:49]),
  sum(unlist(out.winhet)[22]),
  sum(unlist(out.winhet)[13:21]),
  sum(unlist(out.winhet)[1:12]))


# calc all stats per window then loop over windows 

#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 2)
# number of segregating sites
# seg sites - can't just do length(which(nput!="00")), b/c also 1|1, 2|2 in vcf
#> unique(emp.data1[,2])
#[1] "0|0" "1|1" "0|1" "1|2" "2|2" "2|1" "0|2" "3|2"


segs_fun<-function(id){
	(length(which(win.intvl[,id]=="0|1")) + 
	length(which(win.intvl[,id]=="1|2")) +
	length(which(win.intvl[,id]=="2|1")) +
	length(which(win.intvl[,id]=="0|2")) +
	length(which(win.intvl[,id]=="3|2")))
 }

out.winsegs<-lapply(2:ncol(win.intvl), segs_fun)

out.winsegs.perpop<-cbind(sum(unlist(out.winsegs[23:49])),
sum(unlist(out.winsegs[22])),
sum(unlist(out.winsegs[13:21])),
sum(unlist(out.winsegs[1:12]))
    )


#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 3)
# pairwise Fst
# first convert GTs to freqs
splitgts_fun<-function(id){
do.call('rbind', strsplit(as.character(win.intvl[,id]),'|',fixed=TRUE))	
}

win.intvl.hapl<-lapply(2:ncol(win.intvl), splitgts_fun)

y<-data.frame(matrix(unlist(win.intvl.hapl), ncol= sum(config)))
alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))

# convert from character to numeric
y[,] <- sapply(y[,],as.numeric)

alfrfun<-function(nput) {
   inds<-which(alinds==nput)
   op<-rowSums(y[,inds])/(length(inds))
   return(op)
   }
  alfre<-lapply(1:4,alfrfun)
# frequencies needed for input into the fst function: 
bhatiaFST<-function(pop1f,pop2f,nchr1,nchr2) {
  altAllelePop1 <- round(pop1f*nchr1) # count of alt (1) alleles
  altAllelePop2 <- round(pop2f*nchr2)
  
  h1 <- (altAllelePop1*(nchr1-altAllelePop1))/(nchr1*(nchr1-1))
  h2 <- (altAllelePop2*(nchr2-altAllelePop2))/(nchr2*(nchr2-1))
  nHat <- (pop1f - pop2f)^2  - (h1/nchr1) - (h2/nchr2)
  dHat <- nHat+h1+h2
  fst <- nHat/dHat
  return(fst)
}

fst.all<-function(j,i){
fst.hold<-bhatiaFST(alfre[[j]],alfre[[i]], config[j], config[i])
# add this line into the function for calculating pairwise fst, set negative fst vals to 0 (considered loci with no popn differentiation)
fst.hold<-replace(fst.hold, fst.hold<0, 0) 
mean.fst<-mean(fst.hold, na.rm=TRUE) # remove nan vals
return(mean.fst)
}

#fst.out<-rbind(
#fst.all(1,2), # WL - WC # (4,3)
#fst.all(1,3), # WL - EL # (4,2)
#fst.all(1,4), # WL - EM # (4,1)
#fst.all(2,3), # WC - EL # (3,2)
#fst.all(2,4), # WC - EM # (3,1)
#fst.all(3,4) # EL - EM # (2,1)
#)

fst.out<-rbind(
fst.all(4,3),
fst.all(4,2),
fst.all(4,1),
fst.all(3,2),
fst.all(3,1),
fst.all(2,1)
)
#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 4)
# pi
win.intvl.tgts<-emp.data.tgts[intvl[1] <= pos & pos < intvl[2] ,]
# split into haploids
splitgts_fun_tgts<-function(id){
do.call('rbind', strsplit(as.character(win.intvl.tgts[,id]),'|',fixed=TRUE))    
}

win.intvl.hapl.tgts<-lapply(2:ncol(win.intvl.tgts), splitgts_fun_tgts)

pi.out<-rbind(nuc.div(as.DNAbin(win.intvl.hapl.tgts[23:49])),
nuc.div(as.DNAbin(win.intvl.hapl.tgts[22])),
nuc.div(as.DNAbin(win.intvl.hapl.tgts[13:21])),
nuc.div(as.DNAbin(win.intvl.hapl.tgts[1:12])))

# pi for 28:
# [1] NaN # gives NaN here - whereas was not obtaining NaNs for CR pi for the simulated data
# here is perhaps b/c homozygous positions only - may not be major problem

#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 5)
# tajima's d

# here single values rather than sds as well - not looping through simn reps 
tajima.out<-rbind(
tajima.test(as.DNAbin(win.intvl.hapl.tgts[23:49]))[[1]],
tajima.test(as.DNAbin(win.intvl.hapl.tgts[13:21]))[[1]],
tajima.test(as.DNAbin(win.intvl.hapl.tgts[1:12]))[[1]])

#-----------------------------------------------------------------------------------------------------------------------
identifier<-cbind(chr,i, intvl[1],intvl[2]) # output identifier for this window - chr, window number, start pos, end pos

 ## output
 return(list(identifier, out.winhet.perpop,
out.winsegs.perpop,
fst.out,
pi.out,
tajima.out))
}
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# me - after trying random numbers - the function gives output for windows 1& 4, but not for windows 2,3,5 
# - whose intervals do not contain snps in the vcf. 
# So the problem is not the function per se, but that I have created window intervals which do not contain
# positions in the vcf.

# MK 
#good! This is actually an important feature. 
#You should get the fraction of each window that you have data for. 
#Then, you only use windows with >3/4 of positions with information to get the summary stats. 
#Otherwise, you would compare a lot of "half-empty" windows to "complete" simulated windows.

# MK: 
#you want to get the summary statistics for windows in the genome that are informative enough to 
#compare to the simulated data. In simulated data, you have 100% of positions with information, 
#but in the real data that is rarely the case. 
#Thats why you need a cutoff for which windows to use for the comparison. 
#If you take windows where only 20% of the 40k positions have information, 
#the number of segregating sites is very likely underestimated. 
#Probably the best is to calculate the summary stats for all windows in the genome, 
#and then subset the data before you calculate the mean and SD.

# MK (re stats.emp.data.clean.R)
#I think you are calculating the fraction covered by segregating sites, 
#which is expected to be very small. 
#What you would need is the fraction of sites where there is data (including 0/0 across individuals). 
#I have done this for the Skov method in 1kbp windows, I will look up if that can be easily adapted to 40kbp windows

# MK 28/6/21
#Calculate the statistics for all windows where you have data
#then filter by informative windows - from the unfiltered vcf - where will assess how many GTs per window 

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# extract identifier & intervals of non-empty windows
hold.windows=list()
for (j in 1:(length(window)-1)){
intvl<-c(window[[j]],window[[j+1]])
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
if (length(win.intvl)!=0)
hold.windows[[j]]<-cbind(j, intvl[1], intvl[2], nrow(win.intvl)/len)
}

# for windows with only 1 snp - the 4th col is empty - b/c nrow gives null (for vector)
# remove null elements of list ie where  (length(win.intvl)==0)
windows.wdata<-hold.windows[-which(sapply(hold.windows, is.null))]

# ie need to filter out windows with only 1 snp as well
fraction.wind=list()
for (j in 1:length(windows.wdata)){
fraction.wind[[j]]<-cbind(windows.wdata[[j]][1],windows.wdata[[j]][4])
}

df.fraction.wind <- data.frame(matrix(unlist(fraction.wind), nrow=length(fraction.wind), byrow=TRUE),stringsAsFactors=FALSE)

#remove rows with nas - ie windows where there was only 1 snp
df1.fraction.wind<-na.omit(df.fraction.wind)

# ie 19 has been removed from here
#head(df1.fraction.wind)
#  X1       X2
#1  1 0.000075
#2  4 0.000075
#4 22 0.000075
#5 23 0.000775
#6 24 0.002150
#7 25 0.002025

# apply stats.in.win_function to windows with >1snp
	# & writing out the chr, window number & positions as well
stats.out=list()
for (j in 1:nrow(df1.fraction.wind)){
stats.out[[j]]<-stats.in.win_function(df1.fraction.wind$X1[[j]])
}
#There were 11 warnings (use warnings() to see them)

# may need to find more efficient way than looping..

# need to write out summary stats from empirical data
save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/stats_chr",chr,sep=""))



# gives output despite the warnings
#head(stats.out)
#[[1]]
#[[1]][[1]]
#     chr i                  
#[1,]  22 1 16211755 16251755

#[[1]][[2]]
#             [,1] [,2]         [,3]         [,4]
#[1,] 1.851852e-05    0 2.777778e-06 3.958333e-05
#[2,] 2.360797e-05   NA 8.333333e-06 1.982404e-05

#[[1]][[3]]
#             [,1] [,2]         [,3]         [,4]
#[1,] 1.851852e-05    0 2.777778e-06 3.958333e-05

#[[1]][[4]]
#          [,1]
#[1,] 0.1698113
#[2,] 0.1088597
#[3,] 0.0502461
#[4,] 0.0000000
#[5,] 0.2318841
#[6,] 0.1937681

#[[1]][[5]]
#           [,1]
#[1,] 0.16049383
#[2,]        NaN
#[3,] 0.03703704
#[4,] 0.16919192

#[[1]][[6]]
#            [,1]
#[1,]  1.81108200
#[2,] -1.08822734
#[3,]  0.07176555

