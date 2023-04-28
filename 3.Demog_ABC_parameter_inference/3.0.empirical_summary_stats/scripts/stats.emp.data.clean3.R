#!/usr/bin/env
# Tue 29 Jun 2021 15:12:06 CEST
# calc summary statistics for the empirical data in the same way (in 40kb windows) as for the simulated data in test.abc.model.v4.R
#  heterozygosity, number of segregating sites, pairwise Fst, pi, Tajima's D)
# module load gcc/6.3.0 R/3.4.2 xz/5.2.2  BCFTOOLS/1.6

#-----------------------------------------------------------------------------------------------------------------------
# Tue 17 Aug 2021 17:06:13 CEST
#bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_21.vcf.gz  # order of individuals in the vcf != order of individuals simulated   

# order in vcf - subspecies - order of pop in simns 
#1:12 - b.b e_moun (4)
#13:21 - b.g e_lowl (3)
#22:22 - g.d w_cros (2)
#23:49 - g.g w_lowl (1)

# => I should output stats in order of : 23:49, 22:22 , 13:21, 1:12 (in order to have the same order of the output as the simulations)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
require(data.table)
options(stringsAsFactors=F)
library(ape)
library(pegas)
chr=args[1]

# read in alleles directly from bcftools : %TGT = Translated genotype (e.g. C/A)
emp.data<-system(paste("bcftools query -f '%POS [%GT %TGT ]\n' /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chr,".vcf.gz",sep=""),intern=T)
emp.data<-do.call(rbind,strsplit(emp.data,split=" ")) 
pos<-emp.data[,1] 
pos<-as.numeric(pos)
window<-seq(min(pos), max(pos), 40000) # splits chromosome into 40kb non-overlapping regions
emp.data1<-cbind(pos, emp.data[,2:ncol(emp.data)])

# split into 2 dfs one where GTs as 0|1, the other with tgts as C|A
even_indexes<-seq(2,ncol(emp.data1),2)
odd_indexes<-seq(1,ncol(emp.data1),2)
emp.data.gts<-emp.data1[,c(1,even_indexes)]
emp.data.tgts<-emp.data1[,odd_indexes]

window<-seq(min(pos), max(pos), 40000) 

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

het_fun<-function(id){ 
	(length(which(win.intvl[,id]=="0|1")) + 
	length(which(win.intvl[,id]=="1|2")) +
	length(which(win.intvl[,id]=="2|1")) +
	length(which(win.intvl[,id]=="0|2")) +
	length(which(win.intvl[,id]=="3|2")))/len }

out.winhet<-lapply(2:ncol(win.intvl), het_fun) 
out.winhet.perpop<-cbind(sum(unlist(out.winhet)[23:49]),
  sum(unlist(out.winhet)[22]),
  sum(unlist(out.winhet)[13:21]),
  sum(unlist(out.winhet)[1:12]))


# calc all stats per window then loop over windows 

#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 2)
# number of segregating sites

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
# set negative fst vals to 0 (considered loci with no popn differentiation)
fst.hold<-replace(fst.hold, fst.hold<0, 0) 
mean.fst<-mean(fst.hold, na.rm=TRUE) # remove nan vals
return(mean.fst)
}

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

#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 5)
# tajima's d

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

# extract identifier & intervals of non-empty windows
hold.windows=list()
for (j in 1:(length(window)-1)){
intvl<-c(window[[j]],window[[j+1]])
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
if (length(win.intvl)!=0)
hold.windows[[j]]<-cbind(j, intvl[1], intvl[2], nrow(win.intvl)/len)
}

windows.wdata<-hold.windows[-which(sapply(hold.windows, is.null))]

fraction.wind=list()
for (j in 1:length(windows.wdata)){
fraction.wind[[j]]<-cbind(windows.wdata[[j]][1],windows.wdata[[j]][4])
}

df.fraction.wind <- data.frame(matrix(unlist(fraction.wind), nrow=length(fraction.wind), byrow=TRUE),stringsAsFactors=FALSE)

#remove rows with nas - ie windows where there was only 1 snp
df1.fraction.wind<-na.omit(df.fraction.wind)

stats.out=list()
for (j in 1:nrow(df1.fraction.wind)){
stats.out[[j]]<-stats.in.win_function(df1.fraction.wind$X1[[j]])
}

save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/stats_chr",chr,sep=""))


