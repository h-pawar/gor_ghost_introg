# Mon 27 Sep 2021 10:16:39 CEST
# building from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/het.emp.31aug.R

# new statistics to calculate in the empirical data which were introduced in the new simulation batch:
    #full_segsites is the number of population-wise fixed sites and the number of population-wise segregating sites; 
    # fix_op is the fixed sites per individual (same format as heterozygotes in your script)

#-----------------------------------------------------------------------------------------------------------------------
# NB - new simns generated with
#/scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R # new version supercedes test.abc.model.v4.R (used to generate prev 1-2000 simns)
#/scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.arr # bash wrapper script
#-----------------------------------------------------------------------------------------------------------------------

# module load gcc/6.3.0 R/3.4.2 xz/5.2.2  BCFTOOLS/1.6
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
# Tue 28 Sep 2021 19:24:40 CEST
# but issue - in segs.emp.27sep.R
# looking at 0s,1s but not 2nd,3rd alt alleles ** 
#could deal with this by assigning the 2/3 to 1s - to see the fixed sites per population


etest<-gsub("2", "1", emp.data.gts[,c(2:ncol(emp.data.gts))], fixed=TRUE)
etest<-gsub("3", "1", etest[,c(1:ncol(etest))], fixed=TRUE)

#ncol(etest)
#[1] 49

#ncol(emp.data.gts)
#[1] 50

emp.data.gts<-cbind(emp.data.gts[,1],etest)

# ncol(emp.data.gts)
#[1] 50
#------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# function to calculate summary statistics per window:

het.in.win_function<-function(i) {

# testing i=1
intvl<-c(window[[i]],window[[i+1]]) # interval to extract GTs from (40kb region)
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# Fri 10 Sep 2021 16:54:46 CEST
# MK:
#lets try to improve the simulations.
#I have the following suggestions.

#1. Improvements in the statistics. I found that there seem to be no fully fixed sites in the simulated data. 
#I suggest to calculate the following statistics: 
#population-wise fixed sites & segregating sites, individual fixed sites (like heterozygosity). 
#That could be done this way, building on your code:

## you need to define beforehands that cofig2 is:
cofig2<-cbind(c(0,cumsum(config)[-length(config)]),cumsum(config))

#cofig2
#     [,1] [,2]
#[1,]    0   24
#[2,]   24   42
#[3,]   42   44
#[4,]   44   98


 #it makes sense to take the population-wise segregating sites as non-fixed sites. This seems like a more informative measure. 
 #This would exclude the CRG population, because this value is the same as heterozygosity with just one individua

#-----------------------------------------------------------------------------------------------------------------------


#Tue  5 Oct 2021 18:33:48 CEST
#MK: What I suggest is the following:
#a) For the empirical data, you do not exclude differences to human (1/1 across all individuals). 
#Of course, then the sites that are fixed per individual and per population are inflated, and they appear very high. 
#The "fixedsites" and "fixedsitesperid" are almost the same for each population. You need to adjust this, for example by doing
#  y[,] <- sapply(y[,],as.numeric)
#  fxd<-which(rowSums(y)<sum(config))
#  y<-y[fxd,]

#and then re-calculate the segfun across chromosomes. The same is a problem for fixed sites per population, 
#you need to modify the function:

#fixedid_fun<-function(id){
#    length(which(win.intvl[fxd,id]=="1|1")) }

#As a consequence the fixedsitesperid will go down by an order of magnitude (mas o menos), 
#you need to see how that fits with the empirical data, but it may end somewhere within but at the lower end. 
#For the fixedsites per population, the same is expected. For the segsites, nothing changes.
#-----------------------------------------------------------------------------------------------------------------------


# 1) #population-wise fixed sites and the number of population-wise segregating sites

splitgts_fun<-function(id){
do.call('rbind', strsplit(as.character(win.intvl[,id]),'|',fixed=TRUE)) 
}

win.intvl.hapl<-lapply(2:ncol(win.intvl), splitgts_fun)

y<-data.frame(matrix(unlist(win.intvl.hapl), ncol= sum(config)))
alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))

# convert from character to numeric
y[,] <- sapply(y[,],as.numeric)
fxd<-which(rowSums(y)<sum(config))
y<-y[fxd,]

segfun<-function(y) {
      idx<-y
      pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]);op<-table(rowSums(idx[,pop,drop=F]))
      op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)
      op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
      c(op1,op2)}
      apply(cofig2,1,pofun,idx=idx)
  }

#segfun(y)
#     [,1] [,2] [,3] [,4]
#[1,]    0    0    0    0
#[2,]    2    2    0    3

# this is in the order of : EM, EL, WC, WL

# to match simld data may be best to output as WL,WC,EL,EM (could do this here or at end?)

# output segfun(y)
#eseg<-segfun(y) # EM, EL, WC, WL
eseg<-segfun(y)[,c(4,3,2,1)] # in order of WL,WC,EL,EM
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
## 2) fixed sites per individual (same format as heterozygotes in your script).
# go to next section - this one had errors


# for this need diploids rather than haploids
# gsub("|", "", win.intvl[,c(2:ncol(win.intvl))], fixed=TRUE)

#dipl<- gsub("|", "", win.intvl[,c(2:ncol(win.intvl))], fixed=TRUE)

# dipl[1:24]
# [1] "00" "00" "00" "01" "00" "01" "01" "00" "01" "00" "00" "00" "00" "00" "00"
#[16] "01" "00" "01" "00" "00" "00" "00" "00" "00"

# ie need to do this per popn
#length(which(dipl[1:24]=="11"))

#count_fun_f<-function(y0,y1) {
#o<-length(which(dipl[y0:y1]=="11"))
#return(o)
#}

# cofig 2 doesn't work here b/c now diploids instead of haploids
#cofig[i,1]+1, cofig[i,2]

# count_fun_f(cofig[2,1]+1,cofig[2,2])
#[1] 0

#count_fun_f<-function(i) {
#o<-length(which(dipl[cofig[i,1]+1:cofig[i,2]]=="11"))
#return(o)
#}

#https://stackoverflow.com/questions/37704470/error-in-match-funfun
# sapply(1:4,count_fun_f) # output this
#[1] 0 0 0 0

#eop<-sapply(1:4,count_fun_f) # EM, EL, WC, WL

#eop1<-c(eop[4],eop[3],eop[2],eop[1]) # WL, WC, EL, EM

# ie output this as counts
# then will need to calc over windows
  # & standardise as with heterozygosity
    # ie divide by len (taking into account missingness of data & filtering performed) then *1000

#-----------------------------------------------------------------------------------------------------------------------

# Tue 28 Sep 2021 20:10:45 CEST
# going back to how i was calculating heterozygosity in het.emp.31aug.R

# ncol(win.intvl)
#[1] 50


# this is the equivalent for  fixed sites per individual 

#fixedid_fun<-function(id){ 
#  length(which(win.intvl[,id]=="1|1")) }

#you need to modify the function:
fixedid_fun<-function(id){
    length(which(win.intvl[fxd,id]=="1|1")) }

out.fixedid<-lapply(2:ncol(win.intvl), fixedid_fun)

fixedsitesperid<-list(
unlist(out.fixedid)[23:49],
unlist(out.fixedid)[22],
unlist(out.fixedid)[13:21],
unlist(out.fixedid)[1:12]
  )

# no fixed sites in this window

# fixedsitesperid
#[[1]]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

#[[2]]
#[1] 0

#[[3]]
#[1] 0 0 0 0 0 0 0 0 0

#[[4]]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0

# this is then outputting as counts

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
identifier<-cbind(chr,i, intvl[1],intvl[2]) # output identifier for this window - chr, window number, start pos, end pos

 ## output
 # ie output 1) window identifier
                # 2) number of population-wise fixed sites and the number of population-wise segregating sites; 
                # 3) fixed sites per individual

# return(list(identifier, eseg, eop))
# return(list(identifier, eseg, eop1))
# amended how i was calculating the stats:
return(list(identifier, eseg, fixedsitesperid))
}
#------------------------------------------------------------------------------------------------------------------------

# see if this works for chr22

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

# apply het.in.win_function to windows with >1snp
  # & writing out the chr, window number & positions as well
stats.out=list()
for (j in 1:nrow(df1.fraction.wind)){
stats.out[[j]]<-het.in.win_function(df1.fraction.wind$X1[[j]])
}

#-----------------------------------------------------------------------------------------------------------------------
# head(stats.out)
#[[1]]
#[[1]][[1]]
#     chr i                  
#[1,]  22 1 16211755 16251755

#[[1]][[2]]
#     [,1] [,2] [,3] [,4]
#[1,]    0    0    0    0
#[2,]    3    0    2    2

#[[1]][[3]]
#[[1]][[3]][[1]]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

#[[1]][[3]][[2]]
#[1] 0

#[[1]][[3]][[3]]
#[1] 0 0 0 0 0 0 0 0 0

#[[1]][[3]][[4]]
# [1] 0 0 0 0 0 0 0 0 0 0 0 0
#-----------------------------------------------------------------------------------------------------------------------


# may need to find more efficient way than looping..

# write out summary stats from empirical data

# mkdir to store this & change name of output files to segs
# mkdir -p /scratch/devel/hpawar/admix/abc/emp.data/test/segs

#save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr",chr,sep=""))
save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/seg_chr",chr,sep=""))

