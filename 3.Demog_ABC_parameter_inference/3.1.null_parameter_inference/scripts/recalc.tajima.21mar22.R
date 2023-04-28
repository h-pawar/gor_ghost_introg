#!/usr/bin/r
# Mon 21 Mar 2022 16:56:16 CET
# found a bug in how I was calculating the tajimas d statistic =>  recalculate this statistic for all the ~35,700 null demography simulations.
#------------------------------------------------------------------------------------------------------------------------
#module load gcc/6.3.0 R/3.4.2
# for troubleshooting
#mkdir /dev/shm/mydata
#cd /dev/shm/mydata
#------------------------------------------------------------------------------------------------------------------------
require(data.table)
options(stringsAsFactors=F)
library(mgcv)
library(gap)
library(parallel)
library('pegas')
# ptype = array number
ptype=(commandArgs(TRUE))
ptype=as.numeric(as.character(ptype))
print(ptype)
options(scipen=100)
iter=51

# random numbers
randtyp<-sample(100000:500000,10000,replace=F)
#------------------------------------------------------------------------------------------------------------------------

mura<-1.235e-08 # mu, mutation rate
rbm<-9.40e-09 # R recombination rate
gent<-19 # generation time
slen<-40000 #length of region to simulate 40kb
# simulate 40kb b/c calculated S* on a genome-wide scale with a window size of 40kb

mr<-mura*slen*4000 ## scaled mutation rate
## scaled mutation rate. 4 * N * u * length 
  # reference N = 1000
  # => mu * length of fragment * 4 * 1000
sdv=mr*0.233
recv<-rbm*4000*slen
# scaled recombination rate. R * 4 * N * length of fragment

tetaro<-cbind(rnorm(20000,mr,sdv),rnbinom(20000,mu=recv,size=0.5))

tetaro[,1]<-ifelse(tetaro[,1]<0.001,0.001,tetaro[,1])
tetaro[,2]<-ifelse(tetaro[,2]<0.001,0.001,tetaro[,2])

types <- data.table(
  out_pop = c(1,1,3,3)-1,
  in_pop = c(3,4,1,2)-1,
  outgroup = c("WL","WL","EL","EL"),
  ingroup = c("EL","EM","WL","WC"),
  out_chr = c(54,54,18,18),
  in_chr = c(18,24,54,2)
)

types$sum <- types$out_chr + types$in_chr

comps<-list(c(54,0,18,0),c(54,0,0,24),c(54,0,18,0),c(0,2,18,0))

# this number of haploid chr to simulate: w_lowl, w_cros, e_lowl, e_moun
ivec<-c(54,2,18,24);isu=sum(ivec);ivec=paste(ivec,collapse=" ")

#------------------------------------------------------------------------------------------------------------------------
# parameter vals read in from (used to generate the  35543  null simns, for which I have calculated the prev stats & processed in /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/generateabc.6oct.R)

load(file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/null.simns_param"),verbose=T)
inputsets<-simns_param

#------------------------------------------------------------------------------------------------------------------------

config<-c(54,2,18,24)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
#sum(config)
len=40000


#------------------------------------------------------------------------------------------------------------------------

tajima_function<-function(nput) {
msout <-read.ms.output(paste("/dev/shm/mydata/abc.sim_",nput,sep=""))
# tajima's d
# need to recalculate for all null simns -
refun<-function(input) { op=gsub(0,"C",input);gsub(1,"T",op) }
vec1<-seq(1,ncol(msout$gametes[[1]]))
  haploid.loop<-function(y){
    get_gt2<-function(i) {     as.data.frame((msout$gametes[[y]][,c(i)])) }
    lapply(vec1,get_gt2)
  }
  haploids=lapply(1:length(msout$gametes),haploid.loop)
format.function2<-function(y, i){
lapply(haploids[[y]][[i]],refun)
}

loop.taj<-function(x,y){
reps_pop1=list()
for (j in 1:length(msout$gametes)){
pop1=list()
for (i in x: y) {
pop1[[i]]<-unlist(format.function2(j,i))
}
reps_pop1[[j]]<-pop1
}
# add tajima's d to the function
hold.tajima=list()
for (j in 1:length(msout$gametes)){
reps_pop1[[j]]<-reps_pop1[[j]][lengths(reps_pop1[[j]]) != 0]
hold.tajima[[j]]<-tajima.test(as.DNAbin(reps_pop1[[j]]))
}

onlytajima=list()
for (j in 1:length(hold.tajima) ) {
onlytajima[[j]]<-hold.tajima[[j]][1]
}

perwind<-unlist(onlytajima)
out.tajima<-cbind(mean(perwind,na.rm=T), sd(perwind,na.rm=T))
# output both ways
perwind[is.na(perwind)] = 0 
out.tajima1<-cbind(mean(perwind), sd(perwind))

return(list(out.tajima,out.tajima1))

}

# per pop (3: WL, EL, EM), output: mean tajimas d, sd tajimas d, per window tajimas d

allout.taj<-list(
loop.taj(1,54),
loop.taj(57,74),
loop.taj(75,98)
)

## output
return(allout.taj)
}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# give the fixed demog parameter vals
t1=0.0003947368
t2=0.095
t3=0.12
t3.1=t3+0.0025 # 10 generations
t6=0.4473684
t6.1=t6+0.00025 # 1 generation

#-----------------------------------------------------------------------------------------------------------------------

wholefunction=function(input) {
initial_wl=input[1]
initial_wc=input[2]
initial_el=input[3]
initial_em=input[4]
time1=t1
time1_el=input[5]
time2=t2
time2_el=input[6]
time3=t3
time3_em=input[7]
time3.1=time3+0.0025
time3.1_em=input[8]
time4=input[9]
time4_ance=input[10]
time5=input[11]
time5_wl=input[12]
time6=t6
time6_geneflow_w_e=input[13]
time6_geneflow_e_w=input[14]
time6.1=time6+0.00025
time7=input[15]
time7_ancw=input[16]
time8=input[17]
time8_ancew=input[18]

#------------------------------------------------------------------------------------------------------------------------
# generate simulation values, i.e. a list of complete parameters to input to the simulation
simval<-list(
c("-en",time1,paste(3,time1_el,sep=" ")),
c("-en",time2,paste(3,time2_el,sep=" ")),
c("-en",time3,paste(4,time3_em,sep=" ")),
c("-en",time3.1,paste(4,time3.1_em,sep=" ")),
c("-ej",time4,paste(4,3,sep=" ")),
c("-en",time4,paste(3,time4_ance,sep=" ")),
c("-en",time5,paste(1,time5_wl,sep=" ")),
c("-em",time6,paste(3,1,time6_geneflow_w_e,sep=" ")),
c("-em",time6,paste(1,3,time6_geneflow_e_w,sep=" ")),
c("-em",time6.1,paste(3,1,0,sep=" ")),
c("-em",time6.1,paste(1,3,0,sep=" ")),
c("-ej",time7,paste(2,1,sep=" ")),
c("-en",time7,paste(1,time7_ancw,sep=" ")),
c("-ej",time8,paste(3,1,sep=" ")),
c("-en",time8,paste(1,time8_ancew,sep=" "))
)

# get an ordered list and merged part of the simulation parameters - ordering the events is necessary
simval<-do.call(rbind,simval);simval<-simval[order(as.numeric(simval[,2])),];simval=paste(t(simval),collapse=" ")
 
# generate the ms simulation 
# store data temporally - otherwise you flood the cluster - need to create tmpdir
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 2500 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 4 54 2 18 24 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em,"  ",simval," > /dev/shm/mydata/abc.sim_",input[19],sep=""),intern=F)

#------------------------------------------------------------------------------------------------------------------------

tajima_function(input[19])

}



tajresults<-mclapply(inputsets[[ptype]],wholefunction,mc.cores=6,mc.silent=F)

# save in sep dir
rep_param<-inputsets[[ptype]]

save(tajresults,rep_param,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/tajima_21mar22/taj_sim",ptype,sep=""))

for (i in 1:length(inputsets[[ptype]])){
system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[ptype]][[i]][[19]],sep=""),intern=F)
}
