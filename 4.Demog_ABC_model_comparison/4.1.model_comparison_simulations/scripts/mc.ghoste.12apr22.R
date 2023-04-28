#!/usr/bin/r

# Tue 12 Apr 2022 17:03:05 CEST
# generate model comparison simulations using final weighted median posteriors from the final abc parameter inference for ghost gene flow into e_anc

#module load gcc/6.3.0 R/3.4.2
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

iter=50 # to target 50 * 100 -> 5000 reps

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

# sample mutation rate from normal dist (rnorm(n, mean = 0, sd = 1))
# sample recombination rate from -ve binomial dist (sampling parameter vals for input into ms from prob dist using tbs flag)

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

# this number of haploid chr to simulate: w_lowl (1), w_cros (2), e_lowl (3), e_moun (4), ghost (5)
ivec<-c(54,2,18,24,0);isu=sum(ivec);ivec=paste(ivec,collapse=" ")

#------------------------------------------------------------------------------------------------------------------------

# posteriors from generateghosteabc.12apr22.R

# 1) current Ne 
w_lowl_t0=64.7488 
w_cros_t0=14.5582 
e_lowl_t0=20.8943 
e_moun_t0=2.1578
ghost_t0=25
# 2) rapid & severe population decline of E lowland van der valk 2018, plumptr 2016
t1=0.0003947368
e_lowl_t1=0.2431
# 3)  Grauer’s gorillas went through a period of population growth and range expansion  van der valk 2018 & Tocheri 2016
t2=0.095
e_lowl_t2=21.9824
## 4) severe population bottleneck experienced by mountain gorillas
t3=0.12
e_moun_t3=0.1150
t3.1=t3+0.0025 # 10 generations
e_moun_t3.1=3.0090
## 5) E subspecies split # merge eas_moun into eas_lowl 
t4=0.1980
e_anc_t4=5.3253
t_archintrog=0.5037
archaicintrog=99.0813
t_archintrog.1=t_archintrog+0.00025 # 1 generation
## 6) recent contraction in population  w lowland - mcmanus et al. 2015
t5=0.5381 
w_lowl_t5=48.2880 
# 7) migration pulse W > E 
t6=0.4473684
admix_w_e_t6=11.0722
admix_e_w_t6=33.0968 
t6.1=t6+0.00025 # 1 generation
# 8) W subspecies split ie merge wes_cros into wes_lowl
t7=5.9778 
w_anc_t7=98.1352 
# 9) gor split (W - E) # merge eastern into western 
t8=12.7037
gor_anc=14.3641 
t9=45.0179
gor_ghost_anc=35.3814

#------------------------------------------------------------------------------------------------------------------------

inputsets=as.data.frame(t(randtyp[1:iter]))

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


config<-c(54,2,18,24)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
len=40000

all_function<-function(nput) {
  msout <-read.ms.output(paste("/dev/shm/mydata/mc.ghoste.sim_",nput,sep=""))

  # MK: get genotype tables  
  vec<-seq(1,(ncol(msout$gametes[[1]])-1),2)
  ids.loop<-function(y){
    get_gt<-function(i) {     cbind(as.data.frame(paste0(msout$gametes[[y]][,c(i)], msout$gametes[[y]][,c(i+1)]))) }
    lapply(vec,get_gt)
  }
  ids=lapply(1:length(msout$gametes),ids.loop)

  ## MK: output Heterozygosity
  count_fun<-function(y,ids) {
    count_fun2<-function(i,y,ids) { length(which(ids[[y]][[i]]=="01")) + length(which(ids[[y]][[i]]=="10")) }
    unlist(lapply(1:length(vec),count_fun2,ids=ids,y=y))
    }  

  full=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun,ids=ids)
  full_het=rowSums(do.call(cbind,full))/(msout[grep("nreps",names(msout))][[1]]*len)
  mean_fun<-function(typ,fhets) { c(mean(fhets[(typ[1]+1):typ[2]]),sd(fhets[(typ[1]+1):typ[2]])) }
  het_op<-apply(cofig,1,mean_fun,fhets=full_het)*1000
  
  # get heterozygosity per kbp, which is a more accesible value
  

#------------------------------------------------------------------------------------------------------------------------
cofig2<-cbind(c(0,cumsum(config)[-length(config)]),cumsum(config))

#------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# Tue 14 Dec 2021 12:14:33 CET
# population-wise fixed sites & segregating sites, individual fixed sites (like heterozygosity). 

alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))

splitgt<-function(i){
y<-data.frame(matrix(unlist((msout$gametes[[i]])), ncol= sum(config)))

y[,] <- sapply(y[,],as.numeric)
fxd<-which(rowSums(y)<sum(config))
y<-y[fxd,]
return(y)
}

segfun<-function(y,ids) {
      idx<-splitgt(y)
      pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]);op<-table(rowSums(idx[,pop,drop=F]))
      op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)
      op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
      c(op1,op2)}
      apply(cofig2,1,pofun,idx=idx)
  }
  full_segs=do.call(rbind,lapply(1:msout[grep("nreps",names(msout))][[1]],segfun,ids=ids))
 

full_segsites=rbind(colSums(full_segs[seq(1,499,2),]),colSums(full_segs[seq(2,500,2),])) 

#-----------------------------------------------------------------------------------------------------------------------

## individual fixed sites
  count_fun_f<-function(y,ids) {
    count_fun2_f<-function(i,y,ids) { 

y1<-data.frame(matrix(unlist(ids[[y]]), ncol= sum(config)/2))
y0<-data.frame(matrix(unlist((msout$gametes[[y]])), ncol= sum(config)))
y0[,] <- sapply(y0[,],as.numeric)
fxd<-which(rowSums(y0)<sum(config))

length(which(y1[fxd,i]=="11"))
}
    unlist(lapply(1:length(vec),count_fun2_f,ids=ids,y=y))
  } 


  full_f=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun_f,ids=ids)
  full_fix=rowSums(do.call(cbind,full_f))/(msout[grep("nreps",names(msout))][[1]]*len)
  mean_fun_f<-function(typ,ffix) { c(mean(ffix[(typ[1]+1):typ[2]]),sd(ffix[(typ[1]+1):typ[2]])) }
  fix_op<-apply(cofig,1,mean_fun_f,ffix=full_fix)*1000


#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
 
  ##  fst here
  alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))
algam<-do.call(rbind,msout$gametes)
alfrfun<-function(nput) {
   inds<-which(alinds==nput)
   op<-rowSums(algam[,inds])/(length(inds))
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
fst.all(1,2),
fst.all(1,3),
fst.all(1,4),
fst.all(2,3),
fst.all(2,4),
fst.all(3,4)
)
 
#----------------------------------------------------------------------------------------
# pi & tajima's d
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

loop.pi<-function(x,y){

reps_pop1=list()
for (j in 1:length(msout$gametes)){
pop1=list()
for (i in x: y) {
pop1[[i]]<-unlist(format.function2(j,i))
}
reps_pop1[[j]]<-pop1
}

hold.pi=list()
for (j in 1:length(msout$gametes)){
reps_pop1[[j]]<-reps_pop1[[j]][lengths(reps_pop1[[j]]) != 0]
hold.pi[[j]]<-nuc.div(as.DNAbin(reps_pop1[[j]]))
}

out.pi<-cbind(mean(unlist(hold.pi)), sd(unlist(hold.pi)))
return(out.pi)
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
# process tajimas d below
hold.tajima=list()
for (j in 1:length(msout$gametes)){
reps_pop1[[j]]<-reps_pop1[[j]][lengths(reps_pop1[[j]]) != 0]
hold.tajima[[j]]<-tajima.test(as.DNAbin(reps_pop1[[j]]))
}

out.tajima<-cbind(mean(unlist(hold.tajima)), sd(unlist(hold.tajima)))
return(out.tajima)
}

fout.pi<-rbind(
loop.pi(1,54),
loop.pi(55,56),
loop.pi(57,74),
loop.pi(75,98)
)


#----------------------------------------------------------------------------------------
# amended Thu 24 Mar 2022 15:28:41 CET


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

perwind[is.na(perwind)] = 0 
out.tajima1<-cbind(mean(perwind), sd(perwind))

return(out.tajima1)
}

# per pop (3: WL, EL, EM), output: mean tajimas d, sd tajimas d, per window tajimas d

allout.taj<-list(
loop.taj(1,54),
loop.taj(57,74),
loop.taj(75,98)
)


#------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

  return(list(het_op,full_segsites,fix_op,fst.out,fout.pi,allout.taj))


}

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

wholefunction=function(input) {

initial_wl=w_lowl_t0
initial_wc=w_cros_t0
initial_el=e_lowl_t0
initial_em=e_moun_t0
initial_ghost=ghost_t0
time1=t1
time1_el=e_lowl_t1
time2=t2
time2_el=e_lowl_t2
time3=t3
time3_em=e_moun_t3
time3.1=time3+0.0025
time3.1_em=e_moun_t3.1
time4=t4
time4_ance=e_anc_t4
timeI=t_archintrog
timeI_archgeneflow=archaicintrog
timeI.1=timeI+0.00025
time5=t5
time5_wl=w_lowl_t5
time6=t6
time6_geneflow_w_e=admix_w_e_t6
time6_geneflow_e_w=admix_e_w_t6
time6.1=time6+0.00025
time7=t7
time7_ancw=w_anc_t7
time8=t8
time8_ancew=gor_anc
time9=t9
time9_allanc=gor_ghost_anc

#------------------------------------------------------------------------------------------------------------------------
# generate simulation values, i.e. a list of complete parameters to input to the simulation
simval<-list(
c("-en",time1,paste(3,time1_el,sep=" ")),
c("-en",time2,paste(3,time2_el,sep=" ")),
c("-en",time3,paste(4,time3_em,sep=" ")),
c("-en",time3.1,paste(4,time3.1_em,sep=" ")),
c("-ej",time4,paste(4,3,sep=" ")),
c("-en",time4,paste(3,time4_ance,sep=" ")),
c("-em",timeI,paste(3,5,timeI_archgeneflow,sep=" ")),
c("-em",timeI.1,paste(3,5,0,sep=" ")),
c("-en",time5,paste(1,time5_wl,sep=" ")),
c("-em",time6,paste(3,1,time6_geneflow_w_e,sep=" ")),
c("-em",time6,paste(1,3,time6_geneflow_e_w,sep=" ")),
c("-em",time6.1,paste(3,1,0,sep=" ")),
c("-em",time6.1,paste(1,3,0,sep=" ")),
c("-ej",time7,paste(2,1,sep=" ")),
c("-en",time7,paste(1,time7_ancw,sep=" ")),
c("-ej",time8,paste(3,1,sep=" ")),
c("-en",time8,paste(1,time8_ancew,sep=" ")),
c("-ej",time9,paste(5,1,sep=" ")),
c("-en",time9,paste(1,time9_allanc,sep=" "))
)

# get an ordered list and merged part of the simulation parameters - ordering the events is necessary
simval<-do.call(rbind,simval);simval<-simval[order(as.numeric(simval[,2])),];simval=paste(t(simval),collapse=" ")

 
# generate the ms simulation 
# store data temporally - otherwise you flood the cluster - need to create tmpdir

system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 250 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 5 54 2 18 24 0 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em," -n 5 ",initial_ghost,"  ",simval," > /dev/shm/mydata/mc.ghoste.sim_",input[1],sep=""),intern=F)


#------------------------------------------------------------------------------------------------------------------------

return(all_function(input[1]))

}


# the following command generates simulated data for the "iter" number of times and interprets the output; it divides this task to 6 cores on the cluster, and collects the data of these "iter" iterations into a single object
simuresults<-mclapply(inputsets,wholefunction,mc.cores=6,mc.silent=F)

save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/final_12apr22/ghoste/test.mc.ghoste_sim",ptype,sep=""))

for (i in 1:iter){
# ABC + ghost
system(paste("rm /dev/shm/mydata/mc.ghoste.sim_",inputsets[[i]][[1]],sep=""),intern=F)
}

#------------------------------------------------------------------------------------------------------------------------
