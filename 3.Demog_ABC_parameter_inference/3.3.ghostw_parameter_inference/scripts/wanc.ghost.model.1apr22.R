#!/usr/bin/r

# Fri  1 Apr 2022 11:51:19 CEST
# ghost gene flow into W model - generate ABC parameter inference simulations 
#-----------------------------------------------------------------------------------------------------------------------

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
#ptype=1 # for troubleshooting purposes
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

#--------------------------------------------------------------------------------------------------------------------
# annotate below with 1) weighted median posteriors 2) visually what makes sense from the posterior curves 

# 1) priors for current Ne 
w_lowl_t0=64.7488 # w_lowl_t0 # could fix?
w_cros_t0=14.5582 # w_cros_t0 # fix
e_lowl_t0=20.8943 # e_lowl_t0 # fix
e_moun_t0=2.1578 # e_moun_t0 # fix
ghost_t0=runif(1000,min=0.1,max=100)
# 2) rapid & severe population decline of E lowland over last 20y (0.2kya) of ~80% - van der valk 2018, plumptr 2016
t1=0.0003947368
e_lowl_t1=0.2431 # e_lowl_t1 # fix
# 3)  Grauer’s gorillas went through a period of population growth and range expansion 5,000–10,000 years ago [22, 23]' van der valk 2019 citing van der valk 2018 & Tocheri 2016
t2=0.095
e_lowl_t2=21.9824 # e_lowl_t2 # fix
## 4) severe population bottleneck experienced by mountain gorillas - xue et al - this was likely after their split from eastern lowland
t3=0.12
e_moun_t3=0.1150 # e_moun_t3 # fix
# MK: 148 You may model this bottleneck as a temporal change, where at t3 it goes down to this value, and t3 minus 5 or 10 generations it goes back to the previous/terminal Ne.
t3.1=t3+0.0025 # 10 generations
e_moun_t3.1=runif(1000,min=2.9191,max=19.1638)
#e_moun_t3.1=9.9114 # e_moun_t3.1 - priors - take credible interval # 2.9191, 19.1638
## 5) E subspecies split # merge eas_moun into eas_lowl ** (adding parameters from xue et al. 2015 psmc into the mcmanus model)
t4=runif(1000,min=0.1446,max=0.2420)
#t4=0.1839 # t4  - priors t4 # 0.1446, 0.2420 
e_anc_t4=5.3253 # e_anc_t4 # fix 
# t7=5.9778  WL split time, = 5.9778  + 0.0025 = 5.9803
t_archintrog=runif(1000,min=5.9803,max=14.3616)
# if take the lower bound woudl be min=0.1471
archaicintrog=runif(1000,min=0,max=100)
t_archintrog.1=t_archintrog+0.00025 # 1 generation
## 6) recent contraction in population size 23 ka in w lowland - mcmanus et al. 2015
t5=0.5381 # t5 # fix
w_lowl_t5=48.2880 # w_lowl_t5 # fix
# 7) migration pulse W > E 3% - for admixture need to give as a proportion not percent** 
t6=0.4473684
admix_w_e_t6=runif(1000,min=6.7955,max=99.2159)
#admix_w_e_t6=61.0343 # admix_w_e_t6 - priors #6.7955, 99.2159
admix_e_w_t6=33.0968 # admix_e_w_t6 # fix
t6.1=t6+0.00025 # 1 generation
# 8) W subspecies split ie merge wes_cros into wes_lowl
t7=5.9778 # t7 # fix
w_anc_t7=98.1352 # w_anc_t7 # fix
# 9) gor split (W - E) # merge eastern into western at time 261,000/(4 * 1000 * 19) = 3.434211
t8=12.7037 # t8 # fix
gor_anc=14.3641 # gor_anc # fix
# 11) split between common anc of all gorillas + ghost pop
t9=runif(1000,min=15,max=50)
gor_ghost_anc=runif(1000,min=10,max=100)

#------------------------------------------------------------------------------------------------------------------------

inputsets=
  as.data.frame(t(cbind(
  sample(ghost_t0,iter,replace=T),
  sample(e_moun_t3.1,iter,replace=T),
  sample(t4,iter,replace=T),
  sample(t_archintrog,iter,replace=T),
  sample(archaicintrog,iter,replace=T),
  sample(admix_w_e_t6,iter,replace=T),
  sample(t9,iter,replace=T),
  sample(gor_ghost_anc,iter,replace=T),randtyp[1:iter])))

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

config<-c(54,2,18,24)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
len=40000


all_function<-function(nput) {
  msout <-read.ms.output(paste("/dev/shm/mydata/wanc.ghost.abc.sim_",nput,sep=""))
  
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
  # heterozygosity per kbp, which is a more accesible value
#------------------------------------------------------------------------------------------------------------------------

cofig2<-cbind(c(0,cumsum(config)[-length(config)]),cumsum(config))

 ## population fixed sites & segregating sites
segfun<-function(y,ids) {
      idx<-msout$gametes[[y]]
      pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]);op<-table(rowSums(idx[,pop,drop=F]))
      op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)
      op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
      c(op1,op2)}
      apply(cofig2,1,pofun,idx=idx)
  }
  full_segs=do.call(rbind,lapply(1:msout[grep("nreps",names(msout))][[1]],segfun,ids=ids))
  full_segsites=rbind(colSums(full_segs[seq(1,4999,2),]),colSums(full_segs[seq(2,5000,2),])) 

  ## individual fixed sites
  count_fun_f<-function(y,ids) {
    count_fun2_f<-function(i,y,ids) { length(which(ids[[y]][[i]]=="11")) }
    unlist(lapply(1:length(vec),count_fun2_f,ids=ids,y=y))
  } 
 
  full_f=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun_f,ids=ids)
  full_fix=rowSums(do.call(cbind,full_f))/(msout[grep("nreps",names(msout))][[1]]*len)
  mean_fun_f<-function(typ,ffix) { c(mean(ffix[(typ[1]+1):typ[2]]),sd(ffix[(typ[1]+1):typ[2]])) }
  fix_op<-apply(cofig,1,mean_fun_f,ffix=full_fix)*1000

# full_segsites is the # of population-wise fixed sites and the # of population-wise segregating sites; 
 # fix_op is the fixed sites per individual (same format as heterozygotes in your script).

#------------------------------------------------------------------------------------------------------------------------
  
  ## fst here
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
#  set negative fst vals to 0 (considered loci with no popn differentiation)
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


# Wed 23 Mar 2022 14:24:36 CET # updated tajimas d stat below

fout.pi<-rbind(
loop.pi(1,54),
loop.pi(55,56),
loop.pi(57,74),
loop.pi(75,98)
)

#----------------------------------------------------------------------------------------
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

# output both ways
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
initial_ghost=input[1]
time1=t1
time1_el=e_lowl_t1
time2=t2
time2_el=e_lowl_t2
time3=t3
time3_em=e_moun_t3
time3.1=time3+0.0025
time3.1_em=input[2]
time4=input[3]
time4_ance=e_anc_t4
timeI=input[4]
timeI_archgeneflow=input[5]
timeI.1=timeI+0.00025
time5=t5
time5_wl=w_lowl_t5
time6=t6
time6_geneflow_w_e=input[6]
time6_geneflow_e_w=admix_e_w_t6
time6.1=time6+0.00025
time7=t7
time7_ancw=w_anc_t7
time8=t8
time8_ancew=gor_anc
time9=input[7]
time9_allanc=input[8]
#------------------------------------------------------------------------------------------------------------------------
# generate simulation values, i.e. a list of complete parameters to input to the simulation
simval<-list(
c("-en",time1,paste(3,time1_el,sep=" ")),
c("-en",time2,paste(3,time2_el,sep=" ")),
c("-en",time3,paste(4,time3_em,sep=" ")),
c("-en",time3.1,paste(4,time3.1_em,sep=" ")),
c("-ej",time4,paste(4,3,sep=" ")),
c("-en",time4,paste(3,time4_ance,sep=" ")),
c("-em",timeI,paste(1,5,timeI_archgeneflow,sep=" ")),
c("-em",timeI.1,paste(1,5,0,sep=" ")),
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

# ABC + ghost
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 2500 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 5 54 2 18 24 0 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em," -n 5 ",initial_ghost,"  ",simval," > /dev/shm/mydata/wanc.ghost.abc.sim_",input[9],sep=""),intern=F)

#------------------------------------------------------------------------------------------------------------------------

return(all_function(input[9]))

}

# the following command generates simulated data for the "iter" number of times and interprets the output; it divides this task to 6 cores on the cluster, and collects the data of these "iter" iterations into a single object
simuresults<-mclapply(inputsets,wholefunction,mc.cores=6,mc.silent=F)

save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/ghost/final_1mar22/wanc/wanc.ghost.abc_sim",ptype,sep=""))


for (i in 1:iter){
# ABC + ghost
system(paste("rm /dev/shm/mydata/wanc.ghost.abc.sim_",inputsets[[i]][[9]],sep=""),intern=F)
}
