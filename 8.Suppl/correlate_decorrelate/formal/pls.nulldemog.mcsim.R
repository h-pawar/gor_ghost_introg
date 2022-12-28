#!/usr/bin/r

#Mon 21 Feb 2022 12:33:39 CET
# generate model comparison simulations for PLS ABC posteriors ('reduced dimensionality window ABC')

# amending below structure from 
# model A) null demography
#/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/test.null.demog.v1.R
#/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/abc+ghost/modelcomp/test.null.demog.arr


# Wed 24 Nov 2021 10:37:57 CET
# for null demog model (only the 4 extant populations)
#-----------------------------------------------------------------------------------------------------------------------

# Mon 22 Nov 2021 15:26:23 CET
# for model comparison need to generate x simulations per model, obtain summary stats & do the comparison

# MK
#You should rather simulate a larger number (maybe 25000) for each demographic model. 
#In principle, you can re-use parts of the code you used to generate the models: 
#Take the posterior parameters for each model as the input (instead of random values), obtain the same summary stats as before as output. 
#This would be again a lot of simulated data, so you may want to use perhaps 250 or 500 instead of 2500 windows per simulation 
#(just make sure the way you calulate statistics is not affected).

# me
#Would it be ok for me to do 2500 reps of 2500 windows (for consistency with the previous scripts) for each model, rather than 25,000 reps of 250 windows? 

# MK
#the power really comes with the number of replicates. But do 2500 reps first and see how it looks like, probably it will be fine already.

# MK
#you run the same simulation again and again to get a distribution of values. You could still run the iterations per job to collect the data, 

#-----------------------------------------------------------------------------------------------------------------------
# following postabc.sim.model.v1.R
# & adding calcn of summary stats (equivalent to test.ghoste.demog.v1.R & test.ghostw.demog.v1.R) & equivalent to prev scripts using all function & whole function
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
## see below what ptype stands for
  # ptype = array number
ptype=(commandArgs(TRUE))
#ptype=1 # for troubleshooting purposes
ptype=as.numeric(as.character(ptype))
print(ptype)
options(scipen=100)

# number of iterations per run of script # Q - how to choose this? * # MK: recommend that you may use iter=102 instead to target ~100,000 simulations.  
#iter=1 # run as iter=1 to get sense of time, then for actual jobs return to iter=102

# MK: Then, try to run the script with iter=6 and measure how long that takes. 
#iter=6

# Thu  3 Jun 2021 12:49:00 CEST MK - If half the data takes 40min, you can launch jobs with 51 iterations on 4 CPUs, which should take ~8.5h, so you can set 9h for the time. 
#iter=51
iter=50
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
# standard deviation of the mutation rate - whether this should indeed be *0.233 for gor? ask ** use as a starting point
recv<-rbm*4000*slen
# scaled recombination rate. R * 4 * N * length of fragment

# sample mutation rate from normal dist (rnorm(n, mean = 0, sd = 1))
# sample recombination rate from -ve binomial dist (sampling parameter vals for input into ms from prob dist using tbs flag)

tetaro<-cbind(rnorm(20000,mr,sdv),rnbinom(20000,mu=recv,size=0.5))
  # rnorm(n, mean = 0, sd = 1) n=sample size, ie here number of reps to simulate - 20,000 reps
  # rnbinom(n, size, prob, mu) # where prob here is probability of success in each trial - https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html
      # size = dispersion parameter
      # mu = mean value for the -ve binomial

tetaro[,1]<-ifelse(tetaro[,1]<0.001,0.001,tetaro[,1])
tetaro[,2]<-ifelse(tetaro[,2]<0.001,0.001,tetaro[,2])

types <- data.table(
#  out_pop = c(1,1,3,3),
#  in_pop = c(3,4,1,2),
  out_pop = c(1,1,3,3)-1,
  in_pop = c(3,4,1,2)-1,
  outgroup = c("WL","WL","EL","EL"),
  ingroup = c("EL","EM","WL","WC"),
  out_chr = c(54,54,18,18),
  in_chr = c(18,24,54,2)
)
# should be 27wl instead of 22 - Thu 10 Jun 2021 14:54:59 CEST
# 27*2 = 54, changing 44 to 54 & 22 to 27

types$sum <- types$out_chr + types$in_chr

comps<-list(c(54,0,18,0),c(54,0,0,24),c(54,0,18,0),c(0,2,18,0))

# this number of haploid chr to simulate: w_lowl, w_cros, e_lowl, e_moun
ivec<-c(54,2,18,24);isu=sum(ivec);ivec=paste(ivec,collapse=" ")

#------------------------------------------------------------------------------------------------------------------------
#> summary(myabc_6oct)
# => use the weighted medians from posteriors for demographic parameters
#                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
#Weighted Median:         72.3108   10.1232   13.7748    5.2451   12.3644
#                       e_lowl_t2 e_moun_t3 e_moun_t3.1       t4 e_anc_t4
#Weighted Median:         20.6392    0.2168      3.6443   0.2214  15.0006
#                             t5 w_lowl_t5 admix_w_e_t6 admix_e_w_t6       t7
#Weighted Median:         0.5072   28.9897      79.9123      41.8949   4.1949
#                       w_anc_t7       t8  gor_anc
#Weighted Median:        58.4250  12.6360  17.3891
#------------------------------------------------------------------------------------------------------------------------

# fixed parameters from test.abc.model.v5.R
t1=0.0003947368
t2=0.095
t3=0.12
t6=0.4473684

# posteriors from generateabc.6oct.R
# 1) current Ne of populations 
#initial_wl=72.3108
#initial_wc=10.1232
#initial_el=13.7748
#initial_em=5.2451
# 2) rapid & severe population decline of E lowland over last 20y (0.2kya) of ~80% - van der valk 2018, plumptr 2016
#time1=t1
#time1_el=12.3644
# 3)  Grauer’s gorillas went through a period of population growth and range expansion 5,000–10,000 years ago [22, 23]' van der valk 2019 citing van der valk 2018 & Tocheri 2016
#time2=t2
#time2_el=20.6392
# 4) severe population bottleneck experienced by mountain gorillas - xue et al - this was likely after their split from eastern lowland
#time3=t3
#time3_em=0.2168
# MK: 148 You may model this bottleneck as a temporal change, where at t3 it goes down to this value, and t3 minus 5 or 10 generations it goes back to the previous/terminal Ne.
#time3.1=time3+0.0025
#time3.1_em=3.6443
# 5) E subspecies split # merge eas_moun into eas_lowl
#time4=0.2214
#time4_ance=15.0006
# 6) recent contraction in population size 23 ka in w lowland - mcmanus et al. 2015
#time5=0.5072
#time5_wl=28.9897
# 7) migration pulse W > E 3% - for admixture need to give as a proportion not percent
#time6=t6
#time6_geneflow_w_e=79.9123
#time6_geneflow_e_w=41.8949
# The migration pulse should happen for only one generation and then stop (that is, adding 19/4000 years as one generation)
#time6.1=time6+0.00025
# 8) W subspecies split ie merge wes_cros into wes_lowl
#time7=4.1949
#time7_ancw=58.4250
## 9) gor split (W - E) # merge eastern into western 
#time8=12.6360
# population size of gorilla anc (common anc of the 2 species)
#time8_ancew=17.3891


#------------------------------------------------------------------------------------------------------------------------
# posteriors from PLS_ABC (running abc using the decorrelated summary statistics ie the 10 PLS components)
# posteriors from post.pls.abc.18feb22.R
# 1) current Ne of populations 
initial_wl=89.2094
initial_wc=20.2728
initial_el=16.8704
initial_em=2.9892
# 2) rapid & severe population decline of E lowland over last 20y (0.2kya) of ~80% - van der valk 2018, plumptr 2016
time1=t1
time1_el=4.6104
# 3)  Grauer’s gorillas went through a period of population growth and range expansion 5,000–10,000 years ago [22, 23]' van der valk 2019 citing van der valk 2018 & Tocheri 2016
time2=t2
time2_el=18.3987
# 4) severe population bottleneck experienced by mountain gorillas - xue et al - this was likely after their split from eastern lowland
time3=t3
time3_em=1.4849
# MK: 148 You may model this bottleneck as a temporal change, where at t3 it goes down to this value, and t3 minus 5 or 10 generations it goes back to the previous/terminal Ne.
time3.1=time3+0.0025
time3.1_em=10.5757
# 5) E subspecies split # merge eas_moun into eas_lowl
time4=0.2179
time4_ance=22.9690
# 6) recent contraction in population size 23 ka in w lowland - mcmanus et al. 2015
time5=0.4769
time5_wl=45.6855
# 7) migration pulse W > E 3% - for admixture need to give as a proportion not percent
time6=t6
time6_geneflow_w_e=17.0756
time6_geneflow_e_w=39.6015
# The migration pulse should happen for only one generation and then stop (that is, adding 19/4000 years as one generation)
time6.1=time6+0.00025
# 8) W subspecies split ie merge wes_cros into wes_lowl
time7=7.7345
time7_ancw=86.2332
## 9) gor split (W - E) # merge eastern into western 
time8=12.8414
# population size of gorilla anc (common anc of the 2 species)
time8_ancew=35.0938

#------------------------------------------------------------------------------------------------------------------------
  # may need this to retain identity of the iter 
inputsets=as.data.frame(t(randtyp[1:iter]))
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# code from het.simld.15apr.R - MK amended - to be more efficient
# nested function to calculate heterozygosity & seg sites & fst: 
  # bhatiafst function from - readme_fstcomponents (script written with JS)

config<-c(54,2,18,24)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
#sum(config)
len=40000

all_function<-function(nput) {
  msout <-read.ms.output(paste("/dev/shm/mydata/mc.nullpls.sim_",nput,sep=""))
  ## create empty vector
    # bind pairs of cols, ie 1_2, 3_4, 5_6 etc (the haploids -> diploid ids)
    # MK: You would divide the number of hets by all sites for which you have data, in this case 40kbp (times the number of simulated windows)
    # this is for the 1 simn rep - will need to add another layer for going through all the simn reps
    # loop through every i with step size of 2 again
      # bind pairs of cols, ie 1_2, 3_4, 5_6 etc (the haploids -> diploid ids)
      # MK: To make it more efficient, you could even paste without the separator, then you get 00 01 10 and 11.
    # remove null elements of the list - https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list
    
  # generalise to any number of reps : # msout[grep("nreps",names(msout))] # https://stackoverflow.com/questions/4220631/how-do-i-grep-in-r
  ## apply the function for all simn reps from 1:msout[grep("nreps",names(msout))][[1]]
#    return(count_het)
#  }

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
  
  # estimating heterozygosity per simn rep - is this the way to go? # add "*1000" to this, to get heterozygosity per kbp, which is a more accesible value
  
  ## MK: output segsites # new way to calculate below
#  segfun<-function(y,ids) { 
#    idx<-do.call(cbind,ids[[y]])
#    gt_fun<-function(nput) { length(which(nput!="00")) }
#    pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]); sum(apply(idx[,pop,drop=F],1,gt_fun)) }
#    apply(cofig,1,pofun,idx=idx)
#    }
#  full_seg=colSums(do.call(rbind,lapply(1:msout[grep("nreps",names(msout))][[1]],segfun,ids=ids)))

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

 #it makes sense to take the population-wise segregating sites as non-fixed sites. This seems like a more informative measure. 
 #This would exclude the CRG population, because this value is the same as heterozygosity with just one individua

#good that you are checking the things carefully! It seems that you need to say
#        op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
#Only add "as.numeric()", the typical R problem... For safety also do that in op1:         
#op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)

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
  full_segsites=rbind(colSums(full_segs[seq(1,499,2),]),colSums(full_segs[seq(2,500,2),])) 

# MK
#1) The output of each iteration in full_segs is a data frame with two rows, each representing a different statistic. 
#This line takes the sum of each other row to give the final value for each statistic.

  ## individual fixed sites
  count_fun_f<-function(y,ids) {
    count_fun2_f<-function(i,y,ids) { length(which(ids[[y]][[i]]=="11")) }
    unlist(lapply(1:length(vec),count_fun2_f,ids=ids,y=y))
  } 
 
  full_f=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun_f,ids=ids)
  full_fix=rowSums(do.call(cbind,full_f))/(msout[grep("nreps",names(msout))][[1]]*len)
  mean_fun_f<-function(typ,ffix) { c(mean(ffix[(typ[1]+1):typ[2]]),sd(ffix[(typ[1]+1):typ[2]])) }
  fix_op<-apply(cofig,1,mean_fun_f,ffix=full_fix)*1000

# to output 
#population-wise fixed sites & segregating sites, individual fixed sites (like heterozygosity). 
# full_segsites, fix_op (yes not full_segs)
 
# MK: full_segsites is the # of population-wise fixed sites and the # of population-wise segregating sites; 
 # fix_op is the fixed sites per individual (same format as heterozygotes in your script).

#------------------------------------------------------------------------------------------------------------------------
  
  ## you may add fst here
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
# add this line into the function for calculating pairwise fst, set negative fst vals to 0 (considered loci with no popn differentiation)
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
# add another layer  - loop through simn reps
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
# output mean & sd
out.pi<-cbind(mean(unlist(hold.pi)), sd(unlist(hold.pi)))
return(out.pi)
}



loop.taj<-function(x,y){
# add another layer  - loop through simn reps
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

out.tajima<-cbind(mean(unlist(hold.tajima)), sd(unlist(hold.tajima)))
return(out.tajima)
}

fout.pi<-rbind(
loop.pi(1,54),
loop.pi(55,56),
loop.pi(57,74),
loop.pi(75,98)
)

fout.taj<-rbind(
loop.taj(1,54),
#loop.pi(45,46)[[2]],
loop.taj(57,74),
loop.taj(75,98)
)

#----------------------------------------------------------------------------------------

  ## output
  #return(list(het_op,full_seg,fst.out,fout.pi,fout.taj))
  return(list(het_op,full_segsites,fix_op,fst.out,fout.pi,fout.taj))

}

#------------------------------------------------------------------------------------------------------------------------

wholefunction=function(input) {

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

## Mon 22 Nov 2021 15:01:45 CET
# MK - This would be again a lot of simulated data, so you may want to use perhaps 250 or 500 instead of 2500 windows per simulation 
#(just make sure the way you calulate statistics is not affected).
  # Wed Nov 24 09:17:51 CET 2021 -> try adjusting to 250  
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 250 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -I 4 54 2 18 24 -n 1 ",initial_wl," -n 2 ",initial_wc," -n 3 ",initial_el," -n 4 ",initial_em,"  ",simval," > /dev/shm/mydata/mc.nullpls.sim_",input[1],sep=""),intern=F)

#------------------------------------------------------------------------------------------------------------------------
# now calc summary statistics for this ms simn - code from het.simld.15apr.R
# Q - whether this shoudl be part of the 'wholefunction' also? * - or the all_function shoudl be outside the wholefunction? but called within it once?
  # that would make more sense 

all_function(input[1])

#system(paste("rm /dev/shm/mydata/abc.sim_",input[23],sep=""),intern=F) # this was throwing errors here

return(all_function(input[1]))

}

# whether this wholefunction will work..

#remove temporary data (important when doing many!!)
# MK: 4) Also, VERY important, delete the things you create on /dev/shm/ !
  # Put this after all_function(input[25]) is done: system(paste("rm /dev/shm/mydata/abc.sim_",input[25],sep=""),intern=F)
  # If you dont do this, the temporary simulation files will just accumulate and accumulate, flooding the file system! 
#system(paste("rm /dev/shm/mydata/abc.sim_",input[23],sep=""),intern=F)

# everything above is done as one function instead of a continuous process in R
# this allow parallelization within the same launching of the script
# the following command generates simulated data for the "iter" number of times and interprets the output; it divides this task to 6 cores on the cluster, and collects the data of these "iter" iterations into a single object
simuresults<-mclapply(inputsets,wholefunction,mc.cores=6,mc.silent=F)
# this object, (importantly) together with the input parameters, is then stored as one object
# here, x simulations are done and being processed
# save new simulations to new dir: /scratch/devel/hpawar/admix/abc/simul/test/11sep21
# instead of /scratch/devel/hpawar/admix/abc/simul/test/out/
  # change dir - done
 # mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/ghost/3nov21
# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/modelcomp

#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/null
#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/ghoste
#[hpawar@cnb1 ~]$ mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/ghostw

#save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/11sep21/abc_sim",ptype,sep=""))

# CHANGE DIR **
#ls -lhtr /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/
# mkdir -p /scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/pls_mcsimns

save(simuresults,inputsets,file=paste("/scratch/devel/hpawar/admix/abc/simul/test/modelcomp/out/pls/pls_mcsimns/test.mc.nullpls_sim",ptype,sep=""))

# this works - remove input ms files for this run here 
for (i in 1:iter){
#system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[i]][[23]],sep=""),intern=F) # have now fixed 4 more parameters
#system(paste("rm /dev/shm/mydata/abc.sim_",inputsets[[i]][[19]],sep=""),intern=F)
# ABC + ghost
system(paste("rm /dev/shm/mydata/mc.nullpls.sim_",inputsets[[i]][[1]],sep=""),intern=F)
}

#-----------------------------------------------------------------------------------------------------------------------

