#!/usr/bin/r
# adapting simul_archi_corr.R MK script for pan to the gorillas, using demog as detailed in msprime.gordemog.wcem1.py
#module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix - call from simul_gor.arr
#source venv2/bin/activate
#ID=$SLURM_ARRAY_TASK_ID

#------------------------------------------------------------------------------------------------------------------------

require(data.table)
## see below what ptype stands for
  # ptype = array number
ptype=(commandArgs(TRUE))
ptype=as.numeric(as.character(ptype))
print(ptype)

## 1) create matrix of mutation rates and recombination rates

#------------------------------------------------------------------------------------------------------------------------

# Generation time, mutation rate and recombination rate
#gen_time = 19
#rec_rate = 9.40e-9 # rho
#mu = 1.235e-8
# note from the Besenbacher paper (Fig 3). Here, the divergence time between humans-gorillas is inferred at 13mya, in this case should I adjust the gorilla parameters to those for the mcanus gphocs using a 12mya divergence (the closest to besenbacher's)?
  # => change parameters to the 12 mya gphocs estimates

## this mu is closest to mcmanus et al. 2015, table 1 ghocs parameters estimated using 8mya divergence times
## but note the mcmanus mt rate is excluding CpGs => an underestimate of the true mu 
## & an underestimate relative to our estimate of mu used in sprime
  # => use all parameters from mcmanus 10mya divergence time b/n human-gor 
  # which gives mt rate per gen w/out cpg as 1.169
#------------------------------------------------------------------------------------------------------------------------

options(scipen=100)
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

#head(tetaro)
#         [,1]  [,2]
#[1,] 2.035100 0.001
#[2,] 1.802431 2.000
#[3,] 1.941010 0.001
#[4,] 1.482309 2.000
#[5,] 1.905777 0.001
#[6,] 2.270180 1.000

  # mkdir -p /scratch/devel/hpawar/admix/sstar/simul
#write.table(tetaro,file="/project/devel/mkuhlwilm/sstar/simuls/tetaroAr_corr.txt",sep="\t",row.names=F,col.names=F,quote=F)
#write.table(tetaro,file="/scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt",sep="\t",row.names=F,col.names=F,quote=F)
  # has already been written out

#cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | wc -l
#20000

## create this object once to store the data - need to run
#gmval<-list()
#save(gmval,file=paste("/project/devel/mkuhlwilm/sstar/simulplot/NCgmval_corr",sep=""))


# this package is needed for calculation
library(mgcv)

#-----------------------------------------------------------------------------------------------------------------------

# create data for this range of values (numbers of mutations in each window), here: 15 to 800 in steps of 5
# ptype is the iteration of this, i.e. the step for which a simulation is performed, so it is parallelized
#val2=seq(15,800,5)[ptype]
# below, creating 20,000 short windows, each containing "val2" number of segregating sites (forcefully), but following the assumptions of the demographic model

val2=seq(15,800,5)[ptype]


types <- data.table(
  out_pop = c(1,1,3,3), 
  in_pop = c(3,4,1,2), 
  outgroup = c("WL","WL","EL","EL"),
  ingroup = c("EL","EM","WL","WC"),
  out_chr = c(44,44,18,18),
  in_chr = c(18,24,44,2)
)

types$sum <- types$out_chr + types$in_chr

#types
#   out_pop in_pop outgroup ingroup out_chr in_chr sum
#1:       1      3       WL      EL      44     18  62
#2:       1      4       WL      EM      44     24  68
#3:       3      1       EL      WL      18     44  62
#4:       3      2       EL      WC      18      2  20

# sample configuration: 4 populations : 44 wes_lowl (pop1), 2 wes_cros (pop2), 18 eas_lowl (pop3), 24 eas_moun (pop4) 
# add col of how many haploid chr this corresponds to

comps<-list(c(44,0,18,0),c(44,0,0,24),c(44,0,18,0),c(0,2,18,0))
# comps
#[[1]]
#[1] 44  0 18  0

#[[2]]
#[1] 44  0  0 24

#[[3]]
#[1] 44  0 18  0

#[[4]]
#[1]  0  2 18  0

#comps[[1]]
#[1] 44  0 18  0


#outputs<-list()
#for (i in 1:4){
#isu<-types[i,7]
#print(isu)
#ivec<-comps[[i]]
#print(ivec)
#}
#  sum
#1:  62
#[1] 44  0 18  0
#   sum
#1:  68
#[1] 44  0  0 24
#   sum
#1:  62
#[1] 44  0 18  0
#   sum
#1:  20
#[1]  0  2 18  0



# test this works
#------------------------------------------------------------------------------------------------------------------------

outputs<-list()
for (i in 1:4){
isu<-types[i,7]
#print(isu)
ivec<-comps[[i]]
#print(ivec)
# generate the ms simulation 
# send 1 test job - see if output is generated
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 1 -seeds ",paste(sample(2000:6000,1),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",ivec," -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/test/sim_",ptype,sep=""),intern=F)
# apply S* to simulated vcf
op<-system(paste("python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py  -vcf /scratch/devel/hpawar/admix/sstar/simul/test/sim_",ptype," -ms -target-pops ",types[i,2],"  -ref-pops ",types[i,1]," -p 5000 -s-star -winlen 40000 -winstep 40000 --ms-pop-sizes 4 ",ivec," --ms-num-diploid-inds ",isu/2," --no-pvalues -mssimlen 40000 | cut -f 1-2,7,9,4 --output-delimiter=' '",sep=""),intern=T)
outputs[[types[i,4]]]<-do.call(rbind,strsplit(op,split=" "))
}
#------------------------------------------------------------------------------------------------------------------------

# is it permissible to access & activate a virtual environment from within R?

## testing purpose before running for everything
save(outputs,file="/scratch/devel/hpawar/admix/sstar/simul/test/test")
#load(file="/scratch/devel/hpawar/admix/sstar/simul/test")
## question would be: given a certain no of sites, what is the expectation for S*?

q()
#------------------------------------------------------------------------------------------------------------------------

# need to adapt the postprocessing **

# process the output for the different comparisons
cmvec<-list()
for (cm in (1:length(comps))) {
  cn=cns[[comps[[cm]][1]]][comps[[cm]][3]]
  staroutpa<-outputs[[cn]][-1,]
  cn=cns[[comps[[cm]][2]]][comps[[cm]][4]]
  staroutve<-outputs[[cn]][-1,]
  staroutpa[,1]<-paste(staroutpa[,1],staroutpa[,2])
  staroutve[,1]<-paste(staroutve[,1],staroutve[,2])
  starperind<-list()
  indiv<-unique(staroutve[,4])
  for (ind in (1:length(indiv))) { starperind[[ind]]<- merge(staroutpa[which(staroutpa[,4]==indiv[ind]),c(1,5,3)],staroutve[which(staroutve[,4]==indiv[ind]),c(1,5,3)],by.x=1,by.y=1,all=F) }
  allstars<-do.call(rbind,starperind)[,c(1,2,4,3,5)]
  allstars[,2]<-as.numeric(as.character(allstars[,2]));    allstars[,3]<-as.numeric(as.character(allstars[,3]));    allstars[,4]<-as.numeric(as.character(allstars[,4]));allstars[,5]<-as.numeric(as.character(allstars[,5]))
  if (length(which(allstars[,2]<1))>0) { allstars[which(allstars[,2]<1),2]<-1 }
  if (length(which(allstars[,3]<1))>0) { allstars[which(allstars[,3]<1),3]<-1 }
  
  # save the values for 2,3 and 4 standard deviations from the mean
  # that is, 95% CI, 99% CI, and 99.5% CI
  lva1<-c((1.96*sd(allstars[,2])+mean(allstars[,2])),(2.33*sd(allstars[,2])+mean(allstars[,2])),(2.575*sd(allstars[,2])+mean(allstars[,2])))
  lva2<-c((1.96*sd(allstars[,3])+mean(allstars[,3])),(2.33*sd(allstars[,3])+mean(allstars[,3])),(2.575*sd(allstars[,3])+mean(allstars[,3])))
  cmvec[[cm]]<-list(lva1,lva2)
}

load(file=paste("/project/devel/mkuhlwilm/sstar/simulplot/NCgmval_corr",sep=""))
gmval[[ptype]]<-cmvec
save(gmval,file=paste("/project/devel/mkuhlwilm/sstar/simulplot/NCgmval_corr",sep=""))
# to clean up, remove the large simulated data files (they should be replicable because the number of windows is large)
system(paste("rm /scratch/devel/mkuhlwilm/sstar/simuls/sim_arch",ptype,sep=""),intern=F)

print("nicely done!")
print(Sys.time())
q()

scp /home/mkuhlwilm/Documents/Programs/archaicadmix/simul_archi_corr.R 172.16.10.20:/home/devel/mkuhlwilm/programs/simul_model.R
scp /home/mkuhlwilm/Documents/Programs/archaicadmix/simul_mod.sh 172.16.10.20:/home/devel/mkuhlwilm/programs/simul_mod.sh
ssh 172.16.10.20
mnsubmit /home/devel/mkuhlwilm/programs/simul_mod.sh

## make a gam object from these quantiles, so you can later see if there are "real" segments sticking out
# interactive session
#module load PYTHON/2.7.3 R/3.2.0 tabix 
options(scipen=100)
library(mgcv)
cns<-list(c("WC","CW"),c("WB","BW","CB","BC"))
comps<-list(c(2,1,4,1),c(2,1,2,2),c(2,2,3,1))
nams<-list(c("Central to bonobo","Central to Western"),c("Western to bonobo", "Western to Central"),c("Bonobo to Central","Bonobo to Western"))
hdr<-c("Central chimp","Western chimp","Bonobo")
vals=seq(25,700,5)
load(file=paste("/project/devel/mkuhlwilm/sstar/simulplot/NCgmval_corr",sep=""))

# basically, make a full model for the expectation of S* at any given number of segregating sites, considering the XX% CI defined above 
allgms<-list()
for (cm in (1:length(comps))) {
  gmseta<-list();gmsetb<-list()
  cn1=cns[[comps[[cm]][1]]][comps[[cm]][3]]
  cn2=cns[[comps[[cm]][2]]][comps[[cm]][4]]
  for (vl in (1:length(vals))) {
    gmseta[[vl]]<-c(gmval[[vl]][[1]][[1]],vals[vl])
    gmsetb[[vl]]<-c(gmval[[vl]][[1]][[2]],vals[vl])
  }
  allgms[[cm]]<-list(do.call(rbind,gmseta),do.call(rbind,gmsetb))
}
allmods<-list()
for (cm in (1:length(comps))) {
  allmods[[cm]]<-list()
  for (ty in (1:2)) {
    allmods[[cm]][[ty]]<-list()
    for (qn in c(1:3)) {
      sA<-allgms[[cm]][[ty]][,qn];sS<-allgms[[cm]][[ty]][,4]
      allmods[[cm]][[ty]][[qn]]<-gam(sA~te(sS,k=10), data=list(sA,sS), method="GCV.Cp")
    }
  }
}

save(allmods,file="/project/devel/mkuhlwilm/sstar/simulplot/NCgmsum_2")    

###############################################################
### analyze the real data with the expectation gained above
# interactive session
module load PYTHON/2.7.3 R/3.2.0 tabix gcc
cd /project/devel/mkuhlwilm/sstar/
R --vanilla
options(scipen=100)
library(mgcv)

# these combinations are for chimp/bonobo, different in each species
cn<-unlist(read.table("/project/devel/mkuhlwilm/sstar/comball.txt",sep="\t",header=F)[,1])
comas<-list(c("PTT","PTV"),c("PTT","PTE"),c("PTT","PTS"),c("PTS","PTV"),c("PTS","PTE"),c("PTE","PTV"))
cns<-list(cn[1:20][grep("PPA",cn[1:20])],cn[1:20][-grep("PPA",cn[1:20])])
comas<-list(c("PTT","PTV"),c("PTT","PTE"),c("PTT","PTS"),c("PTS","PTV"),c("PTS","PTE"),c("PTE","PTV"))
chroms=1:23
#chroms=23
medrel<-list()
# Using the most conservative CI (99.5%)
load(file=paste("/project/devel/mkuhlwilm/sstar/simulplot/NCgmsum_2",sep=""))    
qn=3

for (chrom in (chroms)) {
  sdvec[[chrom]]<-list();medrel[[chrom]]<-list() }

daset<-list()
sdvec<-list();gna<-c()
##for (xx in (1:length(comas))) {
## some information and names
xx=1
medrel[[xx]]<-list();sdvec[[xx]]<-list()
sus<-comas[[xx]]
suna<-gsub("PTT","Central",sus);  suna<-gsub("PTV","Western",suna);  suna<-gsub("PTS","Eastern",suna);  suna<-gsub("PTE","NigCam",suna);  gnam<-paste(suna[1],suna[2],sep="_");  print(gnam);gna[xx]<-gnam
comps<-list(c(1,2,which(cns[[1]]==paste("PPA-",sus[1],sep="")),which(cns[[2]]==paste(sus[2],"-",sus[1],sep=""))),c(1,2,which(cns[[1]]==paste("PPA-",sus[2],sep="")),which(cns[[2]]==paste(sus[1],"-",sus[2],sep=""))),c(1,1,which(cns[[1]]==paste(sus[1],"-PPA",sep="")),which(cns[[1]]==paste(sus[2],"-PPA",sep=""))))
nams<-list(c(paste(suna[1]," to Bonobo",sep=""),paste(suna[1]," to ",suna[2],sep="")),c(paste(suna[2]," to Bonobo",sep=""), paste(suna[2]," to ",suna[1],sep="")),c(paste("Bonobo to ",suna[1],sep=""),paste("Bonobo to ",suna[2],sep="")))

# read in all comparisons for all chromosomes, and process them
# in this case, there are two S* calculations because there are three populations
# WC CC and BO, and in each case S* is calculated with using either of the other as outgroups
for (cm in (1:length(comps))) {
  print(cm)
  # first part
  cn=cns[[comps[[cm]][1]]][comps[[cm]][3]]
  starout<-list()
  for (chrom in (chroms)) {
    starout[[chrom]]<-read.table(paste("/scratch/devel/mkuhlwilm/sstar/all/mac_",ifelse(chrom==23,"","output_"),cn,"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] }
  staroutpa<-do.call(rbind,starout)
  # second part (identical with differnet outgroup)
  cn=cns[[comps[[cm]][2]]][comps[[cm]][4]]
  starout<-list()
  for (chrom in (chroms)) {
    starout[[chrom]]<-read.table(paste("/scratch/devel/mkuhlwilm/sstar/all/mac_",ifelse(chrom==23,"","output_"),cn,"_",chrom,".star",sep=""),sep="\t",header=T)[,c(9,4,10,52,12,13,7,1,2)] }
  staroutve<-do.call(rbind,starout)
  # merge them and get values per individual
  staroutpa[,8]<-paste(staroutpa[,8],staroutpa[,9])
  staroutve[,8]<-paste(staroutve[,8],staroutve[,9])
  starperind<-list()
  indiv<-unique(staroutve[,7])
  # here, I only want to use windows where there is enough data (33,333 bp out of 40,000, with at least 30 segregating sites)
  for (ind in (1:length(indiv))) {
    starperind[[ind]]<- merge(staroutpa[which(staroutpa[,7]==indiv[ind] & staroutpa[,4]>=33333 &staroutpa[,2]>=30),c(1,8,2)],staroutve[which(staroutve[,7]==indiv[ind] & staroutve[,4]>=33333 &staroutve[,2]>=30),c(1,8,2)],by.x=2,by.y=2,all=F) }
  allstars<-do.call(rbind,starperind)[,c(1,2,4,3,5)]
  if (length(which(allstars[,2]<1))>0) { allstars[which(allstars[,2]<1),2]<-1 }
  if (length(which(allstars[,3]<1))>0) { allstars[which(allstars[,3]<1),3]<-1 }
  # jitter to make data less stepwise
  allstars[,4]<-as.numeric(jitter(allstars[,4]))
  allstars[,5]<-as.numeric(jitter(allstars[,5]))

  ## either choose a subset of windows (for thinning, maybe for plotting etc), or everything
  #    subse=sample(1:nrow(allstars),50000)
  subse=1:nrow(allstars)
  as2<-allstars[subse,2];as3<-allstars[subse,3];as4<-allstars[subse,4];as5<-allstars[subse,5]
  # look at number of segregating sites and S* values
  newdatA=data.frame(as3=allstars[,3]);newdatB=data.frame(as2=allstars[,2])
  newdatA2=data.frame(sS=allstars[,5]);newdatB2=data.frame(sS=allstars[,4])

  # get the predicted S* for each window given the segregating sites
  glA2<-predict(allmods[[cm]][[1]][[qn]],newdata=newdatA2,type="response")
  glB2<-predict(allmods[[cm]][[2]][[qn]],newdata=newdatB2,type="response")
  
  # get the difference to real S* value
  ldif2<-glA2-allstars[,2]
  ldif3<-glB2-allstars[,3]

  # which windows are outside the expectation for the XX% CI for both comparisons?
  fset3<-which( allstars[,2]>glA2 &allstars[,3]>glB2)
  
  # summary stats
  sdset<-round(c(length(fset1)/nrow(allstars),length(fset2)/nrow(allstars),length(fset3)/nrow(allstars))*100,3)
  mdset<-c(mean(allstars[,2]),mean(allstars[,3]),sd(allstars[,2]),sd(allstars[,3]),mean(allstars[,4]),mean(allstars[,5]),sd(allstars[,4]),sd(allstars[,5]))
  # some values for plotting
  xr=range(allstars[,2]);yr=range(allstars[,3])#[-which(allstars[,2]%in%c("Inf","-Inf")),2],na.rm=T)
  # object to save data
  daset[[cm]]<-list(allstars,fset1,fset2,fset3,xr,yr,sdset,mdset)
}
save(daset,file="/project/devel/mkuhlwilm/sstar/maindataX.Robject")

## plotting and downstream analysis from there
