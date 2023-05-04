#!/usr/bin/r
## Wed 24 Feb 2021 07:30:23 CET
## generate 'high divergence' simulations (& glms) pushing back the divergence time b/n E & W clades of gorilla (to mailund estimate - 429,000)

#------------------------------------------------------------------------------------------------------------------------

require(data.table)
ptype=(commandArgs(TRUE))
ptype=as.numeric(as.character(ptype))
print(ptype)
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
recv<-rbm*4000*slen
# scaled recombination rate. R * 4 * N * length of fragment

# sample mutation rate from normal dist (rnorm(n, mean = 0, sd = 1))
# sample recombination rate from -ve binomial dist (sampling parameter vals for input into ms from prob dist using tbs flag)

tetaro<-cbind(rnorm(20000,mr,sdv),rnbinom(20000,mu=recv,size=0.5))

tetaro[,1]<-ifelse(tetaro[,1]<0.001,0.001,tetaro[,1])
tetaro[,2]<-ifelse(tetaro[,2]<0.001,0.001,tetaro[,2])

val2=seq(15,800,5)[ptype]


types <- data.table(
  out_pop = c(1,1,3,3)-1,
  in_pop = c(3,4,1,2)-1,
  outgroup = c("WL","WL","EL","EL"),
  ingroup = c("EL","EM","WL","WC"),
  out_chr = c(44,44,18,18),
  in_chr = c(18,24,44,2)
)

types$sum <- types$out_chr + types$in_chr

comps<-list(c(44,0,18,0),c(44,0,0,24),c(44,0,18,0),c(0,2,18,0))

#------------------------------------------------------------------------------------------------------------------------


outputs<-list()
for (i in (1:nrow(types))) {
# define number of individuals to simulate in each comparison
isu<-types[i,7]
ivec<-comps[[i]]
# generate the ms simulation 
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 20000  -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4], " -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 5.644737 3 1 -en 5.644737 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/out/high.divergence/sim_",ptype,sep=""),intern=F)
dip<-isu/2
# apply S* to simulated vcf
op<-system(paste("python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py  -vcf /scratch/devel/hpawar/admix/sstar/simul/out/high.divergence/sim_",ptype," -ms -target-pops ",types[i,2],"  -ref-pops ",types[i,1],"    -p 5000 -s-star -winlen 40000 -winstep 40000 --ms-pop-sizes 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4]," --ms-num-diploid-inds ",dip," --no-pvalues -mssimlen 40000 | cut -f 1-2,7,9,4 --output-delimiter=' '",sep=""),intern=T)
outputs[[i]]<-do.call(rbind,strsplit(op,split=" "))
}


#------------------------------------------------------------------------------------------------------------------------


## given a certain no of sites, what is the expectation for S*?

# process the output for the different comparisons

out<-list()
for (i in (1:nrow(types))) {
# 95% CI, 99% CI, and 99.5% CI
out[[i]]<-c((1.96*sd(as.numeric(outputs[[i]][-1,5]))+mean(as.numeric(outputs[[i]][-1,5]))),(2.33*sd(as.numeric(outputs[[i]][-1,5]))+mean(as.numeric(outputs[[i]][-1,5]))),(2.575*sd(as.numeric(outputs[[i]][-1,5]))+mean(as.numeric(outputs[[i]][-1,5]))))
}

proc<-as.data.frame(rbind(out[[1]],out[[2]],out[[3]],out[[4]]))

# only once
#savelist<-list()
#save(savelist,file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/high.divergence/sstar_dist",sep=""))


# in every job
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/high.divergence/sstar_dist",sep=""))
savelist[[ptype]]<-proc;names(savelist[[ptype]])<-val2
save(savelist,file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/high.divergence/sstar_dist",sep=""))
