#!/usr/bin/r
## Tue  2 Feb 2021 12:20:13 CET
# troubleshooting - ms simulated data & S*
#------------------------------------------------------------------------------------------------------------------------

# load the modules & activate venv2 in bash first
#module load gcc/6.3.0 openssl/1.0.2q PYTHON/2.7.17 R/4.0.1 tabix 
#source venv2/bin/activate

# try to run the gorilla demography fromthe cd line - this runs fine & generates output
#cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms 62 1 -seeds 1 -t tbs -r tbs 40000.0 -s 15 -I 4 44 0 18 0 -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/test/sim_1

#(venv2) [hpawar@cnb1 ~]$ head /scratch/devel/hpawar/admix/sstar/simul/test/sim_1
#/home/devel/hpawar/ms/msdir/ms 62 1 -seeds 1 -t tbs -r tbs 40000.0 -s 15 -I 4 44 0 18 0 -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 
#1 0 2

#//	2.03510047898326	0.001
#segsites: 15
#positions: 0.0236 0.0826 0.1109 0.1552 0.2368 0.3540 0.4359 0.4600 0.5961 0.7515 0.7809 0.7844 0.7938 0.9095 0.9868 
#011111010010010
#000000000000100
#011111010010010
#000000000001000

# => the issue must be from how i'm doing something in R - perhaps a syntax error? yes
# go back to R
#------------------------------------------------------------------------------------------------------------------------

require(data.table)
## see below what ptype stands for
  # ptype = array number
ptype=(commandArgs(TRUE))
#ptype=1 # for troubleshooting purposes
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


val2=seq(15,800,5)[ptype]


types <- data.table(
#  out_pop = c(1,1,3,3),
#  in_pop = c(3,4,1,2),
  out_pop = c(1,1,3,3)-1,
  in_pop = c(3,4,1,2)-1,
  outgroup = c("WL","WL","EL","EL"),
  ingroup = c("EL","EM","WL","WC"),
  out_chr = c(44,44,18,18),
  in_chr = c(18,24,44,2)
)

types$sum <- types$out_chr + types$in_chr

comps<-list(c(44,0,18,0),c(44,0,0,24),c(44,0,18,0),c(0,2,18,0))


# ensure works for 1 rep before looping - yes this does
#i=1
#isu<-types[i,7]
#print(isu)
#ivec<-comps[[i]]
#print(ivec)
# generate the ms simulation 
# send 1 test job - see if output is generated
#system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 1 -seeds ",paste(sample(2000:6000,1),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4],	" -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/test/testsim_",ptype,sep=""),intern=F)

#dip<-isu/2
# apply S* to simulated vcf - breaks down here - no output generated
#op<-system(paste("python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py  -vcf /scratch/devel/hpawar/admix/sstar/simul/test/testsim_",ptype," -ms -target-pops ",types[i,2],"  -ref-pops ",types[i,1],"    -p 5000 -s-star -winlen 40000 -winstep 40000 --ms-pop-sizes 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4]," --ms-num-diploid-inds ",dip," --no-pvalues -mssimlen 40000 | cut -f 1-2,7,9,4 --output-delimiter=' '",sep=""),intern=T)

# op
#character(0)

#------------------------------------------------------------------------------------------------------------------------

# from cd line - ms output generated
#[hpawar@login2 ~]$ ls -lhtr /scratch/devel/hpawar/admix/sstar/simul/test/
#total 8.0K
#-rw-r--r-- 1 hpawar devel 1.5K Feb  2 12:26 sim_1
#-rw-r--r-- 1 hpawar devel 1.5K Feb  2 13:25 testsim_1

#[hpawar@login2 ~]$ head /scratch/devel/hpawar/admix/sstar/simul/test/testsim_1 
#/home/devel/hpawar/ms/msdir/ms 62 1 -seeds 2293 -t tbs -r tbs 40000.0 -s 15 -I 4 44 0 18 0 -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 
#2293 0 2

#//	2.03510047898326	0.001
#segsites: 15
#positions: 0.1763 0.2621 0.3391 0.3622 0.4032 0.4418 0.4552 0.4788 0.7205 0.7743 0.7789 0.7870 0.8104 0.9299 0.9761 
#000110010101000
#000100011101000
#100000000000100
#000100011101000

#------------------------------------------------------------------------------------------------------------------------
# check if this is indeed the expected output ** yes
#op
# [1] "chrom winstart n_snps ind_id s_star" "ms1 0 15 i22 0"                     
# [3] "ms1 0 15 i23 0"                      "ms1 0 15 i24 0"                     
# [5] "ms1 0 15 i25 0"                      "ms1 0 15 i26 0"                     
# [7] "ms1 0 15 i27 0"                      "ms1 0 15 i28 0"                     
# [9] "ms1 0 15 i29 0"                      "ms1 0 15 i30 0"
# i'm not sure it has worked properly?
# but at least can start amending the postprocessing steps

 #str(op)
# chr [1:10] "chrom winstart n_snps ind_id s_star" "ms1 0 15 i22 0" ...

#------------------------------------------------------------------------------------------------------------------------

# now try to loop

# works for 1 simulated rep for eahc of the 4 scenarios
#outputs<-list()
#for (i in (1:nrow(types))) {
# define number of individuals to simulate in each comparison
#isu<-types[i,7]
#print(isu)
#ivec<-comps[[i]]
#print(ivec)
# generate the ms simulation 
#system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 1 -seeds ",paste(sample(2000:6000,1),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4], " -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/test/testsim_",ptype,sep=""),intern=F)

#dip<-isu/2
# apply S* to simulated vcf - breaks down here - no output generated
#op<-system(paste("python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py  -vcf /scratch/devel/hpawar/admix/sstar/simul/test/testsim_",ptype," -ms -target-pops ",types[i,2],"  -ref-pops ",types[i,1],"    -p 5000 -s-star -winlen 40000 -winstep 40000 --ms-pop-sizes 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4]," --ms-num-diploid-inds ",dip," --no-pvalues -mssimlen 40000 | cut -f 1-2,7,9,4 --output-delimiter=' '",sep=""),intern=T)
#outputs[[i]]<-do.call(rbind,strsplit(op,split=" "))
#}


#str(outputs)
#List of 4
# $ : chr [1:10, 1:5] "chrom" "ms1" "ms1" "ms1" ...
# $ : chr [1:13, 1:5] "chrom" "ms1" "ms1" "ms1" ...
# $ : chr [1:23, 1:5] "chrom" "ms1" "ms1" "ms1" ...
# $ : chr [1:2, 1:5] "chrom" "ms1" "winstart" "0" ...

#outputs
#[[1]]
#      [,1]    [,2]       [,3]     [,4]     [,5]    
# [1,] "chrom" "winstart" "n_snps" "ind_id" "s_star"
# [2,] "ms1"   "0"        "15"     "i22"    "0"     
# [3,] "ms1"   "0"        "15"     "i23"    "42100" 
# [4,] "ms1"   "0"        "15"     "i24"    "0"     
# [5,] "ms1"   "0"        "15"     "i25"    "0"     
# [6,] "ms1"   "0"        "15"     "i26"    "0"     
# [7,] "ms1"   "0"        "15"     "i27"    "0"     
# [8,] "ms1"   "0"        "15"     "i28"    "0"     
# [9,] "ms1"   "0"        "15"     "i29"    "0"     
#[10,] "ms1"   "0"        "15"     "i30"    "0"     

#[[2]]
#      [,1]    [,2]       [,3]     [,4]     [,5]    
# [1,] "chrom" "winstart" "n_snps" "ind_id" "s_star"
# [2,] "ms1"   "0"        "15"     "i22"    "0"     
# [3,] "ms1"   "0"        "15"     "i23"    "0"     
# [4,] "ms1"   "0"        "15"     "i24"    "0"     
# Q - why the winstart is at 0? & i guess the majority of snps associated with sstar score of 0? is that expected?

# now how to postprocess *

## testing purpose before running for everything
#save(outputs,file="/scratch/devel/hpawar/admix/sstar/simul/test/testsim_sstar_1")
#load(file="/scratch/devel/hpawar/admix/sstar/simul/test/testsim_sstar_1")
## question would be: given a certain no of sites, what is the expectation for S*?


#> outputs[[4]][,5]
#[1] "s_star" "35528" 
#> outputs[[4]][1,]
#[1] "chrom"    "winstart" "n_snps"   "ind_id"   "s_star"  
#> outputs[[4]][-1,]
#[1] "ms1"   "0"     "15"    "i0"    "35528"

# these shoudl be for the 4 diff scenarios 
#comps<-list(c(44,0,18,0),c(44,0,0,24),c(44,0,18,0),c(0,2,18,0))

#types
#   out_pop in_pop outgroup ingroup out_chr in_chr sum
#1:       1      3       WL      EL      44     18  62 # WL outgroup, EL ingroup - 9 
#2:       1      4       WL      EM      44     24  68 # WL outgroup, EM ingroup - 12
#3:       3      1       EL      WL      18     44  62 # EL outgroup, WL ingroup - 22 
#4:       3      2       EL      WC      18      2  20 # EL outgroup, WC ingroup - ie 1 diploid as ingroup

# but how will this scale when generating not 1 simulated rep but 20,000 - may need another layer to the nested list? to avoid confusion?
# could try w 5 reps initially & see what the output is like 

#------------------------------------------------------------------------------------------------------------------------

#system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu,"   5 -seeds ",paste(sample(2000:6000,5),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4], "   -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0   -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/test/testsim_mult",ptype,sep=""),intern=F)
# argument should be -2984 ?
# See msdoc.pdf for explanation of these parameters.
#usage: ms nsam howmany 
 
# but if generating 20,000 reps - do not then need 20,000 seeds?? - ask re this *
# no only 3 seeds specified: (from hudson manual)
#-seeds x1 x2 x3 where the x’s indicate the three integer seed values. (this may have been the cause of the error)

# see if this works
#system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 5 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4], "   -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0   -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/test/testsim_mult",ptype,sep=""),intern=F)
# then using these 5 test reps try to write postprocessing steps - ie extract mean & sd of s* score per demog scenario


# test with 5 reps of the simulations - works fine => can scale up to 20,000 reps
outputs<-list()
for (i in (1:nrow(types))) {
# define number of individuals to simulate in each comparison
isu<-types[i,7]
#print(isu)
ivec<-comps[[i]]
#print(ivec)
# generate the ms simulation 
system(paste("cat /scratch/devel/hpawar/admix/sstar/simul/tetaroAr_corr.txt | /home/devel/hpawar/ms/msdir/ms ",isu," 20000 -seeds ",paste(sample(2000:6000,3),collapse=" ")," -t tbs -r tbs 40000.0 -s ",val2," -I 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4], " -n 1 25.161 -n 2 3.080 -n 3 4.280 -n 4 0.800 -ej 0.1973684 4 3 -em 0.4473684 3 1 1200 -em 0.4473684 3 2 16 -em 0.4521184 3 1 0 -em 0.4521184 3 2 0 -ej 0.8947368 2 1 -en 0.8947368 1 30.693 -ej 3.434211 3 1 -en 3.434211 1 39.751 > /scratch/devel/hpawar/admix/sstar/simul/out/sim_",ptype,sep=""),intern=F)

dip<-isu/2
# apply S* to simulated vcf - breaks down here - no output generated
op<-system(paste("python /scratch/devel/mkuhlwilm/arch/freezing-archer-master/bin/windowed_calculations.py  -vcf /scratch/devel/hpawar/admix/sstar/simul/out/sim_",ptype," -ms -target-pops ",types[i,2],"  -ref-pops ",types[i,1],"    -p 5000 -s-star -winlen 40000 -winstep 40000 --ms-pop-sizes 4 ",comps[[i]][1]," ",comps[[i]][2]," ",comps[[i]][3]," ",comps[[i]][4]," --ms-num-diploid-inds ",dip," --no-pvalues -mssimlen 40000 | cut -f 1-2,7,9,4 --output-delimiter=' '",sep=""),intern=T)
outputs[[i]]<-do.call(rbind,strsplit(op,split=" "))
}


# str(outputs)
#List of 4
# $ : chr [1:46, 1:5] "chrom" "ms1" "ms1" "ms1" ...
# $ : chr [1:61, 1:5] "chrom" "ms1" "ms1" "ms1" ...
# $ : chr [1:111, 1:5] "chrom" "ms1" "ms1" "ms1" ...
# $ : chr [1:6, 1:5] "chrom" "ms1" "ms2" "ms3" ...
# still gives list of 4 ie for each scenario has concatenated the output


#[[4]]
#     [,1]    [,2]       [,3]     [,4]     [,5]    
#[1,] "chrom" "winstart" "n_snps" "ind_id" "s_star"
#[2,] "ms1"   "0"        "15"     "i0"     "0"     
#[3,] "ms2"   "0"        "15"     "i0"     "0"     
#[4,] "ms3"   "0"        "15"     "i0"     "61200" 
#[5,] "ms4"   "0"        "15"     "i0"     "0"     
#[6,] "ms5"   "0"        "15"     "i0"     "0"     

# for all rows of given scenario calculate the mean & sd of the s* val

#> outputs[[4]][,5]
#[1] "s_star" "35528" 


#> outputs[[4]][,5]
#[1] "s_star" "35528" 

#outputs[[4]][-1,5]
#[1] "0"     "0"     "61200" "0"     "0" 

#nrow(outputs[[4]])
#[1] 6
#> nrow(outputs[[1]])
#[1] 46

# sum(as.numeric(outputs[[4]][-1,5]))
#[1] 61200

#mean(as.numeric(outputs[[4]][-1,5]))
#[1] 12240

#sd(as.numeric(outputs[[4]][-1,5]))
#[1] 27369.47

#summary(as.numeric(outputs[[4]][-1,5]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      0       0       0   12240       0   61200


#lva1<-c((1.96*sd(as.numeric(outputs[[4]][-1,5]))+mean(as.numeric(outputs[[4]][-1,5]))),(2.33*sd(as.numeric(outputs[[4]][-1,5]))+mean(as.numeric(outputs[[4]][-1,5]))),(2.575*sd(as.numeric(outputs[[4]][-1,5]))+mean(as.numeric(outputs[[4]][-1,5]))))
#lva1
#[1] 65884.17 76010.87 82716.39

# str(outputs)
#List of 4 # ie outputs = list of 4 (1 list for each scenario, with results of the diff simulated reps concatenated)

#------------------------------------------------------------------------------------------------------------------------


## question would be: given a certain no of sites, what is the expectation for S*?
# process the output for the different comparisons

out<-list()
for (i in (1:nrow(types))) {
# define number of individuals to simulate in each comparison
# save the values for 2,3 and 4 standard deviations from the mean
# that is, 95% CI, 99% CI, and 99.5% CI
out[[i]]<-c((1.96*sd(as.numeric(outputs[[i]][-1,5]))+mean(as.numeric(outputs[[i]][-1,5]))),(2.33*sd(as.numeric(outputs[[i]][-1,5]))+mean(as.numeric(outputs[[i]][-1,5]))),(2.575*sd(as.numeric(outputs[[i]][-1,5]))+mean(as.numeric(outputs[[i]][-1,5]))))
}

#proc<-rbind(out[[1]],out[[2]],out[[3]],out[[4]])

#         [,1]     [,2]     [,3]
#[1,] 13314.29 15650.05 17196.70
#[2,]     0.00     0.00     0.00
#[3,] 22981.90 26627.22 29041.01
#[4,] 65884.17 76010.87 82716.39

# write out proc to file
# could add col of how many segregating sites this corresponds to (for when merging b/n the diff steps of segregating sites)


proc<-as.data.frame(rbind(out[[1]],out[[2]],out[[3]],out[[4]]))

#outputs[[1]][-1,3]
#[1] "15" "15" "15" "15" "15"
#m<-outputs[[1]][-1,3]
#m[1]
#[1] "15"

#proc[,4]=m[1]

#proc
#        V1       V2       V3 V4
#1 13314.29 15650.05 17196.70 15
#2     0.00     0.00     0.00 15
#3 22981.90 26627.22 29041.01 15
#4 65884.17 76010.87 82716.39 15

# output this?

# MK By the way, you dont need a column for the # of sites because that is by definition "val2". 
#You save the output to a file that you create before, with an empty list. 
#Each element of this list should be one "proc" table. 
#You can say "savelist[[ptype]]<-proc;names(savelist[[ptype]])<-val2" to keep the information about the number of sites for each object.

## create this object once to store the data - has been run
#savelist<-list()
#save(savelist,file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/sstar_dist",sep=""))


# in every job
load(file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/sstar_dist",sep=""))
savelist[[ptype]]<-proc;names(savelist[[ptype]])<-val2
save(savelist,file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/sstar_dist",sep=""))



#save(savelist,file=paste("/scratch/devel/hpawar/admix/sstar/simul/out/sstar_dist",sep=""))




