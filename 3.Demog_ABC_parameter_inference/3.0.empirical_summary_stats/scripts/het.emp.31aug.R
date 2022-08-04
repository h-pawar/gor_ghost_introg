# Thu 26 Aug 2021 14:28:52 CEST
# try to recalculate heterozygosity for the empirical windows

# following from results of - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/generateabc.24aug.R
# disparity b/n dist of het for simulated data vs empirical data - looks like i'm calculating wrong in the empirical data 
  # where summary stats for empirical data had been calculated using - /scratch/devel/hpawar/admix/abc/simul/scripts/stats.emp.data.arr
      # which calls /scratch/devel/hpawar/admix/abc/simul/scripts/stats.emp.data.clean3.R
    # & /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.by.clbl1jul.2.R - filtering by informative windows
# explanatory notes - make.abc.model.17aug.sh
#-----------------------------------------------------------------------------------------------------------------------
# Tue 31 Aug 2021 11:59:21 CEST
# MK
#do I interpret correctly that you are calculating the heterozygosity assuming full data (40kbp) for each window? 
#We do know that the data is heavily filtered, so that may explain why it is lower than expected. 
#The best way would be the sum of all heterozygous sites divided by the sum of all sites with data in each individual.

#Also, it seems that the segregating sites are rather the number of heterozygous sites in the population, is that correct? 
#DId you also calculate it that way in the simulations? In my understanding, the segregating sites in the population 
#are all sites that are not 0/0 across all individuals and also not 1/1 across all individuals 
#(or, in terms of allele frequency >0 and <1). 
#This is not the same as being heterozygous in any given individual, so you would miss a lot of mutations that are 
#homozygous in some individuals and still segregating in the population. 
#S* does use those, so it is a relevant type of information in this case. 
#The individual heterozygosity and this value will be quite correlated (you may actually want to test how correlated they are?).
#However, if you did the same for the simulations, keep it this way.  
#In the end, this kind of information is still implicit in some of the other stats.

# me
#I will try as you suggest for heterozygosity.
#To calculate segregating sites in the simulated data I was doing length(which(nput!="00")). 
#I should also do it that way then for the empirical data, but what about the 2/2, 3/3 sites?

# MK
#I would treat it as a frequency question, i.e. calculate the frequency at each locus. 
#Or you make a statement similar to this for each line and apply as function:

#if (length(unique(inputline))==1) { uq<-unique(inputline);op=ifelse(uq%ni%c(0/0,1/1,2/2,3/3), 0,1);return(uq) } else {return(uq) }

#-----------------------------------------------------------------------------------------------------------------------
  # ie need to recalculate het & seg sites
#- sum of all heterozygous sites divided by the sum of all sites with data in each individual.
  # need to calculate sum of all sites with data per individual
#- seg sites - go from freq 
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
#------------------------------------------------------------------------------------------------------------------------
# function to calculate summary statistics per window:

het.in.win_function<-function(i) {
intvl<-c(window[[i]],window[[i+1]]) # interval to extract GTs from (40kb region)
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 1)
# count heterozygous GTs per individual & divide by window length

# not just 0s & 1s here, also alternate alleles - 
het_fun<-function(id){ 
  (length(which(win.intvl[,id]=="0|1")) + 
  length(which(win.intvl[,id]=="1|0")) +  
  length(which(win.intvl[,id]=="1|2")) +
  length(which(win.intvl[,id]=="2|1")) +
  length(which(win.intvl[,id]=="0|2")) +
  length(which(win.intvl[,id]=="2|0")) +
  length(which(win.intvl[,id]=="2|3")) +  
  length(which(win.intvl[,id]=="3|2"))) }

out.winhet<-lapply(2:ncol(win.intvl), het_fun) # apply for all 49 individuals, from col 2-50
# output this & then take mean, sd over all the segments 

# subsets of the heterozygosities - split per pop (WL, WC, EL, EM) - as nested list (of 4)
# => I should output stats in order of : 23:49, 22:22 , 13:21, 1:12
                  # WL , WC, EL, EM

#  output the sum of heterozygosity per window
# & take the mean & sds only at the end when I have all the informative windows
# MK: that is right, it should be per population, not per window. 

#out.winhet.perpop<-cbind(sum(unlist(out.winhet)[23:49]),
#  sum(unlist(out.winhet)[22]),
#  sum(unlist(out.winhet)[13:21]),
#  sum(unlist(out.winhet)[1:12]))

# output counts of het sites per id - try this
hetperid<-list(
unlist(out.winhet)[23:49],
unlist(out.winhet)[22],
unlist(out.winhet)[13:21],
unlist(out.winhet)[1:12]
  )

#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 2) segregating sites

# count if site is segregating per individual for locus x -> iterate over loci 

# output this
#seg.per.window=list()
#for (i in 1:(nrow(win.intvl))){
#seg.per.window[[i]]<-ifelse(win.intvl[i,-1]%in%c("0|0","1|1","2|2","3|3"), 0,1)
#}

# sum of segregating sites per individual for this window - may be more manageable output format than seg.per.window?
#out.seg.perwind<-colSums(data.frame(matrix(unlist(seg.per.window), nrow=length(seg.per.window), byrow=TRUE),stringsAsFactors=FALSE))
  # this doesn't make sense - equivalent to what i was obtaining from calc sum of het
#-----------------------------------------------------------------------------------------------------------------------

#MK: treat it as a frequency question, i.e. calculate the frequency at each locus. 
# #(or, in terms of allele frequency >0 and <1). 

splitgts_fun<-function(id){
do.call('rbind', strsplit(as.character(win.intvl[,id]),'|',fixed=TRUE)) 
}

win.intvl.hapl<-lapply(2:ncol(win.intvl), splitgts_fun)

y<-data.frame(matrix(unlist(win.intvl.hapl), ncol= sum(config)))
alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))

# convert from character to numeric
y[,] <- sapply(y[,],as.numeric)


# MK: You may only count non-fixed sites within each population. - Wed  1 Sep 2021 14:43:41 CEST
# MK: The easiest solution is to remove fixed sites in all gorillas. This would be:
y<-y[which(rowSums(y)<sum(config)),]


alfrfun<-function(nput) {
   inds<-which(alinds==nput)
   op<-rowSums(y[,inds])/(length(inds))
   return(op)
   }
alfre<-lapply(1:4,alfrfun)

# where 1=EM, 2=EL, 3=WC, 4=WL

#  > alfre
#[[1]]
#[1] 0.1666667 0.0000000 0.2083333

#[[2]]
#[1] 0.1111111 0.0000000 0.1666667

#[[3]]
#[1] 0 0 0

#[[4]]
#[1] 0.16666667 0.03703704 0.27777778

# if > 0 : seg site

alfre[[1]][alfre[[1]]>0] <- 1
alfre[[2]][alfre[[2]]>0] <- 1
alfre[[3]][alfre[[3]]>0] <- 1
alfre[[4]][alfre[[4]]>0] <- 1

# then take the sum  -> this is then the seg sites per pop for this window 
lapply(alfre, sum)

out.seg.perwind<-unlist(lapply(alfre, sum))
#[1] 2 2 0 3
# output this -sum of segregating sites per popn for this window 
   # note in the order of  1=EM, 2=EL, 3=WC, 4=WL (reverse to rest of stats - will need to rearrange at end)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
identifier<-cbind(chr,i, intvl[1],intvl[2]) # output identifier for this window - chr, window number, start pos, end pos

 ## output
 # ie output 1) window identifier
                # 2) count of het sites per id for this window
                # 3) count of seg sites per id for this window
 return(list(identifier, hetperid, out.seg.perwind))
}
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#het.in.win_function(1)

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

# apply het.in.win_function to windows with >1snp
	# & writing out the chr, window number & positions as well
stats.out=list()
for (j in 1:nrow(df1.fraction.wind)){
stats.out[[j]]<-het.in.win_function(df1.fraction.wind$X1[[j]])
}

# may need to find more efficient way than looping..

# write out summary stats from empirical data
save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr",chr,sep=""))


#mkdir -p /scratch/devel/hpawar/admix/abc/emp.data/test

#------------------------------------------------------------------------------------------------------------------------
# Wed  1 Sep 2021 14:36:55 CEST
# MK
#something still seems off for the segregating sites.
# It seems all populations have exactly the same number of segregating sites, which is way too high for WC, 
# too small for the easterns. 
# The reason for that is that you include fixed sites 
# (all individuals in the population are 1/1, allele frequency = 1), 
# most of which are actually fixed in comparison to humans (i.e. across the whole gorilla population).
#You may only count non-fixed sites within each population.

#alfre[[1]][alfre[[1]]==0 | alfre[[1]]==1] <- NA
#alfre[[2]][alfre[[2]]==0 | alfre[[2]]==1] <- NA
#alfre[[3]][alfre[[3]]==0 | alfre[[3]]==1] <- NA
#alfre[[4]][alfre[[4]]==0 | alfre[[4]]==1] <- NA

#alfre[[1]]<-sum(!is.na(alfre[[1]]))
#alfre[[2]]<-sum(!is.na(alfre[[2]]))
#alfre[[3]]<-sum(!is.na(alfre[[3]]))
#alfre[[4]]<-sum(!is.na(alfre[[4]]))
#out.seg.perwind<-unlist(alfre)


#Then, however, you may miss differences between the populations, and I dont 
#remember how you calculated the segregating sites in the simulations. 
#The easiest solution is to remove fixed sites in all gorillas. This would be (instead of the code above):

#y[,] <- sapply(y[,],as.numeric)
#y<-y[which(rowSums(y)<sum(config)),]

#I was just checking a bit and it seems using the population-based number 
#or the dataset-based number are quite similar. 
#In a random window, 54% of sites are 1/1 across all individuals, 
#so this is the baseline that shifted your segsites so much. 

#Then, again, I dont know how you calculated the segsites in the simulations, 
#because the WC value might still be too high, it shares a lot of 1/1 with WL.
#I would suggest first trying the second strategy (changing "y") and comparing it to the simulated data.
#For the heterozygosity, WL are a bit lower than expected, and I think that might be due to the filtering, 
#removing more data than necessary, potentially removing more often the hets than homozygous sites. 
#That will lead to an underestimate of the Ne (in all populations), but we have to accept that.
#------------------------------------------------------------------------------------------------------------------------
