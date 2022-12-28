#module load R/4.0.1
#!/usr/bin/r

#-----------------------------------------------------------------------------------------------------------------------
# Wed  2 Feb 2022 16:01:55 CET
# MK
#I think the main issue is the way to read in data, change column/row:
#md_geno <- matrix(scan(file="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_22_Egor.biallelic.GT"), ncol=nind, nrow=nsnps, byrow=T)
#Also, you don't need to do a processing step for missing data, as this table only uses complete genotypes (like S*). 
#Hence, obs.sfs<-table(rowSums(md_geno))
#-----------------------------------------------------------------------------------------------------------------------
require(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
a=args[1] # pass from cd line - chr number
#-----------------------------------------------------------------------------------------------------------------------
# nsnps per GT
#for i in {1..22}
#do
#   cat  4_"$i"_Egor.biallelic.filt.GT | wc -l
#done

# number of SNPs per GT file
#nrows<-c(8649941,10590950,8801570,8064032,7794661,7592642,6306813,6322633,1325944,5550879,5110818,5372361,4328829,3658387,3259675,2671532,2406836,3586343,946585,2275219,1367710,782934)

#sum(nrows)
#106767294

# cat GT files per chr in order -> autosomes GT file
#cat  4_{1..22}_Egor.biallelic.filt.GT >> tmp/4_autosomes_Egor.biallelic.filt.GT

#cat tmp/4_autosomes_Egor.biallelic.filt.GT | wc -l
#106767294
#-----------------------------------------------------------------------------------------------------------------------

# number of SNPs
#nsnps <- 106767294
# number of individuals
nind <- 21

# read the genotype matrix for all autosomes
#md_geno <- matrix(scan(file="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_Egor.biallelic.filt.GT"), ncol=nind, nrow=nsnps, byrow=T)

#srun: Exceeded job memory limit
#srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
#srun: error: cnb4: task 0: Killed
#Session closed
# send as a job
    # still exceeding job memory ->

#-----------------------------------------------------------------------------------------------------------------------

# number of SNPs per GT file
nrows<-c(8649941,10590950,8801570,8064032,7794661,7592642,6306813,6322633,1325944,5550879,5110818,5372361,4328829,3658387,3259675,2671532,2406836,3586343,946585,2275219,1367710,782934)

#sum(nrows)
#106767294


input<-paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_",a,"_Egor.biallelic.filt.GT")
searchString <- ' '
replacementString <- ''
input = gsub(searchString,replacementString,input)
# read in the GT matrix for this chr
md_geno <- matrix(scan(file=input), ncol=nind, nrow=nrows[[a]], byrow=T)


obs.sfs<-table(rowSums(md_geno))

obs.sfs # should be 0 -42

# then only normalise once have calculated site categories for all chromosomes

# only once
savelist<-list()
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))


# in every job
load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))
savelist[[a]]<-obs.sfs 
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))


q()
#-----------------------------------------------------------------------------------------------------------------------

# Run the following interactively - 

# write into a function & run interactively
countsites<-function(a){

nind <- 21
# number of SNPs per GT file
nrows<-c(8649941,10590950,8801570,8064032,7794661,7592642,6306813,6322633,1325944,5550879,5110818,5372361,4328829,3658387,3259675,2671532,2406836,3586343,946585,2275219,1367710,782934)

input<-paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_",a,"_Egor.biallelic.filt.GT")
searchString <- ' '
replacementString <- ''
input = gsub(searchString,replacementString,input)
# read in the GT matrix for this chr
md_geno <- matrix(scan(file=input), ncol=nind, nrow=nrows[[a]], byrow=T)

obs.sfs<-table(rowSums(md_geno))
obs.sfs # should be 0 -42

# in every job
load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))
savelist[[a]]<-obs.sfs 
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))

}

countsites(3)

load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""),verbose=T)
savelist # check chr 3 also added - yes worked

for (i in 4:22){
countsites(i)
}    

savelist
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))

# has generated output
# now concatenate all together

#-----------------------------------------------------------------------------------------------------------------------
# should unlist & sum
all<-data.frame(matrix(unlist(savelist), ncol=length(savelist[[1]]), byrow=TRUE),stringsAsFactors=FALSE)

# looks fine - comparable to head(savelist)
 #head(test)
#       X1    X2    X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14
#1 7858089 18443  8343 6947 5972 5295 5462 5957 4967 4709 4558 4146 4213 4062
#2 9628209 22289 10003 7478 7340 5959 6319 6802 5615 6000 5574 5014 5003 5015
#3 8008825 17385  7133 6342 6306 5181 4981 4788 4482 4380 4717 4676 4430 4029

aut.obs.sfs <-colSums(all)

aut.obs.sfs
      X1       X2       X3       X4       X5       X6       X7       X8 
96904188   223702    98017    84862    75737    68045    67045    67220 
      X9      X10      X11      X12      X13      X14      X15      X16 
   61663    60811    59081    57057    56211    54495    55390    54200 
     X17      X18      X19      X20      X21      X22      X23      X24 
   50920    49796    52704    52563    48313    78569    42659    43083 
     X25      X26      X27      X28      X29      X30      X31      X32 
   41525    37105    36454    36669    35321    33104    32689    31330 
     X33      X34      X35      X36      X37      X38      X39      X40 
   31112    30899    29576    30252    29544    27156    28482    29309 
     X41      X42      X43 
   30003    61130  7689303 

sum(aut.obs.sfs)
[1] 106767294

# then only normalise once have calculated site categories for all chromosomes
# normalise 
test<-aut.obs.sfs/sum(aut.obs.sfs)

names(test)<-c(0, seq(1:42))

#sum(test)
#[1] 1

# drop first category
test1<-test[-1]

x<-as.data.table(test1)

x$V1<-as.numeric(names(test1))

write.table(x,file="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfs.input.txt",sep="\t",row.names=FALSE,col.names=FALSE)

#-----------------------------------------------------------------------------------------------------------------------

#load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""),verbose=T)
#Loading objects:
#  savelist

# then only normalise once have calculated site categories for all chromosomes
# normalise 
#test<-obs.sfs/sum(obs.sfs)

#sum(test)
#[1] 1

# drop first category
#test1<-test[-1]

#x<-as.data.table(test1)

#x$V1<-as.numeric(x$V1)

#write.table(x,file="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfs.input.txt",sep="\t",row.names=FALSE,col.names=FALSE)

#-----------------------------------------------------------------------------------------------------------------------
# this approach didnt work

# number of individuals
#nind <- 21
#auto_geno<-c()
#for (i in 1:length(nrows)){
#input<-paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_",i,"_Egor.biallelic.filt.GT")
#searchString <- ' '
#replacementString <- ''
#input = gsub(searchString,replacementString,input)
#auto_geno[[i]] <- matrix(scan(file=input), ncol=nind, nrow=nrows[[i]], byrow=T)
#}

#Read 181648761 items
#srun: Exceeded job memory limit
#srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
#srun: error: cnb4: task 0: Killed
#Session closed

# perhaps send as a job
# or perhaps cat first in bash

#106,767,294 # should be total nsnps
#181,648,761
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# choosing how many blocks to split each chr into **
# R
require(data.table)
options(stringsAsFactors=F)
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))


# hg19 chr in the correct order
hg19[c(1:18,20,19,22,21),]
#sl<-unlist(hg19[c(1:18,20,19,22,21),2])
#names(sl)<-unlist(hg19[c(1:18,20,19,22,21),1])

# how mnay blocks eg if chr 22 split into 10 runs quickly, how many shoudl the rest of the chr be split into
 51304566/10
[1] 5130457 # per block consists of 5,130kb

hg19<-hg19[c(1:18,20,19,22,21),]
hg19$V3<-hg19[,2]/51304566

hg19
       V1        V2        V3
 1:  chr1 249250621 4.8582542   # 50 blocks
 2:  chr2 243199373 4.7403066   # 50 blocks
 3:  chr3 198022430 3.8597428   # 40 blocks
 4:  chr4 191154276 3.7258726   # 40 blocks
 5:  chr5 180915260 3.5262994   # 35 blocks
 6:  chr6 171115067 3.3352795   # 35 blocks
 7:  chr7 159138663 3.1018421   # 30 blocks
 8:  chr8 146364022 2.8528459   # 30 blocks
 9:  chr9 141213431 2.7524535   # 30 blocks
10: chr10 135534747 2.6417677   # 25 blocks
11: chr11 135006516 2.6314717   # 25 blocks
12: chr12 133851895 2.6089665   # 25 blocks
13: chr13 115169878 2.2448271   # 20 blocks
14: chr14 107349540 2.0923974   # 20 blocks
15: chr15 102531392 1.9984847   # 20 blocks
16: chr16  90354753 1.7611445   # 20 blocks
17: chr17  81195210 1.5826118   # 15 blocks
18: chr18  78077248 1.5218382   # 15 blocks
19: chr19  59128983 1.1525092   # 10 blocks
20: chr20  63025520 1.2284583   # 10 blocks
21: chr21  48129895 0.9381211   # 10 blocks
22: chr22  51304566 1.0000000   # 10 blocks
       V1        V2        V3

# ie most should be 10 blocks, but could increase to 50/40 for the larger chr? # think that makes sense
#-----------------------------------------------------------------------------------------------------------------------
