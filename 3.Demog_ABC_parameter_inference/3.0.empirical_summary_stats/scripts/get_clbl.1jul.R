#!/usr/bin/r
#module load TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools java/latest
# adapting get_clbl.MK.R MK script -  
  # The whole genome vcf files are here: /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz. However, I would recommend to take windows with data after filtering for repeats and mapability.
  # I have a script that kind of does that for 1kbp windows (for the Skov method), see attached. Maybe you can try to adjust it for your purpose, since I don't think I will have time to do it very soon?
#-----------------------------------------------------------------------------------------------------------------------

print(Sys.time())
options(scipen=100)
args = commandArgs(trailingOnly=TRUE)
#chrom=(commandArgs(TRUE))
#chrom=22 # testing how this script works with chr 22
#chrom=21
chrom=args[1]
chrom=as.character(chrom)
print(chrom)
chrom<-unlist(strsplit(chrom,split="_"))
print(chrom)
'%ni%' <- Negate('%in%')
step=1

gors<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/gorilla.lst",sep="\t",as.is=T))


#str(gors)
# Named chr [1:49] "Gorilla_beringei_beringei-Bwiruka" ...
# - attr(*, "names")= chr [1:49] "V11" "V12" "V13" "V14" ...


grps<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/groups.lst",sep="\t",as.is=T))

#head(grps)
#                                          V11 
#"/scratch/devel/mkuhlwilm/arch/gorilla_1.lst" 
#                                          V12 
#"/scratch/devel/mkuhlwilm/arch/gorilla_2.lst" 
#                                          V13 

grn<-do.call(rbind,strsplit(do.call(rbind,strsplit(grps,split="/"))[,6],split="\\."))[,1]
inds<-list();finds<-list()
for ( i in (1:length(grps))) { inds[[i]]<-unlist(read.table(grps[i],sep="\t",as.is=T)); inds[[i]]<-inds[[i]][which(inds[[i]]%ni%c("Gorilla_beringei_graueri-Serufuli","Pan_troglodytes_schweinfurthii-Mgbadolite"))]   }


# inds[[1]] -> [[4]] = gorillas split per population
# inds[[1]]
#                                 V11                                  V12 
# "Gorilla_beringei_beringei-Bwiruka"   "Gorilla_beringei_beringei-Imfura" 
#                                 V13                                  V14 
#  "Gorilla_beringei_beringei-Kaboko" "Gorilla_beringei_beringei-Kahungye" 

# inds[[3]]
#                                  V1 
#"Gorilla_gorilla_diehli-B646_Nyango" 


#inds[[5]] # CR + WL
#inds[[6]] # EM + EL
  #length(inds[[6]])+length(inds[[5]])
#[1] 49     

indlist<-c(1:3,unlist(read.table("/scratch/devel/mkuhlwilm/ga/findivs.txt",header=F,as.is=T,sep="\t")))
# individuals for all species

for (i in 1:length(grps)) { finds[[i]]<-which(indlist%in%inds[[i]]) }
fugrp<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/species.lst",sep="\t",as.is=T))
fugrp<-cbind(c(paste("gorilla_",1:6,sep=""),paste("pan_",1:6,sep=""),paste("hum_",1:2,sep=""),paste("ora_",1:4,sep=""),paste("hupa",1:2,sep="")),c(rep(fugrp[1],6),rep(fugrp[2],6),rep(fugrp[3],2),rep(fugrp[4],4),rep(fugrp[5],2)))


#fugrp
#    [,1]        [,2]                                       
#V11 "gorilla_1" "/scratch/devel/mkuhlwilm/arch/gorilla.lst"
#V11 "gorilla_2" "/scratch/devel/mkuhlwilm/arch/gorilla.lst"
#V11 "gorilla_3" "/scratch/devel/mkuhlwilm/arch/gorilla.lst"
#V11 "gorilla_4" "/scratch/devel/mkuhlwilm/arch/gorilla.lst"
#V11 "gorilla_5" "/scratch/devel/mkuhlwilm/arch/gorilla.lst"
#V11 "gorilla_6" "/scratch/devel/mkuhlwilm/arch/gorilla.lst"


minds<-list();funds<-list()
for ( i in (1:length(grps))) { minds[[i]]<-unlist(read.table(fugrp[i,2],sep="\t",as.is=T)); minds[[i]]<-minds[[i]][which(minds[[i]]%ni%c("Gorilla_beringei_graueri-Serufuli","Pan_troglodytes_schweinfurthii-Mgbadolite"))]   }
for (i in 1:length(grps)) { funds[[i]]<-which(indlist%in%minds[[i]] & indlist%ni%inds[[i]]) }
atyp=unique(fugrp[,2])
#write.table(atyp,file="/scratch/devel/mkuhlwilm/arch/agroups.lst",sep="\t",row.names=F,col.names=F,quote=F)  


# minds[1] # the gorillas
#$ : Named chr [1:49] "Gorilla_beringei_beringei-Bwiruka" "Gorilla_beringei_beringei-Imfura" "Gorilla_beringei_beringei-Kaboko" "Gorilla_beringei_beringei-Kahungye" ...
#  ..- attr(*, "names")= chr [1:49] "V11" "V12" "V13" "V14" ...

#head(minds[[1]])
#                                 V11                                  V12 
# "Gorilla_beringei_beringei-Bwiruka"   "Gorilla_beringei_beringei-Imfura" 
#                                 V13                                  V14 
#  "Gorilla_beringei_beringei-Kaboko" "Gorilla_beringei_beringei-Kahungye" 
#                                 V15                                  V16 
# "Gorilla_beringei_beringei-Katungi"   "Gorilla_beringei_beringei-Maisha" 



library("GenomicRanges")
#chilen<-read.table("/home/devel/mkuhlwilm/hg19.chrom.sizes",sep="\t",header=F,nrow=24)[-c(21),]

#Error in file(file, "rt") : cannot open the connection
#In addition: Warning message:
#In file(file, "rt") :
#  cannot open file '/home/devel/mkuhlwilm/hg19.chrom.sizes': Permission denied

# use instead values here
#https://github.com/igvteam/igv/blob/master/genomes/sizes/hg19.chrom.sizes

require(data.table)
# read in lengths of each chr
chilen<-data.table(
#V1=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"),
V1=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
V2= c(249250621,  243199373,  198022430,  191154276,  180915260,  171115067,  159138663,  146364022,  141213431,  135534747,  135006516,  133851895,  115169878,  107349540,  102531392,  90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566) )

# maybe read in as numeric 1-22, instead of with chr?

#chilen[,1]<-unlist(do.call(rbind,strsplit(as.character(chilen[,1]),split="hr"))[,2])
#chilen<-chilen[,c(1:2)];cv<-seq(0,250000000,1000000)

#cv<-seq(0,250000000,1000000) # 1,000,000 = 1kb windows - change this here to 40kb - whether i have been doing htis wrong in the emp data?

# in simulations & the empirical data I have been using:
#len=40000 # for the empirical data shoudl this have been 40,000,000?
  # checkign this w MK **

#-----------------------------------------------------------------------------------------------------------------------

#can I also check for your script, you take intervals as seq(0,250000000,1000000).
#For the empirical data's summary statistics I was taking window<-seq(min(pos), max(pos), 40000), since I was using len=40000  for the simulations. I should instead use seq(min(pos), max(pos), 40000000) right? '

# MK
#sorry, the code is a bit raw and not well enough commented. These intervals are not correct 
#(in fact, not used in the script, a leftover). The way I defined the 1kbp windows is as follows:
#tevals<-chilen[which(chilen[,1]==chr),2]; tevals<-cbind(seq(0,tevals,1000));  tevals<-cbind(tevals,tevals+999)
#This may be modified for non-overlapping 40 kbp windows.
#-----------------------------------------------------------------------------------------------------------------------

#tevals<-chilen[which(chilen[,1]==chrom),2]; tevals<-cbind(seq(0,tevals,1000));  tevals<-cbind(tevals,tevals+999)
#tevals<-chilen[which(chilen[,1]==chrom),2]; tevals<-cbind(seq(0,tevals[[1]],40000));  tevals<-cbind(tevals,tevals+39999)
#tevals<-IRanges(start=tevals[,1],end=tevals[,2])

#> tevals
#         V2
#1: 51304566 # size of specific chr
# head(tevals)
#       [,1]   [,2]
#[1,]      0   3999
#[2,]  40000  43999
#[3,]  80000  83999
#[4,] 120000 123999
#[5,] 160000 163999
#[6,] 200000 203999

# head(tevals)
#IRanges of length 6
#     start    end width
#[1]      0   3999  4000
#[2]  40000  43999  4000
#[3]  80000  83999  4000
#[4] 120000 123999  4000
#[5] 160000 163999  4000
#[6] 200000 203999  4000


  
############################################################
## for each species, get positions that have data, and that fulfill the two filtering criteria
#tfun<-function(intv, ag) { c(intv,sum(ranges(findOverlaps(IRanges(start=intv[1],end=intv[2]),ag),IRanges(start=intv[1],end=intv[2]),ag)@width)/1000) }
#the output of the function (object "avals") from my script should give straightforward the proportion of each window that has data. 
# in order to get the proportion should divide by 40,000 instead of by 1000
tfun<-function(intv, ag) { c(intv,sum(ranges(findOverlaps(IRanges(start=intv[1],end=intv[2]),ag),IRanges(start=intv[1],end=intv[2]),ag)@width)/40000) }
newchr<-c(1:22,"X")


yy=as.numeric(chrom[1])
#typ=atyp[yy] # this doesn't work here

# atyp
#[1] "/scratch/devel/mkuhlwilm/arch/gorilla.lst"
#[2] "/scratch/devel/mkuhlwilm/arch/pans.lst"   
#[3] "/scratch/devel/mkuhlwilm/arch/human.lst"  
#[4] "/scratch/devel/mkuhlwilm/arch/orang.lst"  
#[5] "/scratch/devel/mkuhlwilm/arch/hupa.lst"   

# typ=atyp[yy] # this line doesn't work
#> typ
#[1] NA


# cat /scratch/devel/mkuhlwilm/arch/gorilla.lst | wc -l
#49
#[hpawar@cnb2 ~]$ head /scratch/devel/mkuhlwilm/arch/gorilla.lst  # the gorilla individuals
#Gorilla_beringei_beringei-Bwiruka
#Gorilla_beringei_beringei-Imfura
#Gorilla_beringei_beringei-Kaboko
#Gorilla_beringei_beringei-Kahungye



# i don't have fields 2 & 3 of chrom?
#ft=as.numeric(chrom[2])
#chr=as.character(chrom[3])
#if (chr==23) { chr<-"X" }


typ=atyp[1]

spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
#spec
#[1] "gorilla"


#print(c(typ,ft,chr,spec))

#    if(length(nchr)==0) { next }
#    for (chr in (nchr)) {
  # whether this first one that is commented out also need to be run?
  #  agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -S ",typ," /scratch/devel/mkuhlwilm/gvcfs/greatape_",chrom,".vcf.gz | /apps/BCFTOOLS/1.9/bin/bcftools view --exclude 'F_MISSING>0.9' -G | bedtools merge -i stdin",sep=""),intern=T)


# MK
#2)  ft == 1 is for least filtering, ft == 2 for more filtering.
#I recommend using the command from ft == 2, since this is the data we actually used downstream of this.

ft=2
chr=chrom
#if (ft==2) {    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -S ",typ," /scratch/devel/mkuhlwilm/gvcfs/greatapeN_",chr,".vcf.gz -G | intersectBed -wa -v -header -a stdin -b /project/devel/mkuhlwilm/RM.bed | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep \"#|mapability=1\" | bedtools merge -i stdin",sep=""),intern=T) }
#/scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig: line 2: java: command not found
#Error: Unable to open file /project/devel/mkuhlwilm/RM.bed. Exiting.


#the /project/devel/mkuhlwilm/RM.bed file no longer exists? So I need to regenerate it.
#My question is, to run repeatmasker for the gorillas can I use the default options, or did you specify different flags to make this bed file?

# MK: Project folder is invalid. Try adding /scratch/ in front. RM is for human, precomputed.

#/scratch/project/devel/mkuhlwilm/RM.bed # adding /scratch/ in front works

# head /scratch/project/devel/mkuhlwilm/RM.bed
#chr1  16777160  16777470  AluSp 2147  +
#chr1  25165800  25166089  AluY  2626  -
#chr1  33553606  33554646  L2b 626 +
#chr1  50330063  50332153  L1PA10  12545 +
#chr1  58720067  58720973  L1PA2 8050  -
#chr1  75496180  75498100  L1MB7 10586 +
#chr1  83886030  83886750  ERVL-E-int  980 -

if (ft==2) {    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -S ",typ,"  /scratch/devel/mkuhlwilm/gvcfs/greatapeN_",chr,".vcf.gz -G |   intersectBed -wa -v -header -a stdin -b /scratch/project/devel/mkuhlwilm/RM.bed |   /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep \"#|mapability=1\" | bedtools merge -i stdin",sep=""),intern=T) }
# this is running!
#[main] INFO vcfbigwig - Starting JOB at Thu Jul 01 17:50:07 CEST 2021 com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig version=a0e66b2b00d9c1a58047eb9cf2eec26bc5efd4fb  built=2016-09-19:14-09-18
#[main] INFO vcfbigwig - Command Line args : -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/
#[main] INFO vcfbigwig - Executing as hpawar@cnb5 on Linux 2.6.32-696.13.2.el6.Bull.128.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0_261-b12

# bcftools view -S : extract sam
#-S, --samples-file FILE
#see Common Options. Note that it is possible to create multiple subsets simultaneously using the split plugin.


#[main] INFO vcfbigwig - Count: 18325648 Elapsed: 2 hours(89.19%) Remains: 16 minutes(10.81%) Last: chr21:48026998
#[main] INFO vcfbigwig - Count: 18341141 Elapsed: 2 hours(89.19%) Remains: 16 minutes(10.81%) Last: chr21:48073955
#[main] INFO vcfbigwig - done: N=18359758
#[main] INFO vcfbigwig - End JOB  [Fri Jul 02 12:01:42 CEST 2021] VCFBigWig done. Elapsed time: 139.60 minutes

#head(agt)
#[1] "chr21\t9411193\t9411194" "chr21\t9411932\t9411967"
#[3] "chr21\t9411995\t9412030" "chr21\t9412061\t9412100"
#[5] "chr21\t9412110\t9412158" "chr21\t9412205\t9412240"

    #if (ft==1) {     agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -S ",typ," /scratch/devel/mkuhlwilm/gvcfs/greatapeN_",chr,".vcf.gz | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep \"#|mapability=1\" | bedtools merge -i stdin",sep=""),intern=T) }
    print("data loaded")
    agt<-do.call(rbind,strsplit(agt,split="\t"))
  

# head(agt)
#     [,1]    [,2]      [,3]     
#[1,] "chr21" "9411193" "9411194"
#[2,] "chr21" "9411932" "9411967"
#[3,] "chr21" "9411995" "9412030"
#[4,] "chr21" "9412061" "9412100"
#[5,] "chr21" "9412110" "9412158"
#[6,] "chr21" "9412205" "9412240"

# str(agt)
# chr [1:60344, 1:3] "chr21" "chr21" "chr21" "chr21" "chr21" ...
#>

    agtr<-IRanges(start=as.numeric(agt[,2]),end=as.numeric(agt[,3]))
    #tevals<-chilen[which(chilen[,1]==chr),2]; tevals<-cbind(seq(0,tevals,1000));  tevals<-cbind(tevals,tevals+999)

# make tevals & do not convert to iranges object
tevals<-chilen[which(chilen[,1]==chrom),2]; tevals<-cbind(seq(0,tevals[[1]],40000));  tevals<-cbind(tevals,tevals+39999)
   
#head(agtr)
#IRanges of length 6
#      start     end width
#[1] 9411193 9411194     2
#[2] 9411932 9411967    36
#[3] 9411995 9412030    36
#[4] 9412061 9412100    40
#[5] 9412110 9412158    49
#[6] 9412205 9412240    36

#------------------------------------------------------------------------------------------------------------------------

# MK: 

# the "tevals" object should not be an IRanges object. 
#The function is designed to work iteratively on a matrix/data.frame with the standard apply() command, taking one line at a time
# and making a small IRanges object out of it when applying the function. 

#the output of the function (object "avals") from my script should give straightforward the proportion of each window that has data. 
#You should keep that information stored for all windows and chromosomes, because you can select later on but have all information. 
#You don't want to know if there are overlaps (which is what findOverlaps is doing), but the sum of the overlap.

#------------------------------------------------------------------------------------------------------------------------

avals<-apply(tevals,1,tfun,ag=agtr)
avals<-t(avals)

#head(avals)
#       [,1]   [,2] [,3]
#[1,]      0  39999    0
#[2,]  40000  79999    0
#[3,]  80000 119999    0
#[4,] 120000 159999    0
#[5,] 160000 199999    0
#[6,] 200000 239999    0


#tail(avals)
#            [,1]     [,2]   [,3]
#[1199,] 47920000 47959999 15.024
#[1200,] 47960000 47999999 26.569
#[1201,] 48000000 48039999 11.799
#[1202,] 48040000 48079999 16.963
#[1203,] 48080000 48119999 11.755
#[1204,] 48120000 48159999  0.000


    aweight1<-as.data.frame(cbind(chr,avals[,c(1)]))
    aweight2<-as.character(as.double(avals[,3]))
    aweight3<-cbind(aweight1,format(aweight2, digits=3, scientific=F))
  

#head(aweight3)
#  chr     V2 format(aweight2, digits = 3, scientific = F)
#1  21      0                                       0     
#2  21  40000                                       0     
#3  21  80000                                       0     
#4  21 120000                                       0     
#5  21 160000                                       0     
#6  21 200000                                       0     
#> tail(aweight3)
#     chr       V2 format(aweight2, digits = 3, scientific = F)
#1199  21 47920000                                       15.024
#1200  21 47960000                                       26.569
#1201  21 48000000                                       11.799
#1202  21 48040000                                       16.963
#1203  21 48080000                                       11.755
#1204  21 48120000                                       0    


  #  load(file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
  #  alweig[[chr]]<-aweight3 # this is looping through all chr?
      # doesnt make sense b/c have only read in chr 21
        # if send as array job, may be best to write out as sep files?
  #  save(alweig,file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  

# write out per chr
save(aweight3,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/",spec,"_weight",ft,".",chrom,".Robject",sep=""))

# mkdir -p /scratch/devel/hpawar/admix/abc/emp.data/window.info

q()
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

 #   avals<-apply(tevals,1,tfun,ag=agtr) # error here **
      # error in applying the tfun - instead of as function, do step by step
   
#Error in apply(tevals, 1, tfun, ag = agtr) : 
#  dim(X) must have a positive length
# which are the x which is not having a positve length?
  # read the forums re this error message
#tfun
#function(intv, ag) { c(intv,sum(ranges(findOverlaps(IRanges(start=intv[1],end=intv[2]),ag),IRanges(start=intv[1],end=intv[2]),ag)@width)/1000) }


#avals<-lapply(tevals,1,tfun,ag=agtr)
#Error in match.fun(FUN) : '1' is not a function, character or symbol

# ok, need to figure this out
# maybe the function doesn't work for this data?
  # may need to rewrite to findoverlaps



  tevals
#IRanges of length 1204
#          start      end width
#[1]           0     3999  4000
#[2]       40000    43999  4000
#[3]       80000    83999  4000


agtr
#IRanges of length 60344
#           start      end width
#[1]      9411193  9411194     2
#[2]      9411932  9411967    36
#[3]      9411995  9412030    36

# changing from apply -> lapply/sapply doesnt help here
#avals<-lapply(tevals,1,tfun,ag=agtr)
#Error in match.fun(FUN) : '1' is not a function, character or symbol
#>  avals<-sapply(tevals,1,tfun,ag=agtr)
#Error in match.fun(FUN) : '1' is not a function, character or symbol


#avals<-apply(tevals,tfun,ag=agtr)
#Error in match.fun(FUN) : argument "FUN" is missing, with no default
#> avals<-apply(tevals,0,tfun,ag=agtr)
#Error in apply(tevals, 0, tfun, ag = agtr) : 
#  dim(X) must have a positive length

# check dim of both

# both have null dimensions - typical of iranges objects
#> dim(tevals)
#NULL
#> dim(agtr)
#NULL

# both null

#str(tevals)
#Formal class 'IRanges' [package "IRanges"] with 6 slots
#  ..@ start          : int [1:1204] 0 40000 80000 120000 160000 200000 240000 280000 320000 360000 ...
#  ..@ width          : int [1:1204] 4000 4000 4000 4000 4000 4000 4000 4000 4000 4000 ...
#  ..@ NAMES          : NULL
#  ..@ elementType    : chr "integer"
#  ..@ elementMetadata: NULL
#  ..@ metadata       : list()

#-----------------------------------------------------------------------------------------------------------------------

#the apply function takes only entire data tables as input.

# read into overlaps b/n the 2 
# keep in mind that the apply function can not be used for one-dimensional data.

 # perhaps use the function, for only the first interval?
#tfun(tevals[1],agtr)

#Error in ranges(findOverlaps(IRanges(start = intv[1], end = intv[2]),  : 
#  error in evaluating the argument 'x' in selecting a method for function 'ranges': Error in findOverlaps(IRanges(start = intv[1], end = intv[2]), ag) : 
 # error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': Error in NSBS(i, x, exact = exact, upperBoundIsStrict = !allow.append) : 
 # subscript contains NAs or out-of-bounds indices


#tfun(tevals[1],ag=agtr)
#Error in ranges(findOverlaps(IRanges(start = intv[1], end = intv[2]),  : 
#  error in evaluating the argument 'x' in selecting a method for function 'ranges': Error in findOverlaps(IRanges(start = intv[1], end = intv[2]), ag) : 
#  error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': Error in NSBS(i, x, exact = exact, upperBoundIsStrict = !allow.append) : 
#  subscript contains NAs or out-of-bounds indices

#tfun(tevals,ag=agtr)
#Error in ranges(findOverlaps(IRanges(start = intv[1], end = intv[2]),  : 
#  error in evaluating the argument 'x' in selecting a method for function 'ranges': Error in findOverlaps(IRanges(start = intv[1], end = intv[2]), ag) : 
#  error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': Error in IRanges(start = intv[1], end = intv[2]) : 
#  'end' and 'width' must be NULLs when 'start' is a Ranges object

#-----------------------------------------------------------------------------------------------------------------------

# need to read into findoverlaps & iranges **

#http://www.biostat.jhsph.edu/~khansen/IRangesLecture.pdf
#Each integer only occur in a single range and there are as few ranges as possible. In addition, it is ordered.


findOverlaps(tevals,agtr)
Hits object with 6436 hits and 0 metadata columns:
         queryHits subjectHits
         <integer>   <integer>
     [1]       237         171
     [2]       237         172
     [3]       237         173
     [4]       237         174
     [5]       237         175
     ...       ...         ...
  [6432]      1202       60232
  [6433]      1203       60278
  [6434]      1203       60279
  [6435]      1203       60280
  [6436]      1203       60281
  -------
  queryLength: 1204
  subjectLength: 60344


# whether i want both hits in the query & the subject? 

# need to go from here back to the positions
  ov<-findOverlaps(tevals,agtr)

tevals[ov@queryHits,]
IRanges of length 6436
          start      end width
[1]     9440000  9443999  4000
[2]     9440000  9443999  4000
[3]     9440000  9443999  4000
[4]     9440000  9443999  4000
[5]     9440000  9443999  4000
...         ...      ...   ...
[6432] 48040000 48043999  4000
[6433] 48080000 48083999  4000
[6434] 48080000 48083999  4000
[6435] 48080000 48083999  4000
[6436] 48080000 48083999  4000
# is this the desired output?
  # & filter the summary stats windows by these positions?


#overlapsRanges(query,subject,hits) 
#overlapsRanges(tevals,agtr,ov) 
#Error: could not find function "overlapsRanges"


#Compute percent overlap and filter the hits:

overlaps <- pintersect(tevals[queryHits(ov)], agtr[subjectHits(ov)])

 overlaps
IRanges of length 6436
          start      end width
[1]     9440000  9440016    17
[2]     9440021  9440062    42
[3]     9440127  9440195    69
[4]     9440201  9440253    53
[5]     9440258  9440266     9
...         ...      ...   ...
[6432] 48043601 48043999   399
[6433] 48080000 48081856  1857
[6434] 48082021 48082672   652
[6435] 48082756 48083558   803
[6436] 48083741 48083999   259



percentOverlap <- width(overlaps) / width(agtr[subjectHits(ov)]) # whether am dividing by the right denominator here?
hits <- ov[percentOverlap > 0.5]


hits
Hits object with 5980 hits and 0 metadata columns:
         queryHits subjectHits
         <integer>   <integer>
     [1]       237         172
     [2]       237         173
     [3]       237         174
     [4]       237         175
     [5]       237         176
     ...       ...         ...
  [5976]      1202       60230
  [5977]      1202       60231
  [5978]      1203       60278
  [5979]      1203       60279
  [5980]      1203       60280
  -------
  queryLength: 1204
  subjectLength: 60344

# is this step necessary?  
tevals[hits@queryHits,]
tevals[hits@queryHits,]
IRanges of length 5980
          start      end width
[1]     9440000  9443999  4000
[2]     9440000  9443999  4000
[3]     9440000  9443999  4000
[4]     9440000  9443999  4000
[5]     9440000  9443999  4000
...         ...      ...   ...
[5976] 48040000 48043999  4000
[5977] 48040000 48043999  4000
[5978] 48080000 48083999  4000
[5979] 48080000 48083999  4000
[5980] 48080000 48083999  4000
# or output this, after have filtered by percentoverlap?


# Q - whether to use intersect() instead of findoverlaps? 
  # here fewer regions picked up than in findoverlaps - presumably here only the unique regions?
#intersect" them, finding the set of regions formed by the overlaps 
  intersect(tevals,agtr)
IRanges of length 6266
          start      end width
[1]     9440000  9440016    17
[2]     9440021  9440062    42
[3]     9440127  9440195    69
[4]     9440201  9440253    53
[5]     9440258  9440266     9
...         ...      ...   ...
[6262] 48043601 48043999   399
[6263] 48080000 48081856  1857
[6264] 48082021 48082672   652
[6265] 48082756 48083558   803
[6266] 48083741 48083999   259


#-----------------------------------------------------------------------------------------------------------------------


# try this - https://www.biostars.org/p/95062/#95098
# now find the overlaps (if doesnt work will need to compile qs)
#o <- findOverlapsdf2.ir, df1.tree, type = "within")
#df2[o@queryHits,]

# https://bioconductor.org/packages/devel/bioc/manuals/IRanges/man/IRanges.pdf
# coudl try this as well:
# overlapsRanges(query,subject,hits) returns a IntegerRanges (or IntegerRangesList) object of the 
# same shape as hits holding the regions of intersec- tion between the overlapping ranges in 
# objects query and subject, which should be the same query and subject used in the call to findOverlaps that generated hits. 
# Same shape means same length when hits is a Hits object, and same length and same elementNROWS when hits is a HitsList object.

# or try this?
#mergeByOverlaps computes the overlap between query and subject according to the arguments in .... 
#It then extracts the corresponding hits from each object and returns a DataFrame containing one column for the 
#query and one for the subject, as well as any mcols that were present on either object. 
#The query and subject columns are named by quoting and deparsing the corresponding argument.

#-----------------------------------------------------------------------------------------------------------------------

# whether to filter by percent overlap
#https://support.bioconductor.org/p/72656/
#refGR <- makeGRangesFromDataFrame(REF)
#testGR <- makeGRangesFromDataFrame(TEST)
#Find overlaps:

#hits <- findOverlaps(refGR, testGR)
#Compute percent overlap and filter the hits:

#overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
#percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
#hits <- hits[percentOverlap > 0.5]
# try this 
#-----------------------------------------------------------------------------------------------------------------------


#ranges(ov)
#Error in is(query, "Ranges") : 
 # argument "query" is missing, with no default

# would intersect give the positons?

# tutorial re iranges
#https://kasperdanielhansen.github.io/genbioconductor/html/IRanges_Basic.html
# hits =  two-column matrix of indicies into the two IRanges.
#The elements of unique(queryHits) gives you the indices of the query ranges which actually had an overlap; you need unique because a query range may overlap multiple subject ranges.

# further tutorial
# http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/plyrangesWorkshops__common-tasks/

# Find overlaps
#over <- foverlaps(groupA, groupB, nomatch = 0)
# Extract exact regions
#over2 <- data.table(
 #   chr = over$chr,
  #  start = over[, ifelse(start > i.start, start, i.start)],
   # end = over[, ifelse(end < i.end, end, i.end)])

#tfun<-function(intv, ag) { c(intv,sum(ranges(findOverlaps(IRanges(start=intv[1],end=intv[2]),ag),IRanges(start=intv[1],end=intv[2]),ag)@width)/40000) }

#apply(tevals,1,tfun,ag=agtr)
#Error in apply(tevals, 1, tfun, ag = agtr) : 
 # dim(X) must have a positive length


 #tfun(tevals, ag=agtr) 
#Error in ranges(findOverlaps(IRanges(start = intv[1], end = intv[2]),  : 
#  error in evaluating the argument 'x' in selecting a method for function 'ranges': Error in findOverlaps(IRanges(start = intv[1], end = intv[2]), ag) : 
#  error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': Error in IRanges(start = intv[1], end = intv[2]) : 
#  'end' and 'width' must be NULLs when 'start' is a Ranges object


# i don't think this is the desired output..
# tfun(tevals[[1]], ag=agtr) 
#   [1]    0    1    2    3    4    5    6    7    8    9   10   11   12   13
#  [15]   14   15   16   17   18   19   20   21   22   23   24   25   26   27
#  [29]   28   29   30   31   32   33   34   35   36   37   38   39   40   41
#  [43]   42   43   44   45   46   47   48   49   50   51   52   53   54   55

#tfun(tevals[[1]], ag=agtr[[1]]) 
#Error in ranges(findOverlaps(IRanges(start = intv[1], end = intv[2]),  : 
#  error in evaluating the argument 'x' in selecting a method for function 'ranges': Error in (function (classes, fdef, mtable)  : 
#  unable to find an inherited method for function ‘findOverlaps’ for signature ‘"IRanges", "integer"’

#>  tfun(tevals, ag=agtr[1]) 
#Error in ranges(findOverlaps(IRanges(start = intv[1], end = intv[2]),  : 
#  error in evaluating the argument 'x' in selecting a method for function 'ranges': Error in findOverlaps(IRanges(start = intv[1], end = intv[2]), ag) : 
#  error in evaluating the argument 'query' in selecting a method for function 'findOverlaps': Error in IRanges(start = intv[1], end = intv[2]) : 
#  'end' and 'width' must be NULLs when 'start' is a Ranges object

# is this the way to go?
#intersect" them, finding the set of regions formed by the overlaps 
#  intersect(tevals,agtr)
#IRanges of length 6266
#          start      end width
#[1]     9440000  9440016    17
#[2]     9440021  9440062    42
#[3]     9440127  9440195    69
#[4]     9440201  9440253    53
#[5]     9440258  9440266     9
#...         ...      ...   ...
#[6262] 48043601 48043999   399
#[6263] 48080000 48081856  1857
#[6264] 48082021 48082672   652
#[6265] 48082756 48083558   803
#[6266] 48083741 48083999   259

#------------------------------------------------------------------------------------------------------------------------

# next steps to continue with, once i get the above to work 
    avals<-t(avals)
    aweight1<-as.data.frame(cbind(chr,avals[,c(1)]))
    aweight2<-as.character(as.double(avals[,3]))
    aweight3<-cbind(aweight1,format(aweight2, digits=3, scientific=F))
    load(file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
    alweig[[chr]]<-aweight3
    save(alweig,file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
#    }
#  }
  print(Sys.time())
  q()

 

## all weigs
  rst<-list()
  for (yy in (1:length(atyp))) {
    typ=atyp[yy];    rst[[yy]]<-list()
    for (ft in c(1,2)) {
      spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1];      nchr<-c(1:22,"X")
      if (length(nchr)>0) { rst[[yy]][[ft]]<-paste(yy,ft,nchr,sep="_") }} }
rst<-unlist(rst)
write.table(rst,"/scratch/devel/mkuhlwilm/weight.txt",sep="\t",row.names=F,col.names=F,quote=F)  
  

##  missing ones for repeat
  rst<-list()
  yy=as.numeric(chrom)
  for (yy in (1:length(atyp))) {
    typ=atyp[yy]
    rst[[yy]]<-list()
    for (ft in c(1,2)) {
      spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
      print(c(spec,ft))
      load(file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
      nchr<-newchr[which(newchr%ni%names(alweig))]
      if (length(nchr)>0) { rst[[yy]][[ft]]<-paste(yy,ft,nchr,sep="_")
      }} }
  
  rst<-unlist(rst)
  write.table(rst,"/scratch/devel/mkuhlwilm/weightrest.txt",sep="\t",row.names=F,col.names=F,quote=F)  

  alli<-list()
  for (i in (1:5)) {  alli[[i]]<-paste(1:23,i,sep="_")  }
  write.table(unlist(alli),"/scratch/devel/mkuhlwilm/mutlist.txt",sep="\t",row.names=F,col.names=F,quote=F)  
  

#####################################################
## merge and write into merged species-specific files
library("GenomicRanges")
options("scipen"=100)
options(scipen=1, digits=8)

for (ft in c(2)) {
  print(ft)
  for (yy in (1:length(atyp))) {
      typ=atyp[yy]; print(typ)
      spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
      load(file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
      nweig<-list()
      for (i in c(1:22,"X")) { tw<-alweig[[which(names(alweig)==i)]]; nweig[[i]]<-tw[-c((nrow(tw)-1):(nrow(tw))),];colnames(nweig[[i]])<-c("chr","V2","format")      }
      aweight<-do.call(rbind,nweig)
      aweight1<-as.data.frame(matrix(unlist(aweight[,1:2]),ncol=2))
  #    aweight2<-as.character(as.double(as.numeric(as.character(aweight[,3]))))
      aweight3<-cbind(aweight1,format(round(as.numeric(as.character(aweight[,3])), digits=3),nsmall=3, scientific=F))
      write.table(aweight3,paste("/scratch/devel/mkuhlwilm/arch/N",ft,"_",spec,"_weights_float.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)  
    }
  }

