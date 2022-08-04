# Wed 29 Sep 2021 08:43:21 CEST
# rerunning after regenerating the summary stats - to rectify issues - see make.abc.model.27sep.sh

#Tue 28 Sep 2021 09:34:29 CEST
# need to run the following interactively once have output (ie rerun filter.het.by.clbl.26aug.R - here renamed as filter.segs.by.clbl.27sep.R)
#- number of population-wise fixed sites and the number of population-wise segregating sites; 
#- fixed sites per individual 
#ls -lhtr /scratch/devel/hpawar/admix/abc/emp.data/test/segs/
#-----------------------------------------------------------------------------------------------------------------------

# #!/usr/bin/r
# Wed  1 Sep 2021 15:04:19 CEST - need to rerun filter.het.by.clbl.26aug.R
# after recalculating seg sites - removing positions fixed between all gor subspecies

# Tue 31 Aug 2021 12:56:03 CEST
# # rerun filter.het.by.clbl.26aug.R for regenerated het & seg sites stats per window

# Thu 26 Aug 2021 16:34:38 CEST
# filter by informative windows for new het estimates

# Wed 18 Aug 2021 09:56:49 CEST
  # refiltering by informative windows after regenerating summary stats for empirical data
  # using  /scratch/devel/hpawar/admix/abc/simul/scripts/stats.emp.data.clean2.R
  # called by /scratch/devel/hpawar/admix/abc/simul/scripts/stats.emp.data.arr
  # b/c order of individuals in the vcf != order of individuals simulated   
      # waiting for this to regenerate - regenerated Wed 18 Aug 2021 11:19:41 CEST
#-----------------------------------------------------------------------------------------------------------------------
#Wed 18 Aug 2021 14:51:51 CEST
# after sending windows.info.issues.18aug, windows.18aug.png, windows.18aug.1.png to MK
#Filtering the windows by >0.75 may be too strict.

# MK
#yes, it seems there is a small amount of informative sites per window. 
#It seems quite centered around 50%, which is strange, but nothing we can do about it, I guess. 
#In this situation, it would make sense to apply a threshold of 50%, to keep a large enough part of the genome.

#-----------------------------------------------------------------------------------------------------------------------
  # Thu  8 Jul 2021 11:23:18 CEST
# filtering by which are sufficiently informative
# following from: /scratch/devel/hpawar/admix/abc/simul/scripts/get_clbl.1jul.R
              #   /scratch/devel/hpawar/admix/abc/simul/scripts/get_clbl.1jul.arr
#-----------------------------------------------------------------------------------------------------------------------
#/scratch/devel/hpawar/admix/abc/emp.data/window.info/ : info per chr, of proportion of sites per window, after filtering by repeats & mappability

#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)
files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/", pattern="gorilla_"), sep = "")

#load(file=files[[1]],verbose=T)
#x<-aweight3


genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-aweight3
}


# str(genomewide)
#List of 22
#$ :'data.frame':  3389 obs. of  3 variables:
#  ..$ chr                                         : Factor w/ 1 level "10": 1 1 1 1 1 1 1 1 1 1 ...
#  ..$ V2                                          : Factor w/ 3389 levels "0","1000000",..: 1 1723 2834 557 1057 1168 1279 1390 1501 1612 ...
#  ..$ format(aweight2, digits = 3, scientific = F): Factor w/ 3044 levels "0       ","0.000075",..: 1 45 101 1000 1788 786 2751 2757 2845 2205 ...
# $ :'data.frame': 3376 obs. of  3 variables:
#  ..$ chr                                         : Factor w/ 1 level "11": 1 1 1 1 1 1 1 1 1 1 ...
#  ..$ V2                                          : Factor w/ 3376 levels "0","1000000",..: 1 1710 2821 557 1044 1155 1266 1377 1488 1599 ...
#  ..$ format(aweight2, digits = 3, scientific = F): Factor w/ 3019 levels "0       ","0.000125",..: 1 6 26 29 137 2464 1132 2293 223 2200 ...
# $ :'data.frame': 3347 obs. of  3 variables:

# head(genomewide[[1]])
#  chr     V2 format(aweight2, digits = 3, scientific = F)
#1  10      0                                     0       
#2  10  40000                                     0.019475
#3  10  80000                                     0.08385 
#4  10 120000                                     0.437775
#5  10 160000                                     0.542675
#6  10 200000                                     0.405 
# col 3 = the proportion of each window that has data. 

# MK : only use windows with >3/4 of positions with information to get the summary stats. Otherwise, you would compare a lot of "half-empty" windows to "complete" simulated windows.

# assess dist of col3

# find windows which col 3 >= 0.75

#-----------------------------------------------------------------------------------------------------------------------

# write this as a function & apply to chr 1:22

# convert all columns from factor to character then to numeric variables
factortonumeric_fun<-function(x){
genomewide[[x]][,c(1:3)] <- sapply(genomewide[[x]][,c(1:3)] ,as.character)
genomewide[[x]][,c(1:3)] <- sapply(genomewide[[x]][,c(1:3)] ,as.numeric)
return(genomewide[[x]])
}

# loop through
genomewide_num=list()
for (j in 1:22){
genomewide_num[[j]]<-factortonumeric_fun(j)
}

#str(genomewide_num)
#List of 22
# $ :'data.frame': 3389 obs. of  3 variables:
#  ..$ chr                                         : num [1:3389] 10 10 10 10 10 10 10 10 10 10 ...
#  ..$ V2                                          : num [1:3389] 0 40000 80000 120000 160000 200000 240000 280000 320000 360000 ...
#  ..$ format(aweight2, digits = 3, scientific = F): num [1:3389] 0 0.0195 0.0838 0.4378 0.5427 ...
# $ :'data.frame': 3376 obs. of  3 variables:
#  ..$ chr                                         : num [1:3376] 11 11 11 11 11 11 11 11 11 11 ...
#  ..$ V2                                          : num [1:3376] 0 40000 80000 120000 160000 200000 240000 280000 320000 360000 ...
#  ..$ format(aweight2, digits = 3, scientific = F): num [1:3376] 0 0.000725 0.007975 0.029925 0.2069 ...


test<-data.frame(do.call(rbind.data.frame, genomewide_num))
colnames(test)<-c('chr','startpos','proportion')

#head(test)
#  chr startpos proportion
#1  10        0   0.000000
#2  10    40000   0.019475

#tail(test)
#      chr  startpos proportion
#72031   9 141000000   0.443375
#72032   9 141040000   0.128375

#informative<-test[ which(test$proportion > 0.75), ]

#head(informative)
#   chr startpos proportion
#11  10   400000   0.817450
#14  10   520000   0.791800

# nrow(test)
#[1] 72036
# nrow(informative)
#[1] 2006

# nrow(informative)/nrow(test)
#[1] 0.02784719
# quite a loss of data

#------------------------------------------------------
# MK: apply a threshold of 50%, to keep a large enough part of the genome
informative<-test[ which(test$proportion > 0.5), ]

 nrow(test)
#[1] 72036
 nrow(informative)
#[1] 32561

 nrow(informative)/nrow(test)
 #[1] 0.4520101
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

# read in windows of the summary stats for empirical data
#/scratch/devel/hpawar/admix/abc/emp.data/stats_chr

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# run the function:
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


statsfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", pattern="seg_chr"), sep = "")

all_function<-function(nput) {
load(file=statsfiles[[nput]],verbose=T)
x<-stats.out
# first just read in the window chr, number, start pos, end pos - for which I have calculated summary statistics 
first_windows=list()
for (j in 1:length(stats.out)){
first_windows[[j]]<-head(x[[j]][[1]])
}

# convert all columns to numeric
first_windows_df<-data.frame(do.call(rbind.data.frame, first_windows))
first_windows_df<- sapply(first_windows_df[,c(1:4)] ,as.numeric)

# find which windows are informative for this chromosome
chrom=as.numeric(x[[1]][[1]][[1]])
chr1.inform<-informative[ which(informative$chr == chrom), ]

# compare the 2 dfs, ie find which summary stats windows overlap with informative window (after filtering for repeats & mappability)
match_fun1<-function(x){
out<-chr1.inform[ which(chr1.inform$startpos >= first_windows_df[x,3] & chr1.inform$startpos <= first_windows_df[x,4]), ]
return(list(out,x)) # add identifier, which row of first_windows_df has an informative window, to reconstruct below
}

# ignore error when looping
out_test=list()
for (i in 1:nrow(first_windows_df)) {
  
  skip_to_next <- FALSE
  
  # Note that print(b) fails since b doesn't exist
  
  tryCatch(

out_test[[i]]<-cbind(match_fun1(i)[[1]], match_fun1(i)[[2]])

, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}


match_hold<-data.frame(do.call(rbind.data.frame, out_test))

colnames(match_hold)[4]<-'identifier'

#head(match_hold)
#      chr startpos proportion identifier
#25982   1   840000   0.841050          4
#25988   1  1080000   0.796950          7


# ie subset lists of x by match_hold$identifier
ns <-match_hold[,4]

# loop over numbers in vector - https://stackoverflow.com/questions/14592383/using-non-sequential-vector-as-input-for-a-loop
stats_inform=list()
for (n in ns) {
stats_inform[[n]]<-x[[n]]
}

 # => remove null elements of the list
#mylist[lengths(mylist) != 0] # https://stackoverflow.com/questions/33004238/r-removing-null-elements-from-a-list

stats_inform1<-stats_inform[lengths(stats_inform) != 0]
# these are the ones (sufficiently informative windows + summary statistics) to retain - ie write this out to file **
  # & need to generalise so can apply a large function with these steps to the rest of the chromosomes


#mkdir -p /scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/
#chrom=as.numeric(stats_inform1[[1]][[1]][[1]])
#save(stats_inform1,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr",chrom,"_informativewindows",sep=""))
save(stats_inform1,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/seg_chr",chrom,"_informativewindows",sep=""))
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

all_function(1)
 # has generated output check this is correct
 # load(file="/scratch/devel/hpawar/admix/abc/emp.data/test/segs/seg_chr1_informativewindows",verbose=T)
#head( stats_inform1) # against:
# seems fine

for (x in 2:22){
all_function(x)
}

# have all generated
# ls -lhtr /scratch/devel/hpawar/admix/abc/emp.data/test/segs/seg_chr*_informativewindows
