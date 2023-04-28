
#-----------------------------------------------------------------------------------------------------------------------
# filtering by windows which are sufficiently informative
# following from: /scratch/devel/hpawar/admix/abc/simul/scripts/get_clbl.1jul.R
              #   /scratch/devel/hpawar/admix/abc/simul/scripts/get_clbl.1jul.arr
#-----------------------------------------------------------------------------------------------------------------------
#/scratch/devel/hpawar/admix/abc/emp.data/window.info/ : info per chr, of proportion of sites per window, after filtering by repeats & mappability

#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)
files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/", pattern="gorilla_"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-aweight3
}

#-----------------------------------------------------------------------------------------------------------------------

factortonumeric_fun<-function(x){
genomewide[[x]][,c(1:3)] <- sapply(genomewide[[x]][,c(1:3)] ,as.character)
genomewide[[x]][,c(1:3)] <- sapply(genomewide[[x]][,c(1:3)] ,as.numeric)
return(genomewide[[x]])
}

genomewide_num=list()
for (j in 1:22){
genomewide_num[[j]]<-factortonumeric_fun(j)
}

test<-data.frame(do.call(rbind.data.frame, genomewide_num))
colnames(test)<-c('chr','startpos','proportion')

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

out_test=list()
for (i in 1:nrow(first_windows_df)) {
  
  skip_to_next <- FALSE

  
  tryCatch(

out_test[[i]]<-cbind(match_fun1(i)[[1]], match_fun1(i)[[2]])

, error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
}


match_hold<-data.frame(do.call(rbind.data.frame, out_test))

colnames(match_hold)[4]<-'identifier'


# ie subset lists of x by match_hold$identifier
ns <-match_hold[,4]

stats_inform=list()
for (n in ns) {
stats_inform[[n]]<-x[[n]]
}

stats_inform1<-stats_inform[lengths(stats_inform) != 0]
# these are the ones (sufficiently informative windows + summary statistics) to retain - ie write this out to file **

save(stats_inform1,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/seg_chr",chrom,"_informativewindows",sep=""))
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

all_function(1)

for (x in 2:22){
all_function(x)
}

