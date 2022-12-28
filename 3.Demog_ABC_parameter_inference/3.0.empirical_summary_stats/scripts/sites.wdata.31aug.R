# following filter.het.by.clbl.26aug.R -> to determine proportion of sites in the informative windows

# 1) calculate sum of all sites with data in each individual -> the below is rather all sites with data for the informative windows (rather than per individual) 
# 2) identify which windows had infinite fsts -> need to deduct their proportion*40000 from the total
# 3) sum of all heterozygous sites divided by the sum of all sites with data
#-----------------------------------------------------------------------------------------------------------------------

#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------

files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/", pattern="gorilla_"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-aweight3
}

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

test<-data.frame(do.call(rbind.data.frame, genomewide_num))
colnames(test)<-c('chr','startpos','proportion')
# MK: apply a threshold of 50%, to keep a large enough part of the genome
informative<-test[ which(test$proportion > 0.5), ]


#-----------------------------------------------------------------------------------------------------------------------


statsfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/", pattern="het_chr"), sep = "")

#statsfiles[[28]]
#-----------------------------------------------------------------------------------------------------------------------
#[1] "/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr22"
#nput=28

proportion_function<-function(nput) {
#all_function<-function(nput) {
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

#> length(stats_inform1)
#[1] 386
#> length(ns)
#[1] 386
#> length(match_hold)
#[1] 4
#> nrow(match_hold)
#[1] 386
# ie go from the match_hold dfs to obtain # snps per window
#head(match_hold$proportion)
#[1] 0.570375 0.619275 0.505325 0.520325 0.582400 0.700125

match_hold$proportion*40000

return(sum(match_hold$proportion*40000))
#[1] 9347158

# run this per chr 
  # but this will include windows with infinite fsts 

#mkdir -p /scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/
#chrom=as.numeric(stats_inform1[[1]][[1]][[1]])
#save(stats_inform1,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr",chrom,"_informativewindows",sep=""))
}

ny <-c(1,2,4, 6, 8, 10,12,14,16, 18, 20,23,24,26, 28, 31,33,35,37,39,41,43)

hold_bp=list()
for (x in ny){
hold_bp[[x]]<-proportion_function(x)
}


hold_bp1<-hold_bp[lengths(hold_bp) != 0]

#unlist(hold_bp1)
# [1] 62554110 42484251 37366865 34833686 32906709 25718818 23430970 20901450
# [9] 22645167 26784431  8203157 77279374 16854681 11497984  9347158 58631793
#[17] 51936429 53501132 53442018 44807229 41973768 12138431

sum(unlist(hold_bp1))
[1] 769239611
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# estimasting het from the total (pre minus infinite fst windows)

 (colSums(hetperid_wl)/769239611)*1000
       X1        X2        X3        X4        X5        X6        X7        X8 
0.9714268 0.9468961 0.9627182 0.9432419 0.9931704 0.9671356 0.9387439 0.9938464 
       X9       X10       X11       X12       X13       X14       X15       X16 
0.9338872 0.8781009 0.9305891 0.9178389 0.8995988 0.9169354 0.8694079 0.9272572 
      X17       X18       X19       X20       X21       X22       X23       X24 
0.9423761 0.7262315 0.9182600 0.8857760 0.7536013 0.8308919 0.9277369 0.9327770 
      X25       X26       X27 
0.8742101 0.9292566 0.8969013 


mean((colSums(hetperid_wl)/769239611))*1000
[1] 0.9114375
> sd((colSums(hetperid_wl)/769239611))*1000
[1] 0.06202543



cbind(mean((colSums(hetperid_el)/769239611))*1000,sd((colSums(hetperid_el)/769239611))*1000)
cbind(mean((colSums(hetperid_em)/769239611))*1000,sd((colSums(hetperid_em)/769239611))*1000)

(sum(hetperid_wc)/769239611)*1000

> cbind(mean((colSums(hetperid_el)/769239611))*1000,sd((colSums(hetperid_el)/769239611))*1000)
          [,1]       [,2]
[1,] 0.4120228 0.03061125
> cbind(mean((colSums(hetperid_em)/769239611))*1000,sd((colSums(hetperid_em)/769239611))*1000)
          [,1]       [,2]
[1,] 0.4298602 0.03338122
> 
> (sum(hetperid_wc)/769239611)*1000
[1] 0.700753



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2) identify which windows had infinite fsts -> need to deduct their proportion*40000 from the total

files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/subset/", pattern="_informativewindows"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-stats_inform1
}

#replace NaN values with 0
find.nas_fun1<-function(x) {
for (j in 1:length(genomewide[[x]])){
    # convert window info from character -> numeric - so do not lose this info when replacing non-finite numbers with 0
genomewide[[x]][[j]][[1]]<-as.numeric(genomewide[[x]][[j]][[1]])
}
out_1<-list()
for (j in 1:length(genomewide[[x]])){
out_1[[j]]<-lapply(genomewide[[x]][[j]], function(y) replace(y, !is.finite(y), 0))
}
return(out_1)
}

#MK: If they are NaN (not a number), they are due to division by zero (I guess actually often 0/0) and might be corrected to 0 instead. 

# run for all 22 chr
finalwindows=list()
for (i in 1:length(genomewide)){
finalwindows[[i]]<-find.nas_fun1(i)
}

#> str(genomewide_out)
#List of 22
# $ :List of 1735
#  ..$ :List of 6
#-----------------------------------------------------------------------------------------------------------------------

format_stats_fun<-function(x,y) {
otest=list()
for (i in 1:length(finalwindows[[x]])){
otest[[i]]<-t(as.data.frame(finalwindows[[x]][[i]][[y]]))
}
outest<-data.frame(matrix(unlist(otest), nrow=length(otest), byrow=TRUE),stringsAsFactors=FALSE)
return(outest)
}

stats_asdf_fun<-function(x){
stats_o=list()
for (i in 1:22){
stats_o[[i]]<-format_stats_fun(i,x)
}
ptest<-do.call(rbind, stats_o)
return(ptest)
}

# summary stats for all retained windows (genome-wide)
#het_df<-stats_asdf_fun(2)
#segsites_df<-stats_asdf_fun(3)
fst_df<-stats_asdf_fun(4)
#pi_df<-stats_asdf_fun(5)
#tajima_df<-stats_asdf_fun(6)

inf_fsts=list()
for (i in 1:ncol(fst_df)){
inf_fsts[[i]]<-which(fst_df[,i] > 1)
}

#sort(unlist(inf_fsts))
# [1]   551   741   984  1820  2339  2781  3202  3216  4047  4047  4461  5887
#[13]  6161  7976  9234  9234  9888 10136 10235 10304 10743 12926 12973 13803
#[25] 14435 15478 16690 16899 17010 17010 17369 17369 17402 18040 18275 19526
#[37] 19781 19860 20050 20353 20745 20878 20977 22584 22734 22779 22842 23572
#[49] 24226 24383 24901 24901 25196 25276 27006 27006 27104 28020 28474 28488
#[61] 28661 29120 29614 29624 29645 30666 30822 31228

#-----------------------------------------------------------------------------------------------------------------------

fst_df1<-fst_df[-c(sort(unlist(inf_fsts))),]


# extract window identifiers for these infinite fst windows

format_stats_fun<-function(x,y) {
otest=list()
for (i in 1:length(finalwindows[[x]])){
otest[[i]]<-t(as.data.frame(finalwindows[[x]][[i]][[y]]))
}
outest<-data.frame(matrix(unlist(otest), nrow=length(otest), byrow=TRUE),stringsAsFactors=FALSE)
return(outest)
}

stats_asdf_fun<-function(x){
stats_o=list()
for (i in 1:22){
stats_o[[i]]<-format_stats_fun(i,x)
}
ptest<-do.call(rbind, stats_o)
return(ptest)
}

wind_df<-stats_asdf_fun(1)

ny<-c(sort(unlist(inf_fsts)))

> wind_df[551,]
    X1   X2       X3       X4
551 10 1132 45338713 45378713


hold_wind=list()
for (x in ny){
hold_wind[[x]]<- wind_df[x,]
}

hold_wind1<-hold_wind[lengths(hold_wind) != 0]

hold_wind1_df<-data.frame(matrix(unlist(hold_wind1), nrow=length(hold_wind1), byrow=TRUE),stringsAsFactors=FALSE)
#Ω these are the windows with infinite fsts - extract proportion of sites in these windows 


#-----------------------------------------------------------------------------------------------------------------------




 #statsfiles[2]
#[1] "/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr10"


output_matchhold_function<-function(nput) {
#all_function<-function(nput) {
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

return(match_hold)
}


#ny <-c(1,2,4, 6, 8, 10,12,14,16, 18, 20,23,24,26, 28, 31,33,35,37,39,41,43)


ny <-c(1, 
23, 
31,                    
33,                   
35,                    
37,                    
39,                    
41,                    
43,                   
2,                   
4,                   
6,                   
8,                   
10,                   
12,                   
14,                   
16,                   
18,                   
20,                                     
24,                   
26,                   
28  )
# read in the informative windows in order of chr - still not in order of chr

hold_matchhold=list()
for (x in ny){
hold_matchhold[[x]]<-output_matchhold_function(x)
}

hold_matchhold1<-hold_matchhold[lengths(hold_matchhold) != 0]

#-----------------------------------------------------------------------------------------------------------------------

#hold_wind1_df # windows with infinite fsts

#hold_matchhold1[[2]][which(hold_matchhold1[[2]]$identifier==hold_wind1_df[1,2]),] 

#hold_matchhold1[[2]][which(hold_matchhold1[[2]]$identifier==hold_wind1_df[1,2]),3]*40000
#[1] 25554

#prop_inf_fsts_perchr_fun<-function(x,y){
#z<-hold_matchhold1[[x]][which(hold_matchhold1[[x]]$identifier==hold_wind1_df[y,2]),3]*40000
#return(z)
#}

#prop_inf_fsts_perchr_fun(2,1)

# for chr 10

#identifiers<-hold_wind1_df[,2]
#chrs<- hold_wind1_df[,1]

#chr10_prop_inffst=list()
#for (y in 1:3){
#chr10_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(2,y)
#}


#chr11_prop_inffst=list()
#for (y in 4:8){
#chr11_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(3,y)
#}

#hold_matchhold1[[3]][which(hold_matchhold1[[3]]$identifier==hold_wind1_df[4,2]),3]*40000
  # q - why these identifier numbers don't correspond??

# hold_wind1_df[4,]
#  X1  X2      X3      X4
#4 11 155 6361245 6401245


#hold_matchhold1[[3]][ which(hold_matchhold1[[3]]$startpos >= hold_wind1_df[4,3] & hold_matchhold1[[3]]$startpos <= hold_wind1_df[4,4]), ]
#chr startpos proportion identifier
#3550  11  6400000    0.52345        146

# maybe go through like this instead?

#hold_matchhold1[[2]][ which(hold_matchhold1[[2]]$startpos >= hold_wind1_df[1,3] & hold_matchhold1[[2]]$startpos <= hold_wind1_df[1,4]), ]

#-----------------------------------------------------------------------------------------------------------------------

prop_inf_fsts_perchr_fun<-function(x,y){
y1<-hold_matchhold1[[x]][ which(hold_matchhold1[[x]]$startpos >= hold_wind1_df[y,3] & hold_matchhold1[[x]]$startpos <= hold_wind1_df[y,4]), ]
y2<-y1[3]*40000
return(y2)
}

chr1_prop_inffst=list()
for (y in 20:22){
chr1_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(1,y)
}

chr10_prop_inffst=list()
for (y in 1:3){
chr10_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(2,y)
}

chr11_prop_inffst=list()
for (y in 4:8){
chr11_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(3,y)
}

chr12_prop_inffst=list()
for (y in 9:10){
chr12_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(4,y)
}

# 5-13
for (i in 5:22){
print(head(hold_matchhold1[[i]][1,1]))
}
#[1] 13
#[1] 14
#[1] 15
#[1] 16
#[1] 17
#[1] 18
#[1] 19
#[1] 2
#[1] 20
#[1] 21
#[1] 22
#[1] 3
#[1] 4
#[1] 5
#[1] 6
#[1] 7
#[1] 8
#[1] 9

chr13_prop_inffst<-prop_inf_fsts_perchr_fun(5,11)
chr14_prop_inffst<-prop_inf_fsts_perchr_fun(6,12)
chr15_prop_inffst<-prop_inf_fsts_perchr_fun(7,13)
chr17_prop_inffst<-prop_inf_fsts_perchr_fun(9,14)

chr18_prop_inffst=list()
for (y in 15:19){
chr18_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(10,y)
}

chr2_prop_inffst=list()
for (y in 24:31){
chr2_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(12,y)
}

chr20_prop_inffst<-prop_inf_fsts_perchr_fun(13,23)

chr3_prop_inffst=list()
for (y in 32:39){
chr3_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(16,y)
}


chr4_prop_inffst=list()
for (y in 40:43){
chr4_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(17,y)
}

chr5_prop_inffst=list()
for (y in 44:49){
chr5_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(18,y)
}

chr6_prop_inffst=list()
for (y in 50:51){
chr6_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(19,y)
}

chr7_prop_inffst=list()
for (y in 52:56){
chr7_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(20,y)
}

chr8_prop_inffst=list()
for (y in 57:61){
chr8_prop_inffst[[y]]<-prop_inf_fsts_perchr_fun(21,y)
}

chr9_prop_inffst<-prop_inf_fsts_perchr_fun(22,62)

sum(
unlist(chr1_prop_inffst),
unlist(chr2_prop_inffst),
unlist(chr3_prop_inffst),
unlist(chr4_prop_inffst),
unlist(chr5_prop_inffst),
unlist(chr6_prop_inffst),
unlist(chr7_prop_inffst),
unlist(chr8_prop_inffst),
unlist(chr9_prop_inffst),
unlist(chr10_prop_inffst),
unlist(chr11_prop_inffst),
unlist(chr12_prop_inffst),
unlist(chr13_prop_inffst),
unlist(chr14_prop_inffst),
unlist(chr15_prop_inffst),
#unlist(chr16_prop_inffst),
unlist(chr17_prop_inffst),
unlist(chr18_prop_inffst),
#unlist(chr19_prop_inffst),
unlist(chr20_prop_inffst)
#unlist(chr21_prop_inffst),
#unlist(chr22_prop_inffst)
)
[1] 1509532

#-----------------------------------------------------------------------------------------------------------------------

sum(unlist(hold_bp1))
[1] 769239611

# total number of informative sites - (sum of sites of windows with infinite fsts that have been filtered out) 
769239611 - 1509532
1] 767730079


#-----------------------------------------------------------------------------------------------------------------------

# 2) read in heterozygosity files 
    # generated using
        # /scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.31aug.R, 
        #called by - /scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr
        #  #/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.het.by.26aug.R - ran interactively
# only the informative windows

hetfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/", pattern="_informativewindows"), sep = "")

hetgenomewide=list()
for (i in 1:length(hetfiles)){
load(file=hetfiles[[i]],verbose=T)
hetgenomewide[[i]]<-stats_inform1
}
# 2.1) calc heterozygosity per individual

format_hetperid_fun<-function(x,y) {
htest=list()
for (i in 1:length(hetgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(hetgenomewide[[x]][[i]][[2]][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


hetperid_asdf_fun<-function(y){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_hetperid_fun(i,y)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}

hetperid_wl<-hetperid_asdf_fun(1)
hetperid_wc<-hetperid_asdf_fun(2)
hetperid_el<-hetperid_asdf_fun(3)
hetperid_em<-hetperid_asdf_fun(4)

hetperid_wl1<-hetperid_wl[-c(sort(unlist(inf_fsts))),]
hetperid_wc1<-hetperid_wc[-c(sort(unlist(inf_fsts))),]
hetperid_el1<-hetperid_el[-c(sort(unlist(inf_fsts))),]
hetperid_em1<-hetperid_em[-c(sort(unlist(inf_fsts))),]

# 3) sum of all heterozygous sites divided by the sum of all sites with data
# *1000 for het per kb

cbind(mean((colSums(hetperid_wl1)/767730079))*1000,sd((colSums(hetperid_wl1)/767730079))*1000)
cbind(mean((colSums(hetperid_el1)/767730079))*1000,sd((colSums(hetperid_el1)/767730079))*1000)
cbind(mean((colSums(hetperid_em1)/767730079))*1000,sd((colSums(hetperid_em1)/767730079))*1000)

sum(hetperid_wc1/767730079)*1000


> cbind(mean((colSums(hetperid_wl1)/767730079))*1000,sd((colSums(hetperid_wl1)/767730079))*1000)
          [,1]       [,2]
[1,] 0.9108459 0.06199505
> cbind(mean((colSums(hetperid_el1)/767730079))*1000,sd((colSums(hetperid_el1)/767730079))*1000)
          [,1]       [,2]
[1,] 0.4115742 0.03066836
> cbind(mean((colSums(hetperid_em1)/767730079))*1000,sd((colSums(hetperid_em1)/767730079))*1000)
          [,1]       [,2]
[1,] 0.4295022 0.03339686
> 

> sum(hetperid_wc1/767730079)*1000
[1] 0.7008401

# v similar vals w & w/out removing windows w infinite fsts (b/c only small proportion of windows have infinite fsts)

#-----------------------------------------------------------------------------------------------------------------------
# note this is considering proportion of sites per window (rather than per individual)
