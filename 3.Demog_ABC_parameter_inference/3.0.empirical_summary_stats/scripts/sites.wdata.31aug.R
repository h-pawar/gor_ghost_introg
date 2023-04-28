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
# MK: apply a threshold of 50%, to keep a large enough part of the genome
informative<-test[ which(test$proportion > 0.5), ]


#-----------------------------------------------------------------------------------------------------------------------


statsfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/", pattern="het_chr"), sep = "")

proportion_function<-function(nput) {
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


match_hold$proportion*40000

return(sum(match_hold$proportion*40000))

ny <-c(1,2,4, 6, 8, 10,12,14,16, 18, 20,23,24,26, 28, 31,33,35,37,39,41,43)

hold_bp=list()
for (x in ny){
hold_bp[[x]]<-proportion_function(x)
}


hold_bp1<-hold_bp[lengths(hold_bp) != 0]


sum(unlist(hold_bp1))
#[1] 769239611
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
#MK: If they are NaN (not a number), they are due to division by zero (I guess actually often 0/0) and might be corrected to 0 instead. 
find.nas_fun1<-function(x) {
for (j in 1:length(genomewide[[x]])){
genomewide[[x]][[j]][[1]]<-as.numeric(genomewide[[x]][[j]][[1]])
}
out_1<-list()
for (j in 1:length(genomewide[[x]])){
out_1[[j]]<-lapply(genomewide[[x]][[j]], function(y) replace(y, !is.finite(y), 0))
}
return(out_1)
}


finalwindows=list()
for (i in 1:length(genomewide)){
finalwindows[[i]]<-find.nas_fun1(i)
}

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


fst_df<-stats_asdf_fun(4)

inf_fsts=list()
for (i in 1:ncol(fst_df)){
inf_fsts[[i]]<-which(fst_df[,i] > 1)
}

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


hold_wind=list()
for (x in ny){
hold_wind[[x]]<- wind_df[x,]
}

hold_wind1<-hold_wind[lengths(hold_wind) != 0]

hold_wind1_df<-data.frame(matrix(unlist(hold_wind1), nrow=length(hold_wind1), byrow=TRUE),stringsAsFactors=FALSE)
# these are the windows with infinite fsts - extract proportion of sites in these windows 


#-----------------------------------------------------------------------------------------------------------------------

output_matchhold_function<-function(nput) {
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

return(match_hold)
}


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
# read in the informative windows in order of chr 

hold_matchhold=list()
for (x in ny){
hold_matchhold[[x]]<-output_matchhold_function(x)
}

hold_matchhold1<-hold_matchhold[lengths(hold_matchhold) != 0]

#-----------------------------------------------------------------------------------------------------------------------


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

unlist(chr17_prop_inffst),
unlist(chr18_prop_inffst),

unlist(chr20_prop_inffst)

)
#[1] 1509532

#-----------------------------------------------------------------------------------------------------------------------

sum(unlist(hold_bp1))
#[1] 769239611

# total number of informative sites - (sum of sites of windows with infinite fsts that have been filtered out) 
#769239611 - 1509532
#1] 767730079


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


# v similar vals w & w/out removing windows w infinite fsts (b/c only small proportion of windows have infinite fsts)

#-----------------------------------------------------------------------------------------------------------------------
# note this is considering proportion of sites per window (rather than per individual)
