# Wed 29 Sep 2021 09:01:56 CEST
# have regenerated the stats - need to reprocess

# Tue 28 Sep 2021 10:27:35 CEST
# using sites with data per id - from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/sites.wdata.31aug.R
# process new stats (need to normalise fixed sites per id by sites w data)

#-----------------------------------------------------------------------------------------------------------------------


#filter.het.by.clbl.26aug.R -> determined proportion of sites in the informative windows as
# total number of informative sites - (sum of sites of windows with infinite fsts that have been filtered out) 
#769239611 - 1509532
#1] 767730079

# here process new stats 
#- number of population-wise fixed sites and the number of population-wise segregating sites; : may make sense to normalise this as well? 
# - the fixed sites per individual: normalising this by number of sites with data 


# /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/sites.wdata.31aug.R - process stats (determining number of sites with data)
# in sites.wdata.31aug.R these 3 steps performed
# 1) calculate sum of all sites with data in each individual -> the below is rather all sites with data for the informative windows (rather than per individual) 
# 2) identify which windows had infinite fsts -> need to deduct their proportion*40000 from the total
# 3) sum of all heterozygous sites divided by the sum of all sites with data
#-----------------------------------------------------------------------------------------------------------------------

#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# run this section **

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

#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#sum(unlist(hold_bp1))
#[1] 769239611

# total number of informative sites - (sum of sites of windows with infinite fsts that have been filtered out) 
#769239611 - 1509532
#1] 767730079


#-----------------------------------------------------------------------------------------------------------------------
# Tue 28 Sep 2021 11:46:21 CEST
# 2) read in files w newly calculated stats
#- number of population-wise fixed sites and the number of population-wise segregating sites; 
# - the fixed sites per individual
# generated using
    #/scratch/devel/hpawar/admix/abc/simul/scripts/segs.emp.27sep.R  # called by 
    #/scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr
#then ran interactively - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.segs.by.clbl.27sep.R
    # note this is equivalent to /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.het.by.clbl.26aug.R
    # but with paths amended for new statistics

segsfiles <- paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/test/segs/", pattern="_informativewindows"), sep = "")

segsgenomewide=list()
for (i in 1:length(segsfiles)){
load(file=segsfiles[[i]],verbose=T)
segsgenomewide[[i]]<-stats_inform1
}

#-----------------------------------------------------------------------------------------------------------------------
# format: 
# head(segsgenomewide[[1]])
#[[1]]
#[[1]][[1]]
#     chr  i                    
#[1,] "10" "2" "138713" "178713"

#[[1]][[2]]
#     [,1] [,2] [,3] [,4]
#[1,]  160  178  169  166 # fixed sites per popn
#[2,]  173   23  108   98 # seg sites per popn

#[[1]][[3]] # equivalent to what i was outputting for het - but instead here counts of fixed sites per id
#[[1]][[3]][[1]]
# [1] 176 169 169 170 176 174 169 175 171 172 169 174 179 173 179 170 184 177 169
#[20] 180 185 172 172 179 176 176 165

#[[1]][[3]][[2]]
#[1] 178

#[[1]][[3]][[3]]
#[1] 171 175 174 178 174 175 169 178 178

#[[1]][[3]][[4]]
# [1] 171 178 176 178 173 180 174 178 176 180 178 176
#-----------------------------------------------------------------------------------------------------------------------



# 2.1) number of population-wise fixed sites and the number of population-wise segregating sites

#> segsgenomewide[[1]][[1]][[2]]
#     [,1] [,2] [,3] [,4]
#[1,]  159  177  168  166
#[2,]  173   23  108   98
#> segsgenomewide[[1]][[1]][[2]][1,]
#[1] 159 177 168 166

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# A) number of population-wise fixed sites

format_fixedsitesperid_fun<-function(x,y,z) {
htest=list()
for (i in 1:length(segsgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(segsgenomewide[[x]][[i]][[2]][z,][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


fixedsitesperid_asdf_fun<-function(y,z){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_fixedsitesperid_fun(i,y,z)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}



#-----------------------------------------------------------------------------------------------------------------------
# from new simns - /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R 
# full_segsites=rbind(colSums(full_segs[seq(1,4999,2),]),colSums(full_segs[seq(2,5000,2),])) 
#-----------------------------------------------------------------------------------------------------------------------

# output from emp data:
#eseg<-segfun(y)[,c(4,3,2,1)] # in order of WL,WC,EL,EM


fixedsitesperid_wl<-fixedsitesperid_asdf_fun(1,1)
# remove stats from windows with infinite fsts
fixedsitesperid_wl1<-fixedsitesperid_wl[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_wc<-fixedsitesperid_asdf_fun(2,1)
fixedsitesperid_wc1<-fixedsitesperid_wc[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_el<-fixedsitesperid_asdf_fun(3,1)
fixedsitesperid_el1<-fixedsitesperid_el[-c(sort(unlist(inf_fsts))),]
fixedsitesperid_em<-fixedsitesperid_asdf_fun(4,1)
fixedsitesperid_em1<-fixedsitesperid_em[-c(sort(unlist(inf_fsts))),]


# perhaps should normalise by data coverage here also?

rbind(sum(fixedsitesperid_wl1/767730079)*1000,
sum(fixedsitesperid_wc1/767730079)*1000,
sum(fixedsitesperid_el1/767730079)*1000,
sum(fixedsitesperid_em1/767730079)*1000)


# A) number of population-wise fixed sites - before considering 2nd & 3rd alt alleles - after correcting is below
# WL, WC, EL, EM
#        [,1]
#[1,] 6.277046
#[2,] 6.952979
#[3,] 6.759320
#[4,] 6.706769


# A) number of population-wise fixed sites
# WL, WC, EL, EM

#       [,1]
#[1,] 6.318390
#[2,] 6.968396
#[3,] 6.774888
#[4,] 6.723011
#-----------------------------------------------------------------------------------------------------------------------

#Wed  6 Oct 2021 09:30:19 CEST
# after recalculating to  exclude differences to human (1/1 across all individuals). 
# fixedsitesperid will go down by an order of magnitude. For the segsites, nothing changes.
# as MK suggested
    # in make.abc.model.27sep.sh
    # rerunning - #/scratch/devel/hpawar/admix/abc/simul/scripts/segs.emp.27sep.R  # called by 
         #/scratch/devel/hpawar/admix/abc/simul/scripts/het.emp.26aug.arr
    # &/Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/filter.segs.by.clbl.27sep.R

# A) number of population-wise fixed sites
# WL, WC, EL, EM

           [,1]
[1,] 0.07708829
[2,] 0.72709409
[3,] 0.53358597
[4,] 0.48170967

#-----------------------------------------------------------------------------------------------------------------------
# B) number of population-wise segregating sites
# z=2 for B

segsitesperid_wl<-fixedsitesperid_asdf_fun(1,2)
segsitesperid_wl1<-segsitesperid_wl[-c(sort(unlist(inf_fsts))),]
segsitesperid_wc<-fixedsitesperid_asdf_fun(2,2)
segsitesperid_wc1<-segsitesperid_wc[-c(sort(unlist(inf_fsts))),]
segsitesperid_el<-fixedsitesperid_asdf_fun(3,2)
segsitesperid_el1<-segsitesperid_el[-c(sort(unlist(inf_fsts))),]
segsitesperid_em<-fixedsitesperid_asdf_fun(4,2)
segsitesperid_em1<-segsitesperid_em[-c(sort(unlist(inf_fsts))),]

rbind(sum(segsitesperid_wl1/767730079)*1000,
sum(segsitesperid_wc1/767730079)*1000,
sum(segsitesperid_el1/767730079)*1000,
sum(segsitesperid_em1/767730079)*1000)

# B) number of population-wise segregating sites
# WL, WC, EL, EM
#          [,1]
#[1,] 3.9377629
#[2,] 0.6925051
#[3,] 1.2538534
#[4,] 1.4693484

# but issue - in segs.emp.27sep.R
# looking at 0s,1s but not 2nd,3rd alt alleles **  - rectified now - Wed 29 Sep 2021 09:08:50 CEST

# B) number of population-wise segregating sites
# WL, WC, EL, EM
    [,1]
[1,] 3.9409085
[2,] 0.6936383
[3,] 1.2544721
[4,] 1.4698317

#Wed  6 Oct 2021 09:30:19 CEST
# after recalculating to  exclude differences to human (1/1 across all individuals) - segsites stay the same - as expected. 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2.2) fixed sites per individual

#-----------------------------------------------------------------------------------------------------------------------
# from new simns - /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R 
  ## individual fixed sites calculated like this
 # count_fun_f<-function(y,ids) {
 #   count_fun2_f<-function(i,y,ids) { length(which(ids[[y]][[i]]=="11")) }
 #   unlist(lapply(1:length(vec),count_fun2_f,ids=ids,y=y))
 # } 
 
 # full_f=lapply(1:msout[grep("nreps",names(msout))][[1]],count_fun_f,ids=ids)
 # full_fix=rowSums(do.call(cbind,full_f))/(msout[grep("nreps",names(msout))][[1]]*len)
 # mean_fun_f<-function(typ,ffix) { c(mean(ffix[(typ[1]+1):typ[2]]),sd(ffix[(typ[1]+1):typ[2]])) }
 # fix_op<-apply(cofig,1,mean_fun_f,ffix=full_fix)*1000
#-----------------------------------------------------------------------------------------------------------------------

# 2nd issue - not sure i'm outputting the right thing to be equivalen to fix_op **

# eg format is like this - whether should instead be multiple ids per popn? yes - have now rectified
# segsgenomewide[[1]][[1]][[3]]
#[1] 4 3 1 2

#-----------------------------------------------------------------------------------------------------------------------
#now fixed sites are output in the same way as the het per id
#> segsgenomewide[[1]][[1]][[3]]
#[[1]]
# [1] 176 169 169 170 176 174 169 175 171 172 169 174 179 173 179 170 184 177 169
#[20] 180 185 172 172 179 176 176 165

#[[2]]
#[1] 178

#[[3]]
#[1] 171 175 174 178 174 175 169 178 178

#[[4]]
# [1] 171 178 176 178 173 180 174 178 176 180 178 176
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# 2.2) calc fixed sites per individual

# this calculation is equivalent to processing het sites per id in het.emp.31aug.R
format_fixedperid_fun<-function(x,y) {
htest=list()
for (i in 1:length(segsgenomewide[[x]])){
htest[[i]]<-t(as.data.frame(segsgenomewide[[x]][[i]][[3]][[y]]))
}
hutest<-data.frame(matrix(unlist(htest), nrow=length(htest), byrow=TRUE),stringsAsFactors=FALSE)
return(hutest)
}


fixedperid_asdf_fun<-function(y){
stats_h=list()
for (i in 1:22){
stats_h[[i]]<-format_fixedperid_fun(i,y)
}
ptest<-do.call(rbind, stats_h)
return(ptest)
}

fixedperid_wl<-fixedperid_asdf_fun(1)
fixedperid_wc<-fixedperid_asdf_fun(2)
fixedperid_el<-fixedperid_asdf_fun(3)
fixedperid_em<-fixedperid_asdf_fun(4)

fixedperid_wl1<-fixedperid_wl[-c(sort(unlist(inf_fsts))),]
fixedperid_wc1<-fixedperid_wc[-c(sort(unlist(inf_fsts))),]
fixedperid_el1<-fixedperid_el[-c(sort(unlist(inf_fsts))),]
fixedperid_em1<-fixedperid_em[-c(sort(unlist(inf_fsts))),]

# 3) sum of all fixed sites divided by the sum of all sites with data
# *1000 for fixed sites per kb

cbind(mean((colSums(fixedperid_wl1)/767730079))*1000,sd((colSums(fixedperid_wl1)/767730079))*1000)
cbind(mean((colSums(fixedperid_el1)/767730079))*1000,sd((colSums(fixedperid_el1)/767730079))*1000)
cbind(mean((colSums(fixedperid_em1)/767730079))*1000,sd((colSums(fixedperid_em1)/767730079))*1000)
sum(fixedperid_wc1/767730079)*1000


#cbind(mean((colSums(fixedperid_wl1)/767730079))*1000,sd((colSums(fixedperid_wl1)/767730079))*1000)
#         [,1]       [,2]
#[1,] 6.868883 0.03072747
#> cbind(mean((colSums(fixedperid_el1)/767730079))*1000,sd((colSums(fixedperid_el1)/767730079))*1000)
#         [,1]      [,2]
#[1,] 7.121846 0.0156905
#> cbind(mean((colSums(fixedperid_em1)/767730079))*1000,sd((colSums(fixedperid_em1)/767730079))*1000)
#         [,1]       [,2]
#[1,] 7.116987 0.01364289
#> 
#> sum(fixedperid_wc1/767730079)*1000
#[1] 6.968396

# Wed  6 Oct 2021 09:35:24 CEST after recalculating:
# fixedsitesperpop also go down by order of magnitude
 cbind(mean((colSums(fixedperid_wl1)/767730079))*1000,sd((colSums(fixedperid_wl1)/767730079))*1000)
          [,1]       [,2]
[1,] 0.6275814 0.03072747
> cbind(mean((colSums(fixedperid_el1)/767730079))*1000,sd((colSums(fixedperid_el1)/767730079))*1000)
         [,1]      [,2]
[1,] 0.880544 0.0156905
> cbind(mean((colSums(fixedperid_em1)/767730079))*1000,sd((colSums(fixedperid_em1)/767730079))*1000)
          [,1]       [,2]
[1,] 0.8756855 0.01364289
> sum(fixedperid_wc1/767730079)*1000
[1] 0.7270941

#-----------------------------------------------------------------------------------------------------------------------
