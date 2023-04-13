# Mon 20 Mar 2023 10:15:33 CET
# focus on the sEs - & check if they are genic enhancers/promoter-interacting/enhancer-interacting 

#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
#library(pheatmap)
#library(ggplot2)
#library(ggpubr)
library(valr)
library(GenomicFeatures)


# extract sE regions in introgressed regions of MG & intersect with VEP *
# so the objects of reg bp per reg element
#load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/emp_reg_gbb_8feb23"),verbose=T)
# previous script output - population level estimates of mean regulatory bp total per element type
# to get back to which regions are associated with sE element type - will need to run interactively gor_regelements_8feb23.R
# directly query MG regions for sE regions 

#-----------------------------------------------------------------------------------------------------------------------
# empirical overlap b/n sstar-skov
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# read in annotated gorilla regulatory elements (generated in garcia-perez et al, & lifted over to hg19 by PEC)

#awk '{print $7}' /scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed | sort | uniq


gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")


#head(gor_reg)
#    V1        V2        V3     V4 V5 V6 V7   V8                 V9
#1 chr4 131391181 131391381  re_n2 wE wE wE  prE ENSGGOG00000012466
#2 chr4 131357672 131358072  re_n3 wE wE wE <NA>               <NA>

gor_reg<-as.data.frame(gor_reg)

#-----------------------------------------------------------------------------------------------------------------------


# sE regions
gor_sereg<-gor_reg[(which(gor_reg$V7=="sE")),]

# add column of genic/promoter interacting etc
gor_seregranges<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]),pos=(gor_sereg[,8]))


#-----------------------------------------------------------------------------------------------------------------------
# only output the intersect regions
intersect_withgenes_function<-function(empirical_id) {
intersect_id<-list()
counts_id<-list()
for (ind in (1:length(empirical_id))) {
intersect_id[[ind]]<-intersect(gor_seregranges, empirical_id[[ind]], ignore.strand=TRUE)
# output - total number of introgressed bp in a regulatory element;  total introgressed bp in this individual
#counts_id[[ind]]<-cbind(sum(reduce(intersect_id[[ind]])@ranges@width), sum(empirical_id[[ind]]@ranges@width))
}
#return(list(intersect_id,do.call(rbind,counts_id)))
return(intersect_id)
}

mg_se<-intersect_withgenes_function(ov_gbb_99)

# as granges list -> reduce

#> reduce(unlist(GRangesList(mg_se)))
#GRanges object with 1007 ranges and 0 metadata columns:
#         seqnames              ranges strand
#            <Rle>           <IRanges>  <Rle>
 #    [1]     chr1     5447068-5448468      *
  #   [2]     chr1     5954130-5954542      *
  #   [3]     chr1     5971320-5976823      *
  #   [4]     chr1     9865463-9867007      *
  #   [5]     chr1     9874404-9875411      *
  #   ...      ...                 ...    ...
  #[1003]     chr7 147373964-147376153      *
  #[1004]     chr7 147387743-147388743      *
  #[1005]     chr7 147417519-147418325      *
  #[1006]     chr7 147420916-147425364      *
  #[1007]     chr7 147425764-147427168      *
  #-------
  #seqinfo: 34 sequences from an unspecified genome; no seqlengths


mg_se_all<-reduce(unlist(GRangesList(mg_se)))

# from this go back to which genes -> is this info also in the reg files?
# yes V9 is the associated gene -> shoudl annotate with this info * - findoverlaps 

#-----------------------------------------------------------------------------------------------------------------------

# 1) looking at the population level


rtes<-findOverlaps(gor_seregranges,mg_se_all)

y1<-unique(as.data.frame(rtes)[,1])
test<-data.frame(gor_seregranges[y1])

# then assess position of enhancers *

#head(test)
#  seqnames     start       end width strand re  pos
#1     chr1 223189674 223190744  1071      * sE <NA>
#2     chr4 131376499 131377101   603      * sE <NA>
#3     chr4 131387453 131389267  1815      * sE <NA>
#4     chr1   5447068   5448468  1401      * sE   gE
#5     chr1   5954130   5954542   413      * sE   gE

table(test$pos)
EiE  gE   P PiE prE 
 20 547  50  89  80 

 nrow(test)
[1] 1007

sum(table(test$pos))
[1] 786

 sum(table(test$pos))/ nrow(test)
[1] 0.7805362

# minority of MG sEs in introg regions don't have corresponding position annotated

# explanation of the abbreviations from https://www.nature.com/articles/s41467-021-23397-1#Sec12
#EiE  gE   P PiE prE 
	# promoter-interacting enhancers (PiE) and enhancer-interacting enhancers (EiE) 
	# genic promoters (gP), gE, and proximal enhancers (prE).
	# gE = intragenic enhancer


# p = just a promoter??

# > test[which(test$pos=="P"),]
#     seqnames     start       end width strand re pos
#8        chr1   9874404   9875411  1008      * sE   P
#13       chr1  10674973  10678395  3423      * sE   P
#21       chr1  36405648  36407603  1956      * sE   P
#31       chr1  37375329  37381130  5802      * sE   P

# Q - a strong enhancer inside a promoter?? - ask PEC re this * (yes is possible)

# further question - I don't have a good expectation of what this distribution looks like for the rest of the genome...
	# from fig 4b - seems like gE are the majority of annotated elements 
		# 4b Average number of RE across species associated with one-to-one orthologous protein-coding genes 

#-----------------------------------------------------------------------------------------------------------------------

# Q - maybe woudl be best to calculate this per individual as well, rather than per population - yes I think resolution is better per id -> then take means

# 2) position of strong enhancers per individual

proc_perid_function<-function(id) {
rtes<-findOverlaps(gor_seregranges,mg_se[[id]])
y1<-unique(as.data.frame(rtes)[,1])
test<-data.frame(gor_seregranges[y1])
out<-table(test$pos)
return(out)}


pos_perid<-list()
for (i in 1:length(mg_se)){
print(i)
pos_perid[[i]]<-proc_perid_function(i)}


> do.call(rbind,pos_perid)
      EiE  gE  P PiE prE
 [1,]   6  99  7  17   9
 [2,]   4 134 10  26  11
 [3,]   7 156 15  25  24
 [4,]   1 116 10  14  15
 [5,]   8 112 14  26  17
 [6,]   3 129  5  23  16
 [7,]   3 129  9  17  21
 [8,]   8 177 14  28  28
 [9,]   6 122  8  24  20
[10,]   4  42  2  13   6
[11,]   2  78  5  12   9
[12,]   4 113 10  17  13

# i think per individual is a better way to go

# positions of each element coudl be useful to retain
# Q - which enhancer type are more interesting per se?

x<-do.call(rbind,pos_perid)
colMeans(x)

> colMeans(x)
       EiE         gE          P        PiE        prE 
  4.666667 117.250000   9.083333  20.166667  15.750000 

