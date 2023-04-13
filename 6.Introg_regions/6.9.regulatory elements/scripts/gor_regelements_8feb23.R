# Wed  8 Feb 2023 17:25:01 CET
# I would like to check the putative introgressed regions for gorilla regulatory elements which PEC annotated in the LCLs project. 
    # published in https://www.nature.com/articles/s41467-021-23397-1
# Do you have a bed file of positions of gorilla regulatory elements available? & if so would you mind passing me the path?
    # relevant file is  Supplementary Data 1

#-----------------------------------------------------------------------------------------------------------------------
# PEC converted excel to bed file & performed liftover to hg19:     

# PEC -  Tue, 7 Feb, 16:18
#I looked for the PCA Raquel did (shown below). 
#The two samples (grey triangles: G1, G2) belong to western gorillas (most likely western lowland gorillas). 
#At some point it could be interesting to highlight that in your response. 
#Also, although this might seem obvious, we must not forget these are LCLs. 
#This means the regulatory elements we see in these cell lines will be most likely linked to immune response and they will likely change their status in other cell types.

#The coordinates for the regulatory elements are shown in each specific species reference assembly. 
#For the gorilla, they are in gorgor4. Since I have all chain files, I have done the liftover for you. 
#However, there is no chain to go from Gorgor4 to hg19 (at least, I have not seen any). 
#Because of that, the liftovers have been done like this: GorGor4 -> hg38 -> hg19. 



#-----------------------------------------------------------------------------------------------------------------------

# for gorilla annotated regulatory elements

#have asked MK

#for intersecting the introgressed regions with the gorilla annotated regulatory regions, can I check the useful output would be

#- mean of ( number of introgressed bp in a regulatory element / total introgressed bp in this individual) # equivalent to protein coding bp analysis
#- or should the output rather be introgressed base pairs per regulatory element type (weak enhancer, strong promoter etc)?

#MK

#I would suggest to calculate both, as both could be interesting. 
#Of course, the first is a direct comparison to the coding bp, so that may go first.

# -> output both measures bp in reg element/total introg bp per id
# bp per reg element / total introg bp per id 


#-----------------------------------------------------------------------------------------------------------------------
# follow intersecting procedure of  random.gene.density.s*.skov.overlap.21oct22.R
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------------------------------------------------
# empirical overlap b/n sstar-skov
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------

# read in annotated gorilla regulatory elements (generated in garcia-perez et al, & lifted over to hg19 by PEC)

#awk '{print $7}' /scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed | sort | uniq


gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")


#head(gor_reg)
#    V1        V2        V3     V4 V5 V6 V7   V8                 V9
#1 chr4 131391181 131391381  re_n2 wE wE wE  prE ENSGGOG00000012466
#2 chr4 131357672 131358072  re_n3 wE wE wE <NA>               <NA>

gor_reg<-as.data.frame(gor_reg)


# PEC: Only focus of column 7, Column 7 is the consensus RE of the SPECIES --> The one column you want

#nrow(gor_reg)
#[1] 68441
#-----------------------------------------------------------------------------------------------------------------------



# read in annotated gorilla regulatory elements (generated in garcia-perez et al, & lifted over to hg19 by PEC)

#awk '{print $7}' /scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed | sort | uniq


gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")


#head(gor_reg)
#    V1        V2        V3     V4 V5 V6 V7   V8                 V9
#1 chr4 131391181 131391381  re_n2 wE wE wE  prE ENSGGOG00000012466
#2 chr4 131357672 131358072  re_n3 wE wE wE <NA>               <NA>

gor_reg<-as.data.frame(gor_reg)


# PEC: Only focus of column 7, Column 7 is the consensus RE of the SPECIES --> The one column you want

#nrow(gor_reg)
#[1] 68441
#-----------------------------------------------------------------------------------------------------------------------

# retain only those with annotated regulatory function in both replicates *
gor_allreg<-gor_reg[-(which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re"| gor_reg$V7=="Non-re")),]

gor_allregranges<-GRanges(seqnames=gor_allreg[,1],ranges=IRanges(start=as.numeric(gor_allreg[,2]),end=as.numeric(gor_allreg[,3]),names=gor_allreg[,1]),strand=rep("*",length(gor_allreg[,1])),re=(gor_allreg[,7]))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

intersect_withgenes_function<-function(empirical_id) {
intersect_id<-list()
counts_id<-list()
for (ind in (1:length(empirical_id))) {
intersect_id[[ind]]<-intersect(gor_allregranges, empirical_id[[ind]], ignore.strand=TRUE)
# output - total number of introgressed bp in a regulatory element;  total introgressed bp in this individual
counts_id[[ind]]<-cbind(sum(reduce(intersect_id[[ind]])@ranges@width), sum(empirical_id[[ind]]@ranges@width))
}
return(list(intersect_id,do.call(rbind,counts_id)))
}



# for MG

test<-intersect_withgenes_function(ov_gbb_99)

#test
#[[1]][[12]]
#GRanges object with 806 ranges and 0 metadata columns:
#        seqnames            ranges strand
#           <Rle>         <IRanges>  <Rle>
#    [1]     chr4   6883642-6883856      *
#    [2]     chr4   6891501-6891701      *
#    [3]     chr4   6894132-6896317      *
#    [4]     chr4   6910751-6911951      *
#    [5]     chr4   6992895-6997501      *
#    ...      ...               ...    ...
#  [802]    chr15 58649504-58649901      *
#  [803]    chr15 58655110-58668203      *
#  [804]    chr15 58669007-58669211      *
#  [805]    chr19 29485758-29486959      *
#  [806]    chr19 51623395-51624794      *
#  -------
#  seqinfo: 36 sequences from an unspecified genome; no seqlengths

#[[2]]
#         [,1]     [,2]
# [1,] 1478861 28826232
# [2,] 1766000 34622287


# take the ratios of all -> # then mean & sd - pop-level

# test[[2]][,1]/test[[2]][,2]
# [1] 0.05130261 0.05100761 0.05723155 0.04869737 0.05450021 0.05216897
# [7] 0.06100206 0.06212825 0.05190936 0.04587091 0.04490061 0.05240231

#- mean of ( number of introgressed bp in a regulatory element / total introgressed bp in this individual) # equivalent to protein coding bp analysis


# pop level estimate for MG: introg bp in regulatory element - normalised by total introg bp in this individual
mean(test[[2]][,1]/test[[2]][,2])
#[1] 0.05276015
sd(test[[2]][,1]/test[[2]][,2])
#[1] 0.005315919

#-----------------------------------------------------------------------------------------------------------------------


eltest<-intersect_withgenes_function(ov_gbg_99)

#(eltest[[2]][,1]/eltest[[2]][,2])
#[1] 0.04439733 0.04810593 0.03559400 0.04937914 0.04334540 0.04980529 0.04751309
#[8] 0.04372746 0.04964324

mean(eltest[[2]][,1]/eltest[[2]][,2])
sd(eltest[[2]][,1]/eltest[[2]][,2])


# pop level estimate for EL
 mean(eltest[[2]][,1]/eltest[[2]][,2])
#[1] 0.04572343
 sd(eltest[[2]][,1]/eltest[[2]][,2])
#[1] 0.004586294


#-----------------------------------------------------------------------------------------------------------------------
mg_r<-cbind(mean(test[[2]][,1]/test[[2]][,2]),sd(test[[2]][,1]/test[[2]][,2]))
el_r<-cbind(mean(eltest[[2]][,1]/eltest[[2]][,2]),sd(eltest[[2]][,1]/eltest[[2]][,2]))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# introg bp per reg element
# need to annotate overlap with which reg elements identified

# input = test (MG)
# input = eltest (EL) 
# id -> for each id of each scenario

bp_perreg_fun<-function(input,id){

rtes<-findOverlaps(gor_allregranges,input[[1]][[id]])

y1<-unique(as.data.frame(rtes)[,1])

# extract mutations from vep data within range of overlap with the candidate genic region
out_re<-list()
for (ind in (1:length(y1))) {
out_re[[ind]]<-data.frame(gor_allregranges[y1[[ind]]])
}

reg_info_1id<-do.call(rbind,out_re)

paths = by(reg_info_1id, reg_info_1id[,"re"], function(x) x)

counts_perel<-list()
for (ind in (1:length(paths))) {
counts_perel[[ind]]<-cbind(sum(paths[[ind]][,4]-1),  unique(paths[[ind]][[6]]))}


p<-as.data.frame(do.call(rbind, counts_perel))

p[,1]<-as.numeric(p[,1])

p[,3]<-p[,1]/input[[2]][id,2]

return(p[,c(2:3)])
}


#bp_perreg_fun(test,1)# works fine

# run for all individuals per scenario
pop_perreg_fun<-function(input){
o_per<-list()
for (ind in (1:nrow(input[[2]]))) {
o_per[[ind]]<-bp_perreg_fun(input,ind)}
return(o_per)
}


p_mg<-pop_perreg_fun(test)
p_el<-pop_perreg_fun(eltest)

# then take pop-wise means of introg bp per reg element 

process_fun<-function(input){

p1_o<- do.call(cbind, input)[,c(1,seq(2,length(input)*2,2))]
p2_o<-as.data.frame(lapply((p1_o)[,-1],as.numeric))

q_per<-list()
for (ind in (1:nrow(p2_o))) {
q_per[[ind]]<- mean(as.numeric(p2_o[ind,]))}

return(as.data.frame(cbind(p1_o[,1],do.call(rbind,q_per))))
}

process_fun(p_mg)
process_fun(p_el)


 process_fun(p_mg)
   V1                   V2
1  aE  0.00600013520787739
2  aP   0.0012663125426476
3  pE 0.000275917360849282
4 P/E  0.00507354711827817
5  pP  0.00101191812786651
6  sE   0.0218697984960251
7  sP   0.0102888069400528
8  wE  0.00831010278360753
9  wP 0.000574572581698412

process_fun(p_el)
   V1                   V2
1  aE  0.00432584258821753
2  aP   0.0011042879024019
3  pE 0.000327911222412902
4 P/E  0.00500449830916939
5  pP 0.000547419315981555
6  sE   0.0166462621732876
7  sP   0.0105582126654595
8  wE  0.00806234428741233
9  wP 0.000550782213115432

# now calculate equivalent for random regions -> first generate random regions of sufficient callability

#-----------------------------------------------------------------------------------------------------------------------

mg_p <- process_fun(p_mg)
el_p <- process_fun(p_el)

emp_reg_mg<-list(mg_r,mg_p)
emp_reg_el<-list(el_r,el_p)

save(emp_reg_mg,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/emp_reg_gbb_8feb23"))
save(emp_reg_el,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/emp_reg_gbg_8feb23"))

#-----------------------------------------------------------------------------------------------------------------------


q()




# testing for individual 1

#ind=1
#empirical_id=ov_gbb_99
#tesintersect_id<-intersect(gor_allregranges, ov_gbb_99[[1]], ignore.strand=TRUE)
# number of introgressed bp in a regulatory element / total introgressed bp in this individual
#tescounts_id<-cbind(sum(reduce(tesintersect_id)@ranges@width), sum(ov_gbb_99[[1]]@ranges@width))
# number of introgressed bp per regulatory element / total introgressed bp in this individual

#> tesintersect_id
#GRanges object with 687 ranges and 0 metadata columns:
#        seqnames            ranges strand
#           <Rle>         <IRanges>  <Rle>
#    [1]     chr4 36285070-36285671      *
#    [2]     chr4 36430394-36433196      *
#    [3]     chr4 38652000-38652381      *
#    [4]     chr4 38656420-38659608      *
#    [5]     chr4 38700574-38700969      *
#    ...      ...               ...    ...
#  [683]    chr19   6091664-6097853      *
#  [684]    chr19   6098449-6101247      *
#  [685]    chr19   6103448-6104000      *
#  [686]    chr19 29485758-29486959      *
#  [687]    chr22 18097412-18098789      *
#  -------
#  seqinfo: 36 sequences from an unspecified genome; no seqlengths
#> tescounts_id
#        [,1]     [,2]
#[1,] 1478861 28826232

#str(ov_gbb_99)
#List of 12



# E/Non-re, P/Non-re =  an enhancer(or promoter) in one replicate, but non-regulatory in the other replicate?
# to be conservative filter out P/Non-re & E/Non-re

# head(gor_reg[which(gor_reg$V7=="P/Non-re"),])
#which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re")


nrow(gor_reg[which(gor_reg$V7=="P/Non-re"),])
[1] 649

nrow(gor_reg[which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re"),])
[1] 3779

gor_reg_filt<-gor_reg[-(which(gor_reg$V7=="P/Non-re" | gor_reg$V7=="E/Non-re")),]

nrow(gor_reg_filt)
[1] 64662

 unique(gor_reg_filt$V7)
 [1] "wE"     "sE"     "aE"     "wP"     "P/E"    "aP"     "sP"     "pE"    
 [9] "Non-re" "pP"
