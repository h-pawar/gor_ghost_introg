#Mon 27 Mar 2023 15:16:01 CEST
#MK call 
#overlap of vf-introg-sE : something in definition of the intersection that is going wrong → 
#try to retain not just the pure overlap, but the entire region (ie retain the whole overlapping region)
# - a merging 


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# Mon 27 Mar 2023 15:16:01 CEST
# AMEND BELOW - try taking the union of vf with introg regions -> rather than the overlapping

#union(g, g2)
#is the union of the coordinates in g and g2 - https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html


# Fri 24 Mar 2023 11:06:06 CET
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


# extract sE regions in adaptive introgressed regions of MG 
#-----------------------------------------------------------------------------------------------------------------------
# empirical overlap b/n sstar-skov
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------

# 1.1) read in volcanofinder output & extract outliers
vfout<-list()
for (i in 1:22) {
test<-read.table(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/4_",i,"_test",sep=""),sep="\t",header=T)
# chromosomes need to be specified in the same way in both data objects (the skov & the volcanofinder)
test$chr<-paste('chr',i,sep='')
vfout[[i]]<-test
}
vfoutall<-do.call(rbind,vfout)
# extract top subsets of volcanofinder scores
#  > 95% outliers
vf_95<-vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.95)),]
vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))



#-----------------------------------------------------------------------------------------------------------------------
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

gor_seregranges<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]))

#> nrow(gor_sereg)
#[1] 16311
#> nrow(gor_reg)
#[1] 68441
#> nrow(gor_sereg)/nrow(gor_reg)
#[1] 0.2383221
# with genes annotated
gor_seregranges1<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]),gene=(gor_sereg[,9]))

#-----------------------------------------------------------------------------------------------------------------------
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


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# here replicates results of mg_se_recalc_genes_16mar23.R

# taking the union of the volcanofinder outliers & the introgressed regions -> result is keep the introgressed windows
# compare to how many genes identified if perform the length cutoff


overlap_vfoutl_skov3<-function(vf_obj,sk_obj) {
out_overlaps_95<-list()
out_1<-list()
for (i in (1:length(sk_obj))) {
# overlap the volcanofinder quantile with the skov fragments per individual
#out_overlaps_95[[i]]<-subsetByOverlaps(vf_obj,sk_obj[[i]]) 
# output the union of coords rather than the subsbyoverlaps
out_overlaps_95[[i]]<-reduce(union(vf_obj,sk_obj[[i]]))
out_overlaps_95[[i]]<-out_overlaps_95[[i]][width(out_overlaps_95[[i]])>2]
}
return(out_overlaps_95)
}
vfov_gbb2<-overlap_vfoutl_skov3(vf_95_ranges,ov_gbb_99)


#[[12]]
#GRanges object with 259 ranges and 0 metadata columns:
#        seqnames              ranges strand

# as opposed to 
#vfov_gbb[[12]] - if include the snps (ie without the window size cutoff)
#GRanges object with 131612 ranges and 0 metadata columns:


vmg_se2<-intersect_withgenes_function(vfov_gbb2)
vmg_se_all2<-reduce(unlist(GRangesList(vmg_se2)))


rtes2<-findOverlaps(gor_seregranges1,vmg_se_all2)

y12<-unique(as.data.frame(rtes2)[,1])



vtest2<-data.frame(gor_seregranges1[y12])



 nrow((vtest2))
nrow(na.omit(vtest2))

#>  nrow((vtest2))
#[1] 1007
#> nrow(na.omit(vtest2))
#[1] 786

# length(unique(na.omit(vtest2)$gene))
#[1] 384

# ok now have reached the same val as previouslys ** -> ie this number of genes in the introg regions as a whole
  # but has not reduced the number of candidates, b/c is not using the volcanofinder info...

# ie replicates the results of  /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/revisions_feb23/mg_se_recalc_genes_16mar23.R  


unique(unlist(strsplit(unique(na.omit(vtest2)$gene),split=",")))

length(unique(unlist(strsplit(unique(na.omit(vtest2)$gene),split=","))))
#[1] 389

# 389 - b/c some elements assigned to multiple genes


mgse_genes<-unique(unlist(strsplit(unique(na.omit(vtest2)$gene),split=",")))


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2) subset introgressed regions by the adaptive introgressed regions


#subsetByOverlaps(vf_95_ranges,ov_gbb_99[[1]]) 


# width(reduce(subsetByOverlaps(ov_gbb_99[[1]],vf_95_ranges)))
# [1] 180001  60001  62001 183001 174001  57001 386001  96001  83001  60001
#[11] 296001 243001  48001  46001 194001  53001 154001 160001 116001 176001
#[21] 175001 172001 111001  52001 187001  31001 288001 239001  64001  95001
#[31] 119001
# this way around...


overlap_vfoutl_skov3<-function(vf_obj,sk_obj) {
out_overlaps_95<-list()
out_1<-list()
for (i in (1:length(sk_obj))) {
# overlap the volcanofinder quantile with the skov fragments per individual
out_overlaps_95[[i]]<-(reduce(subsetByOverlaps(sk_obj[[i]],vf_95_ranges)))
}
return(out_overlaps_95)
}
vfov_gbb2<-overlap_vfoutl_skov3(vf_95_ranges,ov_gbb_99)


#[[12]]
#GRanges object with 259 ranges and 0 metadata columns:
#        seqnames              ranges strand

# as opposed to 
#vfov_gbb[[12]]
#GRanges object with 131612 ranges and 0 metadata columns:


vmg_se2<-intersect_withgenes_function(vfov_gbb2)
vmg_se_all2<-reduce(unlist(GRangesList(vmg_se2)))


rtes2<-findOverlaps(gor_seregranges1,vmg_se_all2)

y12<-unique(as.data.frame(rtes2)[,1])



vtest2<-data.frame(gor_seregranges1[y12])


 nrow((vtest2))
 nrow(na.omit(vtest2))

 length(unique(na.omit(vtest2)$gene))

#>  nrow((vtest2))
#[1] 207
#>  nrow(na.omit(vtest2))
#[1] 155

#>  length(unique(na.omit(vtest2)$gene))
#[1] 72

# ok this has now reduced the number of genes 
# intersect this list with the orthologs 1:1 file from PEC *

gorids_genes<-unique(na.omit(vtest2)$gene)


gorids_genes<-unique(unlist(strsplit(gorids_genes,split=",")))

#length(gorids_genes)
#[1] 73
#-----------------------------------------------------------------------------------------------------------------------

# 3) subset these genes further by those that have 1:1 ortholog with humans

# extract human identifers - Fri 24 Mar 2023 16:20:15 CET
#/scratch/devel/djuan/LCL_projecn identiferstData/txt_files/9936_one_to_one_orthologous_autosomal_protein_coding_genes.txt 

# subset by these genes -> huma

orth_reg=read.table("/scratch/devel/djuan/LCL_projectData/txt_files/9936_one_to_one_orthologous_autosomal_protein_coding_genes.txt")


orth_reg<-as.data.frame(orth_reg)

#-----------------------------------------------------------------------------------------------------------------------

#> orth_reg[which(orth_reg$Gorilla==gorids_genes[[1]]),c(1,3)]
#              Human            Gorilla
#448 ENSG00000117151 ENSGGOG00000015738



hold<-list()
for (i in (1:length(gorids_genes))) {
hold[[i]]<-orth_reg[which(orth_reg$Gorilla==gorids_genes[[i]]),c(1,3)]
}

gene_hg<-do.call(rbind,hold)

#nrow(gene_hg)
#[1] 45
#> length(hold)
#[1] 73
# now down to 45 candidate genes - these woudl be the candidates for gene set enrichment
  # these are the genes regulated by sEs in adaptive introg regions 

# write out this file
# also write out the background *


gene_h <-as.data.frame(gene_hg[,1])

# write out the MG adaptive introg sE genes (genes regulated by sEs in adaptive introg regions)
write.table(gene_h,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/mg.adap.se.genes.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

#-----------------------------------------------------------------------------------------------------------------------

# write out the MG sE genes  (genes regulated by sEs in introg regions)

hold1<-list()
for (i in (1:length(mgse_genes))) {
hold1[[i]]<-orth_reg[which(orth_reg$Gorilla==mgse_genes[[i]]),c(1,3)]
}

mgse_hg<-do.call(rbind,hold1)

#nrow(mgse_hg)
#[1] 235
#length(hold1)
#[1] 389


# write out this file
# also write out the background *


mgse_h <-as.data.frame(mgse_hg[,1])

write.table(mgse_h,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/mg.se.genes.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

#-----------------------------------------------------------------------------------------------------------------------

# write out the background set - genes regulated by gor sEs


#-----------------------------------------------------------------------------------------------------------------------

# write out the background - genes regulated by sE in gorillas

# head(gor_sereg)
#     V1        V2        V3     V4 V5 V6 V7   V8   V9
#3  chr1 223189674 223190744  re_n4 sE sE sE <NA> <NA>
#4  chr1 223111236 223111830 re_n14 sE sE sE <NA> <NA>

gor_sereg1<-gor_sereg[,c(1:3,7,9)]

#na.omit(gor_sereg1$V9)


#>  nrow((gor_sereg1))
#[1] 16311
#>  nrow(na.omit(gor_sereg1))
#[1] 13135


length(unique(unlist(strsplit(unique(na.omit(gor_sereg1)$V9),split=","))))
#[1] 5523

backgroundse_genes<-(unique(unlist(strsplit(unique(na.omit(gor_sereg1)$V9),split=","))))


hold2<-list()
for (i in (1:length(backgroundse_genes))) {
hold2[[i]]<-orth_reg[which(orth_reg$Gorilla==backgroundse_genes[[i]]),c(1,3)]
}

backgroundse_hg<-do.call(rbind,hold2)


#> nrow(backgroundse_hg)
#[1] 3010
 
#> length(hold2)
#[1] 5523



backgroundse_h <-as.data.frame(backgroundse_hg[,1])

write.table(backgroundse_h,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/background.se.genes.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

#-----------------------------------------------------------------------------------------------------------------------





#-----------------------------------------------------------------------------------------------------------------------

> gene_hg
               Human            Gorilla
448  ENSG00000117151 ENSGGOG00000015738
447  ENSG00000117133 ENSGGOG00000015862
450  ENSG00000171517 ENSGGOG00000011067
5786 ENSG00000173848 ENSGGOG00000006693
5820 ENSG00000148484 ENSGGOG00000028048
5996 ENSG00000107798 ENSGGOG00000003979
6036 ENSG00000187122 ENSGGOG00000005553
6038 ENSG00000052749 ENSGGOG00000014850
6597 ENSG00000070961 ENSGGOG00000014187
6620 ENSG00000111145 ENSGGOG00000009227
6621 ENSG00000059758 ENSGGOG00000001172
6650 ENSG00000171310 ENSGGOG00000016823
7153 ENSG00000054654 ENSGGOG00000002988
7655 ENSG00000068305 ENSGGOG00000013064
7801 ENSG00000166501 ENSGGOG00000010971
8023 ENSG00000166446 ENSGGOG00000006084
3075 ENSG00000039123 ENSGGOG00000001867
3065 ENSG00000213949 ENSGGOG00000022320
3063 ENSG00000151883 ENSGGOG00000001459
3057 ENSG00000112972 ENSGGOG00000001140
3058 ENSG00000151882 ENSGGOG00000012830
3044 ENSG00000132356 ENSGGOG00000027961
3016 ENSG00000056097 ENSGGOG00000016302
3015 ENSG00000150712 ENSGGOG00000009387
8744 ENSG00000175387 ENSGGOG00000004706
9350 ENSG00000105186 ENSGGOG00000008596
2180 ENSG00000185008 ENSGGOG00000011779
2681 ENSG00000083896 ENSGGOG00000021958
2874 ENSG00000109445 ENSGGOG00000006078
2892 ENSG00000170390 ENSGGOG00000014019
8373 ENSG00000126353 ENSGGOG00000016340
3201 ENSG00000112893 ENSGGOG00000001547
3237 ENSG00000072682 ENSGGOG00000006366
3457 ENSG00000113300 ENSGGOG00000003675
3706 ENSG00000112659 ENSGGOG00000005599
3707 ENSG00000146216 ENSGGOG00000000442
3710 ENSG00000171467 ENSGGOG00000012851
3711 ENSG00000124574 ENSGGOG00000026443
3717 ENSG00000112715 ENSGGOG00000012714
3715 ENSG00000172432 ENSGGOG00000014744
3788 ENSG00000112701 ENSGGOG00000006460
3967 ENSG00000164442 ENSGGOG00000016096
4008 ENSG00000130340 ENSGGOG00000003616
4009 ENSG00000078269 ENSGGOG00000012210
4530 ENSG00000047249 ENSGGOG00000000515


#-----------------------------------------------------------------------------------------------------------------------

# this is the object with sEs in adaptive introg MG regions & corresponding genes
#> head(na.omit(vtest2))
#   seqnames    start      end width strand re
#1      chr1 86922377 86924377  2001      * sE
##2      chr1 86926569 86929182  2614      * sE
#7      chr1 87236867 87239456  2590      * sE
#9     chr10  5562517  5564519  2003      * sE
#10    chr10  5571357  5579645  8289      * sE
#11    chr10 17264352 17268358  4007      * sE
#                                                       gene
#1  ENSGGOG00000015738,ENSGGOG00000015862,ENSGGOG00000026782
#2                                        ENSGGOG00000015738
#7                                        ENSGGOG00000011067
#9                                        ENSGGOG00000006693
#10                                       ENSGGOG00000041078
#11                                       ENSGGOG00000028048
#> nrow(na.omit(vtest2))
#[1] 155

#gene_hg
#               Human            Gorilla
#448  ENSG00000117151 ENSGGOG00000015738
#447  ENSG00000117133 ENSGGOG00000015862


f<-na.omit(vtest2)



hold_genes_se<-list()
for (i in (1:nrow(gene_hg))) {
hold_genes_se[[i]]<- f[which(f$gene==gene_hg[i,2]),]}

do.call(rbind,hold_genes_se)


#> do.call(rbind,hold_genes_se)
#    seqnames     start       end width strand re               gene
#2       chr1  86926569  86929182  2614      * sE ENSGGOG00000015738
#7       chr1  87236867  87239456  2590      * sE ENSGGOG00000011067
#9      chr10   5562517   5564519  2003      * sE ENSGGOG00000006693
#11     chr10  17264352  17268358  4007      * sE ENSGGOG00000028048
#17     chr10 120014059 120019267  5209      * sE ENSGGOG00000003979
#30     chr10 128054745 128056946  2202      * sE ENSGGOG00000005553

mggenes_sedf<-do.call(rbind,hold_genes_se)

# which sE corresponds to which gene 

mggenes_se<-GRanges(seqnames=mggenes_sedf[,1],ranges=IRanges(start=as.numeric(mggenes_sedf[,2]),end=as.numeric(mggenes_sedf[,3]),names=mggenes_sedf[,1]),strand=rep("*",length(mggenes_sedf[,1])),re=(mggenes_sedf[,6]),gene=(mggenes_sedf[,7]))
write.table(mggenes_sedf,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/mggenes_se.txt",sep=""),sep="\t",row.names=F,col.names=T,quote=F)

#-----------------------------------------------------------------------------------------------------------------------
