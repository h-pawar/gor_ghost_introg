#Mon 27 Mar 2023 15:16:01 CEST

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
test$chr<-paste('chr',i,sep='')
vfout[[i]]<-test
}
vfoutall<-do.call(rbind,vfout)
# extract top subsets of volcanofinder scores
#  > 95% outliers
vf_95<-vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.95)),]
vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))



#-----------------------------------------------------------------------------------------------------------------------

# read in annotated gorilla regulatory elements (generated in garcia-perez et al, & lifted over to hg19 by PEC)
gor_reg=read.table("/scratch/devel/pesteller/gor_reg_4Harvi/Gorilla.re_annotation.with_Non-RE.hg19.bed")
gor_reg<-as.data.frame(gor_reg)

#-----------------------------------------------------------------------------------------------------------------------


# sE regions
gor_sereg<-gor_reg[(which(gor_reg$V7=="sE")),]

gor_seregranges<-GRanges(seqnames=gor_sereg[,1],ranges=IRanges(start=as.numeric(gor_sereg[,2]),end=as.numeric(gor_sereg[,3]),names=gor_sereg[,1]),strand=rep("*",length(gor_sereg[,1])),re=(gor_sereg[,7]))

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
}
return(intersect_id)
}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 2) subset introgressed regions by the adaptive introgressed regions

overlap_vfoutl_skov3<-function(vf_obj,sk_obj) {
out_overlaps_95<-list()
out_1<-list()
for (i in (1:length(sk_obj))) {
out_overlaps_95[[i]]<-(reduce(subsetByOverlaps(sk_obj[[i]],vf_95_ranges)))
}
return(out_overlaps_95)
}
vfov_gbb2<-overlap_vfoutl_skov3(vf_95_ranges,ov_gbb_99)

vmg_se2<-intersect_withgenes_function(vfov_gbb2)
vmg_se_all2<-reduce(unlist(GRangesList(vmg_se2)))


rtes2<-findOverlaps(gor_seregranges1,vmg_se_all2)

y12<-unique(as.data.frame(rtes2)[,1])



vtest2<-data.frame(gor_seregranges1[y12])


 nrow((vtest2))
 nrow(na.omit(vtest2))

 length(unique(na.omit(vtest2)$gene))

# intersect this list with the orthologs 1:1 file from PEC *

gorids_genes<-unique(na.omit(vtest2)$gene)


gorids_genes<-unique(unlist(strsplit(gorids_genes,split=",")))

#-----------------------------------------------------------------------------------------------------------------------

# 3) subset these genes further by those that have 1:1 ortholog with humans
orth_reg=read.table("/scratch/devel/djuan/LCL_projectData/txt_files/9936_one_to_one_orthologous_autosomal_protein_coding_genes.txt")
orth_reg<-as.data.frame(orth_reg)

#-----------------------------------------------------------------------------------------------------------------------

hold<-list()
for (i in (1:length(gorids_genes))) {
hold[[i]]<-orth_reg[which(orth_reg$Gorilla==gorids_genes[[i]]),c(1,3)]
}

gene_hg<-do.call(rbind,hold)

nrow(gene_hg)

length(hold)


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

mgse_h <-as.data.frame(mgse_hg[,1])

write.table(mgse_h,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/mg.se.genes.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

#-----------------------------------------------------------------------------------------------------------------------

# write out the background set - genes regulated by sE in gorillas


gor_sereg1<-gor_sereg[,c(1:3,7,9)]


length(unique(unlist(strsplit(unique(na.omit(gor_sereg1)$V9),split=","))))


backgroundse_genes<-(unique(unlist(strsplit(unique(na.omit(gor_sereg1)$V9),split=","))))


hold2<-list()
for (i in (1:length(backgroundse_genes))) {
hold2[[i]]<-orth_reg[which(orth_reg$Gorilla==backgroundse_genes[[i]]),c(1,3)]
}

backgroundse_hg<-do.call(rbind,hold2)


#> nrow(backgroundse_hg)


backgroundse_h <-as.data.frame(backgroundse_hg[,1])

write.table(backgroundse_h,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/background.se.genes.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

#-----------------------------------------------------------------------------------------------------------------------

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
