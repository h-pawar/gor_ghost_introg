#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
#library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
library(pheatmap)
library(ggplot2)
library(ggpubr)

load("/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.21jul22",verbose=T)


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

#> head(vf_95)
#     location       LR        alpha        D  chr
#1573  1634282 24.36087 4.814070e-05 0.019332 chr1
#1574  1635282 25.34131 5.027065e-05 0.019332 chr1
# marker every 10kb

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 1.2) convert vf outliers to genomic ranges format 

vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))

# convert all vf windows to granges format - want to know all LR scores for the candidate regions
vfall_ranges<-GRanges(seqnames=vfoutall[,5],ranges=IRanges(start=as.numeric(vfoutall[,1]),end=as.numeric(vfoutall[,1]+1),names=vfoutall[,5]),strand=rep("*",length(vfoutall[,1])),LR=as.numeric(vfoutall[,2]))

#-----------------------------------------------------------------------------------------------------------------------

#> do.call(c,e_cand_fin)
#GRanges object with 7 ranges and 1 metadata column:
#                  seqnames              ranges strand |         gene_id
#                     <Rle>           <IRanges>  <Rle> |     <character>
#  ENSG00000112902     chr5     9035138-9546187      - | ENSG00000112902 # SEMA5A
#  ENSG00000212127    chr12   11090005-11324172      - | ENSG00000212127 # TAS2R14
#  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145
#  ENSG00000087128     chr4   69313167-69363322      + | ENSG00000087128
#  ENSG00000228215     chr1 191844625-191980390      + | ENSG00000228215
#  ENSG00000119927    chr10 113909624-113975135      - | ENSG00000119927
#  ENSG00000105997     chr7   27145803-27192200      - | ENSG00000105997
#  -------
#  seqinfo: 265 sequences from an unspecified genome; no seqlengths

#-----------------------------------------------------------------------------------------------------------------------

# works - but better splitting by gene
genes<-do.call(c,e_cand_fin)

#-----------------------------------------------------------------------------------------------------------------------


# extract vf CLR scores (per gene) - split by genes
# plot these scores across the genomic regions


# extract volcanofinder scores per gene
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# extract volcanofinder scores per gene, in format ready to plot

vf_score_pergene<-function(ind){
x<-findOverlaps(genes[ind],vfall_ranges)
x1<-unique(as.data.frame(x)[,2])

out_scores<-list()
for (ind in (1:length(x1))) {
out_scores[[ind]]<-as.data.frame(vfall_ranges[x1[[ind]]])
}

y<-do.call(rbind,out_scores)
rownames(y)<-NULL
y1<-y[,c(1:3,6)]
y1$seqnames<-as.character(y1$seqnames)
return(y1)
}

# output is in bed graph format

#> x[,c(1:3,6)]
#  seqnames   start     end       LR
#1     chr5 9061766 9061767 26.75713
#2     chr5 9062766 9062767 27.97717

    # these are the LR scores of the intersect regions

# plot these
bedgraph_1<-vf_score_pergene(1)
#-----------------------------------------------------------------------------------------------------------------------

bedgraph_2<-vf_score_pergene(2)


pdf(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/sema5a.tas2r14.22aug.pdf"))

# change axes labels 
ggplot(bedgraph_1, aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 5')  + geom_hline(yintercept=23.628, color="#76a5af")
ggplot(bedgraph_2, aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 12')  + geom_hline(yintercept=23.628, color="#76a5af")

dev.off()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/volcanofinder/output/merged/sema5a.tas2r14.22aug.pdf  /Users/harvi/Downloads/gorilla_volcanofinder

# add titles of the gene names, & arrange on 

#+ ggtitle()

# SEMA5A
# TAS2R14

#+ ggtitle("TAS2R14")

x1<-ggplot(bedgraph_1, aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 5')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("SEMA5A")
x2<-ggplot(bedgraph_2, aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 12')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("TAS2R14")


pdf(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/sema5a.tas2r14.22aug.1.pdf"))

x1
x2
ggarrange(x1,x2, ncol = 2, nrow = 1)

dev.off()


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# https://docs.google.com/spreadsheets/d/1yLq5JE5OQC-yZ6E6aBDyz_4LZZwLf6PRlO0IrLcqCZo/edit#gid=1264452875

# plot rest of the candidate genes

#  ENSG00000160145     chr3 123798870-124445172      + | ENSG00000160145 # KALRN
#  ENSG00000087128     chr4   69313167-69363322      + | ENSG00000087128 # TMPRSS11E
#  ENSG00000228215     chr1 191844625-191980390      + | ENSG00000228215 # RP11-541F9.2
#  ENSG00000119927    chr10 113909624-113975135      - | ENSG00000119927 # GPAM
#  ENSG00000105997     chr7   27145803-27192200      - | ENSG00000105997 # HOXA3/4/5
#  -------


out_beds<-list()
for (ind in (1:length(e_cand_fin))) {
out_beds[[ind]]<-vf_score_pergene(ind)
}
#head(out_beds[[7]]$start)
#[1] 27146647 27147647 27148647 27149647 27150647 27151647


out_beds[[7]]

# plot these
bedgraph_1<-vf_score_pergene(1)
#-----------------------------------------------------------------------------------------------------------------------

bedgraph_2<-vf_score_pergene(2)


pdf(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/7candidategenes.25aug.pdf"))

# change axes labels 
ggplot(out_beds[[1]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 5')  + geom_hline(yintercept=23.628, color="#76a5af")  + ggtitle("SEMA5A")
ggplot(out_beds[[2]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 12')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("TAS2R14")
ggplot(out_beds[[3]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 3')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("KALRN")
ggplot(out_beds[[4]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 4')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("TMPRSS11E")
ggplot(out_beds[[5]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 1')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("RP11-541F9.2")
ggplot(out_beds[[6]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 10')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("GPAM")
ggplot(out_beds[[7]], aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 7')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("HOXA3/4/5")

dev.off()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/volcanofinder/output/merged/7candidategenes.25aug.pdf  /Users/harvi/Downloads/gorilla_volcanofinder
# figures used for the poster

#-----------------------------------------------------------------------------------------------------------------------

# arrange in one panel

# Tue 27 Sep 2022 11:40:50 CEST
# arrange in one panel

# try to make x axis text smaller + theme(axis.text.x=element_text(size=rel(0.5)))


a<-ggplot(out_beds[[1]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 5')  + geom_hline(yintercept=23.628, color="red")  + ggtitle("SEMA5A") + theme(axis.text.x=element_text(size=rel(0.5)))
b<-ggplot(out_beds[[2]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 12')  + geom_hline(yintercept=23.628, color="red") + ggtitle("TAS2R14") + theme(axis.text.x=element_text(size=rel(0.5)))
c<-ggplot(out_beds[[3]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 3')  + geom_hline(yintercept=23.628, color="red") + ggtitle("KALRN") + theme(axis.text.x=element_text(size=rel(0.5)))
d<-ggplot(out_beds[[4]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 4')  + geom_hline(yintercept=23.628, color="red") + ggtitle("TMPRSS11E") + theme(axis.text.x=element_text(size=rel(0.5)))
e<-ggplot(out_beds[[5]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 1')  + geom_hline(yintercept=23.628, color="red") + ggtitle("RP11-541F9.2") + theme(axis.text.x=element_text(size=rel(0.5)))
f<-ggplot(out_beds[[6]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 10')  + geom_hline(yintercept=23.628, color="red") + ggtitle("GPAM") + theme(axis.text.x=element_text(size=rel(0.5)))
g<-ggplot(out_beds[[7]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 7')  + geom_hline(yintercept=23.628, color="red") + ggtitle("HOXA3/4/5") + theme(axis.text.x=element_text(size=rel(0.5)))


#http://bioconductor.org/packages/devel/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
pdf("/scratch/devel/hpawar/admix/volcanofinder/output/merged/7candidategenes.25aug.rearranged.pdf")
ggarrange(b,a,c,d,e,f,g, 
          ncol = 3, nrow = 3)
dev.off()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/volcanofinder/output/merged/7candidategenes.25aug.rearranged.pdf  /Users/harvi/Downloads/gorilla_volcanofinder
# have added this to the google drive