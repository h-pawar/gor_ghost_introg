#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
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
test$chr<-paste('chr',i,sep='')
vfout[[i]]<-test
}
vfoutall<-do.call(rbind,vfout)
# extract top subsets of volcanofinder scores
#  > 95% outliers
vf_95<-vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = 0.95)),]

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))

vfall_ranges<-GRanges(seqnames=vfoutall[,5],ranges=IRanges(start=as.numeric(vfoutall[,1]),end=as.numeric(vfoutall[,1]+1),names=vfoutall[,5]),strand=rep("*",length(vfoutall[,1])),LR=as.numeric(vfoutall[,2]))

#-----------------------------------------------------------------------------------------------------------------------

genes<-do.call(c,e_cand_fin)

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

# output is in bed graph format - plot these
bedgraph_1<-vf_score_pergene(1)
#-----------------------------------------------------------------------------------------------------------------------

bedgraph_2<-vf_score_pergene(2)


x1<-ggplot(bedgraph_1, aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 5')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("SEMA5A")
x2<-ggplot(bedgraph_2, aes(x=start, y=LR)) + geom_line( color="#a5af76", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 12')  + geom_hline(yintercept=23.628, color="#76a5af") + ggtitle("TAS2R14")


pdf(paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/sema5a.tas2r14.22aug.1.pdf"))

x1
x2
ggarrange(x1,x2, ncol = 2, nrow = 1)

dev.off()


#-----------------------------------------------------------------------------------------------------------------------

out_beds<-list()
for (ind in (1:length(e_cand_fin))) {
out_beds[[ind]]<-vf_score_pergene(ind)
}



#-----------------------------------------------------------------------------------------------------------------------

# arrange in one panel

# Tue 27 Sep 2022 11:40:50 CEST
# arrange in one panel

a<-ggplot(out_beds[[1]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 5')  + geom_hline(yintercept=23.628, color="red")  + ggtitle("SEMA5A") + theme(axis.text.x=element_text(size=rel(0.5)))
b<-ggplot(out_beds[[2]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 12')  + geom_hline(yintercept=23.628, color="red") + ggtitle("TAS2R14") + theme(axis.text.x=element_text(size=rel(0.5)))
c<-ggplot(out_beds[[3]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 3')  + geom_hline(yintercept=23.628, color="red") + ggtitle("KALRN") + theme(axis.text.x=element_text(size=rel(0.5)))
d<-ggplot(out_beds[[4]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 4')  + geom_hline(yintercept=23.628, color="red") + ggtitle("TMPRSS11E") + theme(axis.text.x=element_text(size=rel(0.5)))
e<-ggplot(out_beds[[5]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 1')  + geom_hline(yintercept=23.628, color="red") + ggtitle("RP11-541F9.2") + theme(axis.text.x=element_text(size=rel(0.5)))
f<-ggplot(out_beds[[6]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 10')  + geom_hline(yintercept=23.628, color="red") + ggtitle("GPAM") + theme(axis.text.x=element_text(size=rel(0.5)))
g<-ggplot(out_beds[[7]], aes(x=start, y=LR)) + geom_line( color="black", size=2, alpha=0.9) + theme_classic() + xlab('Chromosome 7')  + geom_hline(yintercept=23.628, color="red") + ggtitle("HOXA3/4/5") + theme(axis.text.x=element_text(size=rel(0.5)))


pdf("/scratch/devel/hpawar/admix/volcanofinder/output/merged/7candidategenes.25aug.rearranged.pdf")
ggarrange(b,a,c,d,e,f,g, 
          ncol = 3, nrow = 3)
dev.off()

