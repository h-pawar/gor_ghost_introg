#!/usr/bin/r
#module load TABIX/0.2.6 gcc/6.3.0 xz python/2.7.11 R/3.2.0 perl BCFTOOLS/1.6 SAMTOOLS/1.0 BEDTools/2.26.0 bedops vcftools java/latest
# MK: The whole genome vcf files are here: /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz. recommend to take windows with data after filtering for repeats and mapability.
#-----------------------------------------------------------------------------------------------------------------------

print(Sys.time())
options(scipen=100)
args = commandArgs(trailingOnly=TRUE)
chrom=args[1]
chrom=as.character(chrom)
print(chrom)
chrom<-unlist(strsplit(chrom,split="_"))
print(chrom)
'%ni%' <- Negate('%in%')
step=1

gors<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/gorilla.lst",sep="\t",as.is=T))

grps<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/groups.lst",sep="\t",as.is=T))


grn<-do.call(rbind,strsplit(do.call(rbind,strsplit(grps,split="/"))[,6],split="\\."))[,1]
inds<-list();finds<-list()
for ( i in (1:length(grps))) { inds[[i]]<-unlist(read.table(grps[i],sep="\t",as.is=T)); inds[[i]]<-inds[[i]][which(inds[[i]]%ni%c("Gorilla_beringei_graueri-Serufuli","Pan_troglodytes_schweinfurthii-Mgbadolite"))]   }  

indlist<-c(1:3,unlist(read.table("/scratch/devel/mkuhlwilm/ga/findivs.txt",header=F,as.is=T,sep="\t")))
# individuals for all species

for (i in 1:length(grps)) { finds[[i]]<-which(indlist%in%inds[[i]]) }
fugrp<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/species.lst",sep="\t",as.is=T))
fugrp<-cbind(c(paste("gorilla_",1:6,sep=""),paste("pan_",1:6,sep=""),paste("hum_",1:2,sep=""),paste("ora_",1:4,sep=""),paste("hupa",1:2,sep="")),c(rep(fugrp[1],6),rep(fugrp[2],6),rep(fugrp[3],2),rep(fugrp[4],4),rep(fugrp[5],2)))


minds<-list();funds<-list()
for ( i in (1:length(grps))) { minds[[i]]<-unlist(read.table(fugrp[i,2],sep="\t",as.is=T)); minds[[i]]<-minds[[i]][which(minds[[i]]%ni%c("Gorilla_beringei_graueri-Serufuli","Pan_troglodytes_schweinfurthii-Mgbadolite"))]   }
for (i in 1:length(grps)) { funds[[i]]<-which(indlist%in%minds[[i]] & indlist%ni%inds[[i]]) }
atyp=unique(fugrp[,2])


library("GenomicRanges")

#https://github.com/igvteam/igv/blob/master/genomes/sizes/hg19.chrom.sizes

require(data.table)

chilen<-data.table(
V1=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
V2= c(249250621,  243199373,  198022430,  191154276,  180915260,  171115067,  159138663,  146364022,  141213431,  135534747,  135006516,  133851895,  115169878,  107349540,  102531392,  90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566) )

#-----------------------------------------------------------------------------------------------------------------------

## for each species, get positions that have data, and that fulfill the filtering criteria

tfun<-function(intv, ag) { c(intv,sum(ranges(findOverlaps(IRanges(start=intv[1],end=intv[2]),ag),IRanges(start=intv[1],end=intv[2]),ag)@width)/40000) }
newchr<-c(1:22,"X")


yy=as.numeric(chrom[1])

typ=atyp[1]

spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
ft=2
chr=chrom
if (ft==2) {    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -S ",typ,"  /scratch/devel/mkuhlwilm/gvcfs/greatapeN_",chr,".vcf.gz -G |   intersectBed -wa -v -header -a stdin -b /scratch/project/devel/mkuhlwilm/RM.bed |   /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /scratch/project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep \"#|mapability=1\" | bedtools merge -i stdin",sep=""),intern=T) }
    print("data loaded")
    agt<-do.call(rbind,strsplit(agt,split="\t"))
    agtr<-IRanges(start=as.numeric(agt[,2]),end=as.numeric(agt[,3]))
    tevals<-chilen[which(chilen[,1]==chrom),2]; tevals<-cbind(seq(0,tevals[[1]],40000));  tevals<-cbind(tevals,tevals+39999)
    avals<-apply(tevals,1,tfun,ag=agtr)
    avals<-t(avals)
    aweight1<-as.data.frame(cbind(chr,avals[,c(1)]))
    aweight2<-as.character(as.double(avals[,3]))
    aweight3<-cbind(aweight1,format(aweight2, digits=3, scientific=F))
  
 # write out per chr
save(aweight3,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/",spec,"_weight",ft,".",chrom,".Robject",sep="")) 


q()
#------------------------------------------------------------------------------------------------------------------------

