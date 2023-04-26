#!/usr/bin/r

print(Sys.time())
options(scipen=100)
chrom=(commandArgs(TRUE))
chrom=as.character(chrom)
print(chrom)
chrom<-unlist(strsplit(chrom,split="_"))
print(chrom)
'%ni%' <- Negate('%in%')

## these steps are necessary for creating specific ingroup and outgroup patterns
grps<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/groups.lst",sep="\t",as.is=T))[1:6]
grn<-do.call(rbind,strsplit(do.call(rbind,strsplit(grps,split="/"))[,6],split="\\."))[,1]
inds<-list();finds<-list()
for ( i in (1:length(grps))) { inds[[i]]<-unlist(read.table(grps[i],sep="\t",as.is=T)) }
indlist<-c(1:3,unlist(read.table("/scratch/devel/mkuhlwilm/ga/findivsN.txt",header=F,as.is=T,sep="\t")))
for (i in 1:length(grps)) { finds[[i]]<-which(indlist%in%inds[[i]]) }
fugrp<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/species.lst",sep="\t",as.is=T))
fugrp<-cbind(c(paste("gorilla_",1:6,sep=""),paste("pan_",1:6,sep=""),paste("hum_",1:2,sep=""),paste("ora_",1:4,sep=""),paste("hupa",1:2,sep="")),c(rep(fugrp[1],6),rep(fugrp[2],6),rep(fugrp[3],2),rep(fugrp[4],4),rep(fugrp[5],2)))
minds<-list();funds<-list()
for (i in 1:length(grps)) { funds[[i]]<-which(indlist%in%minds[[i]] & indlist%ni%inds[[i]]) }
atyp=unique(fugrp[,2])

library("GenomicRanges")
chilen<-read.table("/home/devel/mkuhlwilm/hg19.chrom.sizes",sep="\t",header=F,nrow=24)[-c(21),]
chilen[,1]<-unlist(do.call(rbind,strsplit(as.character(chilen[,1]),split="hr"))[,2])
chilen<-chilen[,c(1:2)];cv<-seq(0,250000000,1000000)

# only do once
# ft=2
#  for (yy in (1:length(atyp))) {
#    typ=atyp[yy]
#    spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
#    alweig<-list()
#    save(alweig,file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
#    }
 
############################################################
## get positions that have data, and that fulfill the two filtering criteria
tfun<-function(intv, ag) { c(intv,sum(ranges(findOverlaps(IRanges(start=intv[1],end=intv[2]),ag),IRanges(start=intv[1],end=intv[2]),ag)@width)/1000) }
newchr<-c(1:22,"X")

yy=as.numeric(chrom[1])
typ=atyp[yy]
ft=as.numeric(chrom[2])
chr=as.character(chrom[3])
if (chr==23) { chr<-"X" }
spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
print(c(typ,ft,chr,spec))

    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -S ",typ," /scratch/devel/mkuhlwilm/gvcfs/greatapeN_",chr,".vcf.gz -G | intersectBed -wa -v -header -a stdin -b /project/devel/mkuhlwilm/RM.bed | /scratch/devel/mkuhlwilm/jvarkit/jvarkit/dist/vcfbigwig -T mapability -B /project/devel/mkuhlwilm/wgEncodeDukeMapabilityUniqueness35bp.bigWig /dev/stdin/ | egrep \"#|mapability=1\" | bedtools merge -i stdin",sep=""),intern=T)
    print("data loaded")
    agt<-do.call(rbind,strsplit(agt,split="\t"))
    agtr<-IRanges(start=as.numeric(agt[,2]),end=as.numeric(agt[,3]))
    tevals<-chilen[which(chilen[,1]==chr),2]; tevals<-cbind(seq(0,tevals,1000));  tevals<-cbind(tevals,tevals+999)
    avals<-apply(tevals,1,tfun,ag=agtr)
    avals<-t(avals)
    aweight1<-as.data.frame(cbind(chr,avals[,c(1)]))
    aweight2<-as.character(as.double(avals[,3]))
    aweight3<-cbind(aweight1,format(aweight2, digits=3, scientific=F))
    load(file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
    alweig[[chr]]<-aweight3
    save(alweig,file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
  print(Sys.time())
  q()

  

#####################################################
## merge and write into merged species-specific files
## run this part once to get the final weight file
library("GenomicRanges")
options(scipen=1, digits=8)

  ft=2
  for (yy in (1:length(atyp))) {
      typ=atyp[yy]; print(typ)
      spec=unlist(strsplit(unlist(strsplit(typ,split="/"))[6],split="\\."))[1]
      load(file=paste("/scratch/devel/mkuhlwilm/arch/",spec,"_weight",ft,".Robject",sep=""))  
      nweig<-list()
      for (i in c(1:22,"X")) { tw<-alweig[[which(names(alweig)==i)]]; nweig[[i]]<-tw[-c((nrow(tw)-1):(nrow(tw))),];colnames(nweig[[i]])<-c("chr","V2","format")      }
      aweight<-do.call(rbind,nweig)
      aweight1<-as.data.frame(matrix(unlist(aweight[,1:2]),ncol=2))
      aweight3<-cbind(aweight1,format(round(as.numeric(as.character(aweight[,3])), digits=3),nsmall=3, scientific=F))
      write.table(aweight3,paste("/scratch/devel/mkuhlwilm/arch/N",ft,"_",spec,"_weights_float.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)  
    }

