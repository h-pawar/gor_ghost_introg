#!/usr/bin/r

print(Sys.time())
options(scipen=100)
chrom=(commandArgs(TRUE))
chrom=as.character(chrom)
print(chrom)
'%ni%' <- Negate('%in%')

## preparation: define groups
grps<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/groups.lst",sep="\t",as.is=T))[1:6]
grn<-do.call(rbind,strsplit(do.call(rbind,strsplit(grps,split="/"))[,6],split="\\."))[,1]
inds<-list();finds<-list()
for ( i in (1:length(grps))) { inds[[i]]<-unlist(read.table(grps[i],sep="\t",as.is=T)) }
indlist<-c(1:3,unlist(read.table("/scratch/devel/mkuhlwilm/ga/findivsN.txt",header=F,as.is=T,sep="\t")))
for (i in 1:length(grps)) { finds[[i]]<-which(indlist%in%inds[[i]]) }
funds<-list()
for (i in 1:length(grps)) { funds[[i]]<-which(indlist%ni%inds[[i]]) }


############################################################
## correct and complete files
## here, instead of chromosomes, I use the groups
     
options("scipen"=100)
ip=as.numeric(chrom)

ft=2 

  spec="gorilla"
  tm1<-read.table(paste("/scratch/devel/mkuhlwilm/arch/N",ft,"_",spec,"_weights_float.txt.gz",sep=""),as.is=T,header=F)
  tm<-paste(tm1[,1],":",tm1[,2],sep="")
  for (ind in (funds[[ip]])) {
      print(paste(c(grn[ip],ind)))
      tr<-list()
      for (chr in c(1:22,"X")) { tr[[chr]]<-read.table(paste("/scratch/devel/mkuhlwilm/arch/private/",grn[ip],"/",ft,"/",ind,"_",chr,"_observations.txt",sep=""),as.is=T,fill=T,sep="\t") } 
      tr<-do.call(rbind,tr)
      tr[,1]<-gsub("chr","",tr[,1])
      tr<-cbind(tr,paste(tr[,1],tr[,2],sep=":"))
      tr<-merge(tm,tr,by.x=1,by.y=5,all.x=T,all.y=F)
      tr<-cbind(tr[,-1],do.call(rbind,strsplit(as.character(tr[,1]),split=":")))
      tr<-tr[,c(5:6,3:4)]
      tr[is.na(tr[,3]),3]<-0
      ## some are duplicates without data (from previous batch)
      tr<-tr[-which (tr[,3]>0& tr[,4]==""),]
      tr<-tr[order(as.numeric(as.character(tr[,2]))),]
      tr<-tr[order(as.numeric(as.character(tr[,1]))),]
      tr<-unique(tr)
      tr[which(tr[,3]>0&tm1[,3]==0),4]<-""
      tr[which(tr[,3]>0&tm1[,3]==0),3]<-0
      tr[is.na(tr[,4]),4]<-""
      write.table(tr, file=paste("/scratch/devel/mkuhlwilm/arch/private/",grn[ip],"/fin/",ind,"_observations",ft,".txt",sep=""),sep="\t",col.names=F,row.names=F,quote=F) }
    }

print(Sys.time())
q()

