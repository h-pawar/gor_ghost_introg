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

  spec=unlist(strsplit(grn[ip],split="_"))[1]
  spec<-ifelse(spec=="pan","pans",spec)
  spec<-ifelse(spec=="hum","human",spec)
  spec<-ifelse(spec=="ora","orang",spec)
  spec<-ifelse(spec%in%c("hupa1","hupa2"),"hupa",spec)
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
}

  

############################################################
#### now also as eigenstrat format

ft=2
for (chrom in c(1:22)) {
  print(chrom)
  yy=0;yle=seq(0,250000000,1000000);snp<-list();gtss<-list()
  repeat {
    print (yy)
    yy=yy+1;    intvl<-c(yle[yy],yle[yy+1])
    if (yy==length(yle)) { break}
    tagt<-try(load(file=paste("/scratch/devel/mkuhlwilm/arch/private/sep_segsite/",chrom,"/segsite_filt",ft,"_",chrom,"_",yy,".Robject",sep="")),silent=T)
    if (class(tagt) != "try-error") {  
      allagt<-allagt[,c(1:52,167),drop=F]
      nal<-rowSums(!is.na(allagt[,-c(1:3)]))
      rs<-rowSums(matrix(as.numeric(allagt[,-c(1:3)]),ncol=ncol(allagt)-3))
      apos<-allagt[which(nal>25 & rs!=0 & rs!=100),1:3,drop=F]
      gts<-allagt[which(nal>25& rs!=0 & rs!=100),-c(1:3),drop=F]
      gts[is.na(gts)]<-9
      gtss[[yy]]<-apply(format(gts),1,paste,collapse="")
      snp[[yy]]<-cbind(paste(chrom,apos[,1],sep=":"),chrom,"0.0",apos[,1:3,drop=F])
      }
    }
  snp<-do.call(rbind,snp)
  gtss<-unlist(gtss)
  write.table(snp,file=paste("/scratch/devel/mkuhlwilm/arch/eigen/",chrom,".snp",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  write.table(gtss,file=paste("/scratch/devel/mkuhlwilm/arch/eigen/",chrom,".geno",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
  }

idl<-do.call(rbind,strsplit(indlist[c(4:52,167)],split="-"))[,1]
idl<-gsub("Gorilla_beringei_beringei","MG",idl)
idl<-gsub("Gorilla_beringei_graueri","ELG",idl)
idl<-gsub("Gorilla_gorilla_gorilla","WLG",idl)
idl<-gsub("Gorilla_gorilla_diehli","CRG",idl)
idl<-gsub("pygmaeus_SRS396836","Orang",idl)

tl<-cbind(indlist[c(4:52,167)],"U",idl)
write.table(tl,file="/scratch/devel/mkuhlwilm/arch/eigen/gori.ind",sep="\t",col.names=F,row.names=F,quote=F)

pops<-unique(tl[,3])
fosta<-expand.grid(pops,pops,pops,pops)
fosta<-fosta[which(fosta[,1]!=fosta[,2]&fosta[,2]!=fosta[,3]&fosta[,1]!=fosta[,3]&fosta[,1]!=fosta[,4]&fosta[,2]!=fosta[,4]&fosta[,3]!=fosta[,4]),]
write.table(fosta,file="/scratch/devel/mkuhlwilm/arch/eigen/tosta.txt",sep="\t",row.names=F,col.names=F,quote=F)
results<-list()
save(results,file=paste("/scratch/devel/mkuhlwilm/arch/eigen/f4a",sep=""))

for chrom in {1..22}; do echo $chrom
cat /scratch/devel/mkuhlwilm/arch/eigen/$chrom.geno >> /scratch/devel/mkuhlwilm/arch/eigen/gori.geno
cat /scratch/devel/mkuhlwilm/arch/eigen/$chrom.snp >> /scratch/devel/mkuhlwilm/arch/eigen/gori.snp
done
