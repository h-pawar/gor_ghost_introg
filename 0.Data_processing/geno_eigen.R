#!/usr/bin/r

print(Sys.time())
options(scipen=100)
chrom=(commandArgs(TRUE))
chrom=as.character(chrom)
print(chrom)
'%ni%' <- Negate('%in%')

############################################################
#### now also as eigenstrat format
indlist<-c(1:3,unlist(read.table("/scratch/devel/mkuhlwilm/ga/findivsN.txt",header=F,as.is=T,sep="\t")))

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

q()

# in bash:
for chrom in {1..22}; do echo $chrom
cat /scratch/devel/mkuhlwilm/arch/eigen/$chrom.geno >> /scratch/devel/mkuhlwilm/arch/eigen/gori.geno
cat /scratch/devel/mkuhlwilm/arch/eigen/$chrom.snp >> /scratch/devel/mkuhlwilm/arch/eigen/gori.snp
done
