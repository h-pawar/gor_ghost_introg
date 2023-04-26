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
fugrp<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/species.lst",sep="\t",as.is=T))
fugrp<-cbind(c(paste("gorilla_",1:6,sep=""),paste("pan_",1:6,sep=""),paste("hum_",1:2,sep=""),paste("ora_",1:4,sep=""),paste("hupa",1:2,sep="")),c(rep(fugrp[1],6),rep(fugrp[2],6),rep(fugrp[3],2),rep(fugrp[4],4),rep(fugrp[5],2)))
minds<-list();funds<-list()
for (i in 1:length(grps)) { funds[[i]]<-which(indlist%in%minds[[i]] & indlist%ni%inds[[i]]) }


############################################################
## first, do the filtering for each individual
    ft=2
    ## here, I load the genotypes, coverage, mapping quality and MQ0 for all individuals
    ## we walk through each chromosome in steps of 1M for quick processing
    yy=0;filtfin=list();yle=seq(0,250000000,1000000)
    repeat {
      yy=yy+1
      if (yy==length(yle)) { break}
      intvl<-c(yle[yy],yle[yy+1])
      agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools query -r chr",chrom,":",intvl[1],"-",intvl[2]," -u -f '%POS %REF;%ALT [%GT,] [%DP,] [%MQ,] [%MQ0,]\n' /scratch/devel/mkuhlwilm/arch/private/segsite_filt",ft,"_",chrom,".vcf.gz",sep=""),intern=T)
      if (length(agt)==0) { next }
      print(intvl[1])  
      agt<-do.call(rbind,strsplit(agt,split=" "))
      al<-do.call(rbind,strsplit(agt[,2],split=";"))
      gt<-do.call(rbind,strsplit(agt[,3],split=","))
      dp<-matrix(as.numeric(do.call(rbind,strsplit(agt[,4],split=","))),nrow=nrow(agt))
      mq0<-matrix(as.numeric(do.call(rbind,strsplit(agt[,6],split=","))),nrow=nrow(agt))
      mq<-matrix(as.numeric(do.call(rbind,strsplit(agt[,5],split=","))),nrow=nrow(agt))
    
      # which coverage is >5 and <100, which mapping quality is >20, which mq0/dp is < 10%, which sites are biallelic
      gt[which(dp<6)]<-NA
      gt[is.na(dp)]<-NA
      gt[which(dp>100)]<-NA
      gt[which(mq<20)]<-NA
      gt[is.na(mq)]<-NA
      gt[which((mq0/dp)>0.1)]<-NA
      gt[which(gt%ni%c("0/0","0/1","1/1"))]<-NA
      gt[which(al[,1]%ni%c("A","C","G","T")),]<-NA
      gt[which(al[,2]%ni%c("A","C","G","T")),]<-NA
      gt[which(gt=="0/0")]<-0
      gt[which(gt=="0/1")]<-1
      gt[which(gt=="1/1")]<-2
      nacon<-rowSums(is.na(gt))
    
      ## save biallelic information  
      allagt<-cbind(agt[,1],al,gt)[which(nacon!=ncol(gt)),]
    
      save(allagt,file=paste("/scratch/devel/mkuhlwilm/arch/private/sep_segsite/",chrom,"/segsite_filt",ft,"_",chrom,"_",yy,".Robject",sep=""))
      }
  print(Sys.time())
  q()


  
