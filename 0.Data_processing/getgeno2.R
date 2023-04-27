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
## get the individuals into proper observations files

library("GenomicRanges")
chilen<-read.table("/home/devel/mkuhlwilm/hg19.chrom.sizes",sep="\t",header=F,nrow=24)[-c(21),]
chilen[,1]<-unlist(do.call(rbind,strsplit(as.character(chilen[,1]),split="hr"))[,2])
chilen<-chilen[,c(1:2)];cv<-seq(0,250000000,1000000)
tevals<-chilen[which(chilen[,1]==chrom),2]; tevals<-cbind(seq(0,tevals,1000));  tevals<-cbind(tevals,tevals+999)
tevals<-IRanges(start=tevals[,1],end=tevals[,2])
rn<-sample(1:10000)

tfun<-function(v,v1,v2) { paste(v2[which(v1==v)],collapse="," ) }

  ft=2  
  ##  first time
    ## this function turns the matrix into a vcftools freq type file, for a given individual, and immediately does the work of the Filtervariants script from Skov, just correctly this time; the output is directed to the final observations file
    tfunc<-function(ip,funds,finds,agt,apo,chrom,grn,intv) {
      indp=finds[[ip]]
      op=funds[[ip]]
      sbt<-matrix(as.numeric(agt[,indp]),ncol=length(indp))
      sbt<-cbind(rowSums(sbt,na.rm=T),(length(indp)*2)-rowSums(is.na(sbt))*2)
      ibt=matrix(as.numeric(agt[,op]),ncol=length(op))
      ## select sites that are fixed in outgroup, different in ingroup, have coverage in at least half of the individuals of both sets
        selec<-which((((sbt[,1]==0 & rowSums(ibt,na.rm=T)>0))|(sbt[,1]==sbt[,2] & rowSums(ibt,na.rm=T)<(rowSums(is.na(ibt))*2))) & sbt[,2]>=(length(indp)) & ((length(op)*2)-(rowSums(is.na(ibt))*2))>=(length(op)))
        if(length(selec)==0) { return(NULL) }
        ibt=ibt[selec,,drop=F]
        nop<-cbind(paste("chr",chrom,sep=""),apo[selec,1,drop=F],2,2)
        sbt<-sbt[selec,,drop=F]
      ## for each individual, select positions that are not NA, and apply the script
      nm<-sample(rn,1);write.table(cbind(nop[,1],as.numeric(nop[,2])-1,nop[,2]),paste("/dev/shm/mkdata/tmp",nm,".bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F) 
      ancal<-do.call(rbind,strsplit(system(paste("bedtools getfasta -tab -fi /scratch/devel/mkuhlwilm/pseudoarc/ga_ancestral.fa -bed /dev/shm/mkdata/tmp",nm,".bed -fo stdout",sep=""),intern=T),split="\t"));ancal[,1]<-do.call(rbind,strsplit(ancal[,1],split="-"))[,2]
      tv<-tevals[tevals@start>=intv[1]&tevals@start<intv[2]]
      for (ind in (1:ncol(ibt))) {
        tmp<-cbind(nop[!is.na(ibt[,ind]),2,drop=F],cbind(apo[selec,c(2,3),drop=F],2-ibt[,ind],ibt[,ind])[!is.na(ibt[,ind]),,drop=F],sbt[!is.na(ibt[,ind]),,drop=F],ancal[!is.na(ibt[,ind]),2,drop=F])
        sel<-ifelse(tmp[,2]==tmp[,8] & tmp[,6]==0 & tmp[,5]>0,1,ifelse(tmp[,3]==tmp[,8] & tmp[,6]==tmp[,7] & tmp[,4]>0,1,0))
        myout<-countOverlaps(tv,IRanges(start=as.numeric(tmp[which(sel==1),1]),end=as.numeric(tmp[which(sel==1),1])))
        if (length(myout)>0 &sum(sel)>0 )  {
          myout2<-as.data.frame(findOverlaps(IRanges(start=as.numeric(tmp[which(sel==1),1]),end=as.numeric(tmp[which(sel==1),1])),tv))[,2,drop=F]
          myout2<-unlist(lapply(1:1000,tfun,v1=myout2,v2=tmp[which(sel==1),1,drop=F]))
          write.table(cbind(paste("chr",chrom,sep=""),tv@start,myout,myout2), paste("/scratch/devel/mkuhlwilm/arch/private/",grn[ip],"/",ft,"/",op[ind],"_",chrom,"_observations.txt",sep=""),sep="\t",col.names=F,row.names=F,quote=F,append=T) }
        }
      print(grn[ip])
      system(paste("rm /dev/shm/mkdata/tmp",nm,".bed",sep=""),intern=F)
  
      
    ## run it over the chromosome
    print(chrom)
    yy=0;yle=seq(0,250000000,1000000)
    repeat {
      print (yy)
      yy=yy+1;    intvl<-c(yle[yy],yle[yy+1])
      if (yy==length(yle)) { break}
      tagt<-try(load(file=paste("/scratch/devel/mkuhlwilm/arch/private/sep_segsite/",chrom,"/segsite_filt",ft,"_",chrom,"_",yy,".Robject",sep="")),silent=T)
      if (class(tagt) != "try-error") {  
        apos<-allagt[,1:3]
        no<-lapply(1:length(finds),tfunc,finds=finds,funds=funds,apo=apos,agt=allagt,chrom=chrom,grn=grn,intv=intvl)
      }  
      }
  print(Sys.time())
  q()


  
