#!/usr/bin/r

library(ggplot2)
require(gridExtra) 

options(scipen=100)
'%ni%' <- Negate('%in%')
mura=(1.235e-08)
sdfun<-function(inputtab) {
  vec <- unlist(inputtab, use.names = FALSE);   DIM <- dim(inputtab[[1]]);n <- length(inputtab)
  list.sd <- tapply(vec, rep(1:prod(DIM),times = n), sd);  attr(list.sd, "dim") <- DIM
  as.data.frame(list.sd)
}
tfun<-function(tb) { p=tb[3,2];q=tb[4,1];Tadm<-(p+(2*q))/(2*1000*mura);c(Tadm*19,((p)/((Tadm*2*1000*mura)))) }
ft=2
#tfun<-function(tb) {  ti<-(tb[3,2]+(2*tb[4,1]))/(2*1000*1.2e-08);c(ti*25,(tb[3,2]/(ti*2*1000*1.2e-08))) }

## preparation: define groups
grps<-unlist(read.table("/scratch/devel/mkuhlwilm/arch/groups.lst",sep="\t",as.is=T))[1:6]
grn<-do.call(rbind,strsplit(do.call(rbind,strsplit(grps,split="/"))[,6],split="\\."))[,1]
inds<-list();finds<-list()
for ( i in (1:length(grps))) { inds[[i]]<-unlist(read.table(grps[i],sep="\t",as.is=T)) }
indlist<-c(1:3,unlist(read.table("/scratch/devel/mkuhlwilm/ga/findivsN.txt",header=F,as.is=T,sep="\t")))
for (i in 1:length(grps)) { finds[[i]]<-which(indlist%in%inds[[i]]) }
funds<-list()
for (i in 1:length(grps)) { funds[[i]]<-which(indlist%ni%inds[[i]]) }


#######################################################################################################################
# read in the tables
tabls<-list();itab<-list()
typ=0
for (ty in (1:length(grn))) {
  for (tty in (1:length(grn[[ty]]))) {
    typ=typ+1;tabls[[typ]]<-list();itab[[typ]]<-list()
      for (ind in (1:length(sblists[[ty]][[tty]]))) { 
      t1<-try(as.character(unlist(read.table(paste("/scratch/devel/mkuhlwilm/arch/skov/model1/",grn1[ty],"/","trained_",sblists[[ty]][[tty]][ind],".hmm",sep=""),sep="\t",header=F))))
      if (inherits(t1,"try-error")==T) { next }
      tabls[[typ]][[ind]]<-matrix(as.numeric(rbind(unlist(strsplit(unlist(strsplit(t1[2],split=" = "))[2],split="\\[|\\]|, ")),unlist(strsplit(unlist(strsplit(t1[4],split=" = "))[2],split="\\[|\\]|, ")),t(matrix(unlist(strsplit(unlist(strsplit(t1[3],split=" = "))[2],split="\\[|\\]|, |,"))[c(2:4,6:8)],ncol=2)))[,c(2:3)]),ncol=2)
      itab[[typ]][[ind]]<-t1
      }
    }
  }

colfun<-function(input) { c(tfun(input),round(input[2,2]/(1000*mura)*19),round(input[2,1]/(1000*mura)*19)) }
reman<-list();calval<-list();calval2<-list();crval<-list();crval2<-list();coval<-list()
for (ty in (1:length(tabls))) {
  t2<-Reduce("+",tabls[[ty]])/length(tabls[[ty]]);rownames(t2)<-c("Starting","Emission","Transition_within","Transition_arch")
  reman[[ty]]<-t2
  sdv<-sdfun(tabls[[ty]])*2.58
  crval[[ty]]<-do.call(rbind,lapply(tabls[[ty]],tfun));crval2[[ty]]<-c(colMeans(crval[[ty]])*c(1,100),sd(crval[[ty]][,1])*2.58,sd(crval[[ty]][,2])*2.58*100)
  calval[[ty]]<-c(round(t2[1,2]*100,1),round(t2[2,2]/(1000*mura)*19),round(t2[2,1]/(1000*mura)*19))
  calval2[[ty]]<-c(round((t2[1,2]-sdv[1,2])*100,1),round((t2[1,2]+sdv[1,2])*100,1),round((t2[2,2]-sdv[2,2])/(1000*mura)*19),round((t2[2,2]+sdv[2,2])/(1000*mura)*19),round((t2[2,1]-sdv[2,1])/(1000*mura)*19),round((t2[2,1]+sdv[2,1])/(1000*mura)*19))
  coval[[ty]]<-do.call(rbind,lapply(tabls[[ty]],colfun))
}

calval<-do.call(rbind,calval);calval2<-do.call(rbind,calval2);crval2<-do.call(rbind,crval2)
rownames(calval2)<-gn
calval<-cbind(paste(round(crval2[,1])," (",round(crval2[,1]-crval2[,3]),"-",round(crval2[,1]+crval2[,3]),")",sep=""),paste(round(crval2[,2],4)," (",round(crval2[,2]-crval2[,4],4),"-",round(crval2[,2]+crval2[,4],4),")",sep=""), paste(calval[,2]," (",calval2[,3],"-",calval2[,4],")",sep=""),paste(calval[,3]," (",calval2[,5],"-",calval2[,6],")",sep=""))
colnames(calval)<-c("Admixture time","Archaic percentage","Archaic divergence","Within divergence");rownames(calval)<-gn
write.table(calval, file=paste("/scratch/devel/mkuhlwilm/arch/skov/fin_summary_model1flow_gor",".tsv",sep=""),sep="\t",col.names=T,row.names=T,quote=F)
crval<-coval
for (i in (1:length(crval))) { crval[[i]]<-cbind(gn[i],crval[[i]])}

crrval<-do.call(rbind,crval)
colnames(crrval)<-c("typ","T_admix","P_admix","Tcoal_extern","Tcoal_intern")

dival<-list()
for (i in (1:length(coval))) { tr<-coval[[i]][,3]/coval[[i]][,4];  dival[[i]]<-round(c(mean(tr),median(tr),sd(tr)),2) }
dival<-cbind(gn,do.call(rbind,dival))

#######################################################################################################################
## decoding output
options("scipen"=100)
chromlen<-read.table("/home/devel/mkuhlwilm/hg19.chrom.sizes",sep="\t",header=F)
chromlen<-chromlen[c(1:20,22:24),];gval=sum(as.numeric(as.character(chromlen[,2])))
ft=2
gnn<-gn

gtim=19
rbm=9.40e-09
mod=1;modl<-"model1";rval=0.90;rval2=0.95;typ=0
B=T

archfract<-list()
introtim<-list()
archfract2<-list()
introtim2<-list()
for (ty in (1:length(grn))) {
  print(c(gnn[ty],typ))
  for (tty in (1:length(grn[[ty]]))) {
    typ=typ+1;archfract[[typ]]<-list();introtim[[typ]]<-list();archfract2[[typ]]<-list();introtim2[[typ]]<-list()
    for (ind in (1:length(sblists[[ty]][[tty]]))) { 
      t1<-read.table(paste("/scratch/devel/mkuhlwilm/arch/skov/model1/",grn1[ty],"/","decoded_",sblists[[ty]][[tty]][ind],".Summary.txt.gz",sep=""),sep="\t",header=T)
      ## less stringent cutoff (0.9)
      arcfr<-t1[which(t1[,6]=="external"&as.numeric(as.character(t1[,8]))>rval&as.numeric(as.character(t1[,7]))>=5),c(3,5,7,8,2)]   
      intime<-1/(((1-(sum(arcfr[,2])/gval)))*(rbm)*(median(arcfr[,2])))
      introtim[[typ]][[ind]]<-intime
      archfract[[typ]][[ind]]<-arcfr

    ## more stringent cutoff (0.95)
      arcfr2<-t1[which(t1[,6]=="external"&as.numeric(as.character(t1[,8]))>rval2&as.numeric(as.character(t1[,7]))>=5),c(3,5,7,8,2)]   
      intime2<-1/(((1-(sum(arcfr2[,2])/gval)))*(rbm)*(median(arcfr[,2])))
      introtim2[[typ]][[ind]]<-intime2
      archfract2[[typ]][[ind]]<-arcfr2
      }
    }
  }
save(archfract,introtim,archfract2,introtim2,file=paste("/scratch/devel/mkuhlwilm/arch/skov/",modl,"/","data_decoded_gor.Robject",sep=""))  

## get stats
mod=1;modl<-"model1"
## high or low threshold
res<-1

load(file=paste("/scratch/devel/mkuhlwilm/arch/skov/",modl,"/","data_decoded_gor.Robject",sep=""))
if(res==1) {introtim<-introtim2;archfract<-archfract2 }
tims=list();frc<-list();frper=list();arvg<-c();nums=list();typ=0
for (ty in (1:length(grn))) {
  for (tty in (1:length(grn[[ty]]))) {
    typ=typ+1
    if (length(unlist(introtim[[typ]]))==0) { next }
    tims[[typ]]<-as.numeric(unlist(introtim[[typ]]))
    frc[[typ]]<-list();num=list()
    for (i in (1:length(archfract[[typ]]))) { frc[[typ]][[i]]<-sum(archfract[[typ]][[i]][,2]);num[[i]]<-archfract[[typ]][[i]][,1:2] }
    arvg[typ]<-round(mean(do.call(rbind,archfract[[typ]])[,2]))
    frc[[typ]]<-unlist(frc[[typ]]);frc[[typ]]<-frc[[typ]][which(frc[[typ]]!=0)]
    frper[[typ]]<-frc[[typ]]/gval*100;nums[[typ]]<-unique(do.call(rbind,num))
    }
  }
  
names(frper)<-unlist(grn)

tims2<-list()
for (i in (1:length(tims))) { tims2[[i]]<-tims[[i]]*gtim}

sumval<-list()
for (i in (1:length(frper))) {
  sumval[[i]]<-c(paste(round(mean(frper[[i]]),2)," (",round(min(frper[[i]]),2),"-",round(max(frper[[i]]),2),")",sep=""),paste(round(mean(frc[[i]]))," (",min(frc[[i]]),"-",max(frc[[i]]),")",sep=""),paste(round(mean(tims[[i]])*gtim)," (",round(min(tims[[i]])*gtim),"-",round(max(tims[[i]])*gtim),")",sep=""))
}

sumval<-do.call(rbind,sumval)
colnames(sumval)<-c("Proportion (%)","Fraction (Mbp)","Time (years)")
rownames(sumval)<-gnn
save(sumval,tims2,frper,frc,file=paste("/scratch/devel/mkuhlwilm/arch/skov/sumval_gr",ifelse(res==1,"","L"),"2",sep=""))

# write bed files
load(file=paste("/scratch/devel/mkuhlwilm/arch/skov/",modl,"/","data_decoded_gor.Robject",sep=""))
a1<-archfract[[9]]
a11<-list()
for (j in (1:length(a1))) { a11[[j]]<-cbind(paste("chr",as.character(a1[[j]][,5]),sep=""),as.numeric(as.character(a1[[j]][,1])),as.numeric(as.character(a1[[j]][,1]))+as.numeric(as.character(a1[[j]][,2])),inds[[6]][j])  }
a11<-do.call(rbind,a11)
a2<-archfract2[[9]]
a21<-list()
for (j in (1:length(a2))) { a21[[j]]<-cbind(paste("chr",as.character(a2[[j]][,5]),sep=""),as.numeric(as.character(a2[[j]][,1])),as.numeric(as.character(a2[[j]][,1]))+as.numeric(as.character(a2[[j]][,2])),inds[[6]][j])  }
a21<-do.call(rbind,a21)
#a2<-cbind(paste("chr",as.character(a2[,5]),sep=""),as.numeric(as.character(a2[,1])),as.numeric(as.character(a2[,1]))+as.numeric(as.character(a2[,2])))

write.table(a11,file="/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_all.bed",sep="\t",col.names=F,row.names=F,quote=F)
write.table(a21,file="/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed",sep="\t",col.names=F,row.names=F,quote=F)
