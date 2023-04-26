#!/usr/bin/r
print(Sys.time())
options(scipen=100)
ty=(commandArgs(TRUE))
ty=as.character(ty)
print(ty)


############################################################
## do the filtering for each chromosome
fifu<-function(input) { if(is.na(input)==FALSE) { op<-as.numeric(unlist(strsplit(input,split=",")));min(c(op[1]/(sum(op[1:2])),op[2]/sum(op[1:2]))) } else { NA } }
for (chrom in c(1:22,"X"))  {
  yy=0;yle=seq(0,250000000,1000000);print(chrom)
  repeat {
    yy=yy+1;print(yy)
    if (yy==length(yle)) { break}
    intvl<-c(yle[yy],yle[yy+1])
    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools query -r chr",chrom,":",intvl[1],"-",intvl[2]," -u -f '%POS;[%AD ];[%GT ]\n' /scratch/devel/mkuhlwilm/gvcfs/greatapeN_"$chrom".vcf.gz",sep=""),intern=T)
    if (length(agt)==0) { print("nolen"); next }
    agt<-do.call(rbind,strsplit(agt,split=";"))
    agp<-agt[,1]
    agg<-do.call(rbind,strsplit(agt[,3],split=" "))
    agt<-do.call(rbind,strsplit(agt[,2],split=" "))
    agt[which(agt==".")]<-NA
    agt[which(agt=="0")]<-NA
    rose<-rowSums(is.na(agt))
    agt<-agt[which(rose<(ncol(agt))),,drop=F]
    agg<-agg[which(rose<(ncol(agt))),,drop=F]
    agp<-agp[which(rose<(ncol(agt)))]
    agf<-matrix(as.numeric(unlist(apply(agt,c(1:2),fifu))),nrow=nrow(agt))
    agf[which(agf>0 & agf <=0.1 & agg=="0/1")]<--1
    if (length(agf)==0) { print("noagf"); next }
    agf<-rowSums(agf=="-1")
    if (length(which(agf>0))==0) { print("nopo"); next }
    aggp<-cbind(paste("chr",chrom,sep=""),as.numeric(agp[which(agf>0)])-1,agp[which(agf>0)])
    write.table(aggp,paste("/scratch/devel/mkuhlwilm/arch/filter/Na_",ty,"_",chrom,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F,append=T)
  } 
}


print(Sys.time())
q()
