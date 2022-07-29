#!/usr/bin/r
print(Sys.time())
options(scipen=100)
# commenting out the use of ty
#ty=(commandArgs(TRUE))
#ty=as.character(ty)
#print(ty)
# ty = spec

# 2: filter imbalanced hets
# follows from merging the "new" Tumani with chr9 of the "old sample set" - which used the script - /scratch/devel/hpawar/admix/sstar/scripts/gatk.combinevariants.arr

# Tue 28 Jun 2022 14:12:11 CEST
# changing paths to those of new tumani merged with the old sample set
# amending below from MK: filter_hets.R script
  # the rest of the original script - extracts coverage for chr & seg sites # shoudl run these subsequently *

############################################################
## do the filtering for each chromosome
#fifu<-function(input) { if(is.na(input)==FALSE) { op<-as.numeric(unlist(strsplit(input,split=",")));min(c(op[1]/op[2],op[2]/op[1])) } else { NA } }
fifu<-function(input) { if(is.na(input)==FALSE) { op<-as.numeric(unlist(strsplit(input,split=",")));min(c(op[1]/(sum(op[1:2])),op[2]/sum(op[1:2]))) } else { NA } }
#for (chrom in c(1:22,"X"))  {
 # try(system(paste("rm /scratch/devel/mkuhlwilm/arch/filter/str_",ty,"_",chrom,".bed",sep="")))
for (chrom in c(9))  {
  yy=0;yle=seq(0,250000000,1000000);print(chrom)
  repeat {
    yy=yy+1;print(yy)
    if (yy==length(yle)) { break}
    intvl<-c(yle[yy],yle[yy+1])
#    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools query -r chr",chrom,":",intvl[1],"-",intvl[2]," -u -f '%POS;[%AD ];[%GT ]\n' /scratch/devel/mkuhlwilm/arch/subsets/1_",ty,"_",chrom,".vcf.gz",sep=""),intern=T)
    agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools query -r chr",chrom,":",intvl[1],"-",intvl[2]," -u -f '%POS;[%AD ];[%GT ]\n' /scratch/devel/hpawar/admix/sstar/vcf_24jun22/greatape_chr9.vcf.gz",sep=""),intern=T)
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
#    write.table(aggp,paste("/scratch/devel/mkuhlwilm/arch/filter/Na_",ty,"_",chrom,".bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F,append=T)
    write.table(aggp,paste("/scratch/devel/hpawar/admix/sstar/vcf_24jun22/Na_gor_chr9.bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F,append=T)
  } 
}


print(Sys.time())
q()
