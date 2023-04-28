# Thu 26 Aug 2021 14:28:52 CEST
# recalculate heterozygosity for the empirical windows
#-----------------------------------------------------------------------------------------------------------------------

# module load gcc/6.3.0 R/3.4.2 xz/5.2.2  BCFTOOLS/1.6
#-----------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
require(data.table)
options(stringsAsFactors=F)
library(ape)
library(pegas)
chr=args[1]

# read in alleles directly from bcftools : %TGT = Translated genotype (e.g. C/A)
emp.data<-system(paste("bcftools query -f '%POS [%GT %TGT ]\n' /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chr,".vcf.gz",sep=""),intern=T)
emp.data<-do.call(rbind,strsplit(emp.data,split=" "))
pos<-emp.data[,1] 
pos<-as.numeric(pos)
window<-seq(min(pos), max(pos), 40000) # splits chromosome into 40kb non-overlapping regions
emp.data1<-cbind(pos, emp.data[,2:ncol(emp.data)])

# split into 2 dfs one where GTs as 0|1, the other with tgts as C|A
even_indexes<-seq(2,ncol(emp.data1),2)
odd_indexes<-seq(1,ncol(emp.data1),2)
emp.data.gts<-emp.data1[,c(1,even_indexes)]
emp.data.tgts<-emp.data1[,odd_indexes]

window<-seq(min(pos), max(pos), 40000) # splits chromosome into 40kb non-overlapping regions

config<-c(24,18,2,54)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
sum(config)
len=40000

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# function to calculate summary statistics per window:

het.in.win_function<-function(i) {
intvl<-c(window[[i]],window[[i+1]]) # interval to extract GTs from (40kb region)
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 1)
# count heterozygous GTs per individual & divide by window length

het_fun<-function(id){ 
  (length(which(win.intvl[,id]=="0|1")) + 
  length(which(win.intvl[,id]=="1|0")) +  
  length(which(win.intvl[,id]=="1|2")) +
  length(which(win.intvl[,id]=="2|1")) +
  length(which(win.intvl[,id]=="0|2")) +
  length(which(win.intvl[,id]=="2|0")) +
  length(which(win.intvl[,id]=="2|3")) +  
  length(which(win.intvl[,id]=="3|2"))) }

out.winhet<-lapply(2:ncol(win.intvl), het_fun) 

hetperid<-list(
unlist(out.winhet)[23:49],
unlist(out.winhet)[22],
unlist(out.winhet)[13:21],
unlist(out.winhet)[1:12]
  )

#-----------------------------------------------------------------------------------------------------------------------

# summary statistic 2) segregating sites

# count if site is segregating per individual for locus x -> iterate over loci 

splitgts_fun<-function(id){
do.call('rbind', strsplit(as.character(win.intvl[,id]),'|',fixed=TRUE)) 
}

win.intvl.hapl<-lapply(2:ncol(win.intvl), splitgts_fun)

y<-data.frame(matrix(unlist(win.intvl.hapl), ncol= sum(config)))
alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))

y[,] <- sapply(y[,],as.numeric)


# MK: You may only count non-fixed sites within each population, remove fixed sites in all gorillas

y<-y[which(rowSums(y)<sum(config)),]


alfrfun<-function(nput) {
   inds<-which(alinds==nput)
   op<-rowSums(y[,inds])/(length(inds))
   return(op)
   }
alfre<-lapply(1:4,alfrfun)

# where 1=EM, 2=EL, 3=WC, 4=WL

# if > 0 : seg site
alfre[[1]][alfre[[1]]>0] <- 1
alfre[[2]][alfre[[2]]>0] <- 1
alfre[[3]][alfre[[3]]>0] <- 1
alfre[[4]][alfre[[4]]>0] <- 1

lapply(alfre, sum)

out.seg.perwind<-unlist(lapply(alfre, sum)) # sum of segregating sites per popn for this window # note in the order of  1=EM, 2=EL, 3=WC, 4=WL (reverse to rest of stats - will need to rearrange at end)
#-----------------------------------------------------------------------------------------------------------------------

identifier<-cbind(chr,i, intvl[1],intvl[2]) # output identifier for this window - chr, window number, start pos, end pos

 # ie output 1) window identifier
                # 2) count of het sites per id for this window
                # 3) count of seg sites per id for this window
 return(list(identifier, hetperid, out.seg.perwind))
}
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

# extract identifier & intervals of non-empty windows
hold.windows=list()
for (j in 1:(length(window)-1)){
intvl<-c(window[[j]],window[[j+1]])
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]
if (length(win.intvl)!=0)
hold.windows[[j]]<-cbind(j, intvl[1], intvl[2], nrow(win.intvl)/len)
}


windows.wdata<-hold.windows[-which(sapply(hold.windows, is.null))]

fraction.wind=list()
for (j in 1:length(windows.wdata)){
fraction.wind[[j]]<-cbind(windows.wdata[[j]][1],windows.wdata[[j]][4])
}

df.fraction.wind <- data.frame(matrix(unlist(fraction.wind), nrow=length(fraction.wind), byrow=TRUE),stringsAsFactors=FALSE)

df1.fraction.wind<-na.omit(df.fraction.wind)

# apply het.in.win_function to windows with >1snp
	# & writing out the chr, window number & positions as well
stats.out=list()
for (j in 1:nrow(df1.fraction.wind)){
stats.out[[j]]<-het.in.win_function(df1.fraction.wind$X1[[j]])
}

save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/het_chr",chr,sep=""))

