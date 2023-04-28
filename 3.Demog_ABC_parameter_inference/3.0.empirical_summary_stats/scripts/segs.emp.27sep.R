# Mon 27 Sep 2021 10:16:39 CEST
# new statistics to calculate in the empirical data which were introduced in the new simulation batch:
    #full_segsites is the number of population-wise fixed sites and the number of population-wise segregating sites; 
    # fix_op is the fixed sites per individual (same format as heterozygotes in your script)

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

window<-seq(min(pos), max(pos), 40000) 
config<-c(24,18,2,54)
cofig<-cbind(c(0,cumsum(config/2)[-length(config)]),cumsum(config/2))
sum(config)
len=40000


#------------------------------------------------------------------------------------------------------------------------

etest<-gsub("2", "1", emp.data.gts[,c(2:ncol(emp.data.gts))], fixed=TRUE)
etest<-gsub("3", "1", etest[,c(1:ncol(etest))], fixed=TRUE)

emp.data.gts<-cbind(emp.data.gts[,1],etest)

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
# function to calculate summary statistics per window:

het.in.win_function<-function(i) {

intvl<-c(window[[i]],window[[i+1]]) # interval to extract GTs from (40kb region)
win.intvl<-emp.data.gts[intvl[1] <= pos & pos < intvl[2] ,]

cofig2<-cbind(c(0,cumsum(config)[-length(config)]),cumsum(config))

#-----------------------------------------------------------------------------------------------------------------------

# 1) #population-wise fixed sites and the number of population-wise segregating sites

splitgts_fun<-function(id){
do.call('rbind', strsplit(as.character(win.intvl[,id]),'|',fixed=TRUE)) 
}

win.intvl.hapl<-lapply(2:ncol(win.intvl), splitgts_fun)

y<-data.frame(matrix(unlist(win.intvl.hapl), ncol= sum(config)))
alinds<-c(rep(1,config[1]),rep(2,config[2]),rep(3,config[3]),rep(4,config[4]))

y[,] <- sapply(y[,],as.numeric)
fxd<-which(rowSums(y)<sum(config))
y<-y[fxd,]

segfun<-function(y) {
      idx<-y
      pofun<-function(typ,idx) { pop=c((typ[1]+1):typ[2]);op<-table(rowSums(idx[,pop,drop=F]))
      op1=op[which(as.numeric(names(op))==length(pop))];op1=ifelse(length(op1)==0,0,op1)
      op2=sum(op[which(as.numeric(names(op))>0 & as.numeric(names(op))<length(pop))]);op2=ifelse(length(op2)==0,0,op2)
      c(op1,op2)}
      apply(cofig2,1,pofun,idx=idx)
  }

eseg<-segfun(y)[,c(4,3,2,1)] # in order of WL,WC,EL,EM

#-----------------------------------------------------------------------------------------------------------------------
fixedid_fun<-function(id){
    length(which(win.intvl[fxd,id]=="1|1")) }

out.fixedid<-lapply(2:ncol(win.intvl), fixedid_fun)

fixedsitesperid<-list(
unlist(out.fixedid)[23:49],
unlist(out.fixedid)[22],
unlist(out.fixedid)[13:21],
unlist(out.fixedid)[1:12]
  )


#-----------------------------------------------------------------------------------------------------------------------
identifier<-cbind(chr,i, intvl[1],intvl[2]) # output identifier for this window - chr, window number, start pos, end pos

 ## output
 # ie output 1) window identifier
                # 2) number of population-wise fixed sites and the number of population-wise segregating sites; 
                # 3) fixed sites per individual

return(list(identifier, eseg, fixedsitesperid))
}


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

#remove rows with nas - ie windows where there was only 1 snp
df1.fraction.wind<-na.omit(df.fraction.wind)

# apply het.in.win_function to windows with >1snp
  # & writing out the chr, window number & positions as well
stats.out=list()
for (j in 1:nrow(df1.fraction.wind)){
stats.out[[j]]<-het.in.win_function(df1.fraction.wind$X1[[j]])
}

save(stats.out,file=paste("/scratch/devel/hpawar/admix/abc/emp.data/test/segs/seg_chr",chr,sep=""))

