#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14

library(phangorn) 
library(ape) 
library('pegas')
#id=1

#-----------------------------------------------------------------------------------------------------------------------
# pass which scenario from the cd line
args = commandArgs(trailingOnly=TRUE)
id=args[1] 
SPE=args[2]
nput=args[3]

#-----------------------------------------------------------------------------------------------------------------------


# empirical

load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",nput,".id.",id,".rootedNJploutg.Robj.tree",sep=""),verbose=T)
load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",nput,".id.",id,".Robj.tree",sep=""),verbose=T)
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/",SPE,"/",nput,"/id.",id,".nj.pdf",sep="")) 
plotBS(rootedNJ,out[[2]],type="phylogram", cex=0.7)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------

# random
load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",nput,"/",nput,".id.",id,".random.rootedNJploutg.Robj.tree",sep=""),verbose=T)
randomNJ<-rootedNJ
load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",nput,"/",nput,".id.",id,".random.Robj.tree",sep=""),verbose=T)
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/",SPE,"/",nput,"/plots/id.",id,".random.nj.pdf",sep="")) 
plotBS(randomNJ,random_out[[2]],type="phylogram", cex=0.7)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------

#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/GBB/gbb/plots
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/trees/random/GBG/gbg/plots
