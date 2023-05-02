#module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
# plot nj trees for private regions (introg regions unique to 1 id of pop)

library(phangorn) 
library(ape) 
library('pegas')

#-----------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
id=args[1] 
SPE=args[2]
nput=args[3]

#-----------------------------------------------------------------------------------------------------------------------

# empirical

load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/nj/",nput,".id.",id,".rootedNJploutg.Robj.tree",sep=""),verbose=T)
load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/nj/",nput,".id.",id,".Robj.tree",sep=""),verbose=T)
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/nj/plots/",nput,".id.",id,".private.nj.pdf",sep="")) 
plotBS(rootedNJ,out[[2]],type="phylogram", cex=0.7)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------

# random
load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/nj/",nput,".id.",id,".private.random.rootedNJploutg.Robj.tree",sep=""),verbose=T)
randomNJ<-rootedNJ
load(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/nj/",nput,".id.",id,".private.random.Robj.tree",sep=""),verbose=T)
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/nj/plots/",nput,".id.",id,".private.random.nj.pdf",sep="")) 
plotBS(randomNJ,random_out[[2]],type="phylogram", cex=0.7)
dev.off()



#-----------------------------------------------------------------------------------------------------------------------
