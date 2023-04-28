
#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.14
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(adegenet)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
#identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))


identifiers<-cbind(samples.in.vcf, c(rep("GBB-Bwindi",1),rep("GBB-Virunga",2),rep("GBB-Bwindi",2),rep("GBB-Virunga",1),rep("GBB-Bwindi",2),rep("GBB-Virunga",4),rep("GBG",5),rep("Mt Tshiaberimu",1), rep("GBG",3), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")


#> unique(identifiers[,2])
#[1] "GBB-Bwindi"     "GBB-Virunga"    "GBG"            "Mt Tshiaberimu"
#[5] "GGD"            "GGG"  

# set colours

#91C79B # EL
#3A7E99 # MG

#EF7C69 # CR
#F2BA5B # WL 

# or maybe coudl change shape for bwindi vs virunga 

#4595b5 # virunga
#2f677d # bwindi
#b8dbbe # mt T

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# read in .random.pca

# 1st random rep 
scen=(c(rep("GBG",9)))
# test
SPE="GBG"
spe="gbg"
random_pca_obj=list()
for (i in 1:length(scen)){
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/random/",SPE,"/",spe,".",i,".random.pca",sep=""),verbose=T)
random_pca_obj[[i]]<-tes
}


#  length(random_pca_obj[[1]][[1]])
#[1] 31

random_pca_obj[[1]][[1]][[31]]<-identifiers[,2]


#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


random_plot_pca<-function(input,ind,xax,yax,x1,y1) {

pc1pc2<-input[[ind]][[1]][,c(x1:y1,31)]
colnames(pc1pc2)<-c("pc1","pc2","pop")
xaxis<-input[[ind]][[xax]]
yaxis<-input[[ind]][[yax]]

p<-ggplot(data = pc1pc2, aes(x = pc1, y = pc2, color = pop)) +  geom_point(alpha = 0.8,size=3) + xlab(xaxis) + ylab(yaxis) + 
theme_classic()  +  theme(legend.title=element_blank()) + 
	scale_color_manual(values=c("#2f677d","#68abc5","#91C79B","#EF7C69","#F2BA5B","#a4b4a7"))
return(p)
}


i=1

rand_12<-random_plot_pca(random_pca_obj,i,3,4,1,2)
rand_34<-random_plot_pca(random_pca_obj,i,5,6,3,4)
hold_plots<-ggarrange(rand_12,rand_34, ncol = 2, nrow = 1,common.legend = TRUE,legend="bottom")

SPE="GBG"
spe="gbg"
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/",spe,".1randomreppca.10nov22.1.pdf",sep="")) 
hold_plots
dev.off()
