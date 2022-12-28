#Wed 27 Jul 2022 12:06:18 CEST

# plot pcas of private intersect regions (99_s* - strict skov)

#amending from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/plot_pca_sstarskov_27jul22.R # ran interactively


# amending from /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/plot_pca_skov_16jun22.R & referring back to /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/simul_postabc/plot_pca_skov_23may22.R


# same x axis, y axis scale for pcs of empirical & random reps, ie adding equal axes for comparable plots (emp & random pc for same individual)

#-----------------------------------------------------------------------------------------------------------------------
# MK Wed, 18 May, 21:01 
#for the PCA, one random replicate should be enough for the comparison,
# as it should show the "usual" genome-wide pattern. 
# The best would be side-by-side PC1/2 for introgressed vs random regions,
#  and PC3/4 for both as well on the same page.
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
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

# function to check x axis & y axis limits (to have the same scales for panels of the same comparison, eg pc1pc2 empirical vs random rep)
# where input = empirical_scen, input2 = random_scen
    # input=pca_obj, input2=random_pca_obj # for testing, x1=1, y1=2 (ie testing for PC1 PC2)
check_limits_pca<-function(input,x1,y1, input2) {
hold_min<-list()
hold_max<-list()

for (ind in (1:length(input))) {
hold_min[[ind]]<-cbind(min(input[[ind]][[1]][,c(x1:y1,31)][,1]),min(input[[ind]][[1]][,c(x1:y1,31)][,2]))
hold_max[[ind]]<-cbind(max(input[[ind]][[1]][,c(x1:y1,31)][,1]),max(input[[ind]][[1]][,c(x1:y1,31)][,2]))
}

hold_min2<-list()
hold_max2<-list()

for (ind in (1:length(input2))) {
hold_min2[[ind]]<-cbind(min(input2[[ind]][[1]][,c(x1:y1,31)][,1]),min(input2[[ind]][[1]][,c(x1:y1,31)][,2]))
hold_max2[[ind]]<-cbind(max(input2[[ind]][[1]][,c(x1:y1,31)][,1]),max(input2[[ind]][[1]][,c(x1:y1,31)][,2]))
}


# for each individual of the popn, extract min val for the empirical pc, min for the random rep for both axes (& the same for max)
min_x<-cbind(do.call(rbind, hold_min)[,1],do.call(rbind, hold_min2)[,1])
min_y<-cbind(do.call(rbind, hold_min)[,2],do.call(rbind, hold_min2)[,2])


max_x<-cbind(do.call(rbind, hold_max)[,1],do.call(rbind, hold_max2)[,1])
max_y<-cbind(do.call(rbind, hold_max)[,2],do.call(rbind, hold_max2)[,2])


xlim<-list()
for (i in (1:nrow(min_x))) {
xlim[[i]]<-cbind(min(min_x[i,]), max(max_x[i,]))
}

ylim<-list()
for (i in (1:nrow(min_y))) {
ylim[[i]]<-cbind(min(min_y[i,]), max(max_y[i,]))
}

# axis limits for each individual for PCs 1 (xlim), 2 (ylim)

x<- do.call(rbind, xlim)
y<- do.call(rbind, ylim)

# add margin for the plots, ie each of xlim -10
x[,1]<-x[,1]-10
x[,2]<-x[,2]+10

y[,1]<-y[,1]-10
y[,2]<-y[,2]+10

z<-cbind(x,y)

return(z)
}
#-----------------------------------------------------------------------------------------------------------------------

    # input=pca_obj, input2=random_pca_obj # for testing, x1=1, y1=2 (ie testing for PC1 PC2)
#check_limits_pca(pca_obj,1,2,random_pca_obj) # works

#> check_limits_pca(pca_obj,1,2,random_pca_obj)
#           [,1]     [,2]      [,3]     [,4]
# [1,] -278.8890 176.4082 -211.3981 173.1496
# [2,] -297.2078 191.7174 -236.6585 177.0709
# [3,] -301.2560 189.6148 -234.6338 169.0107
# [4,] -282.2050 179.0733 -217.4041 188.8875
# [5,] -295.8610 192.9288 -232.1909 196.9267
# [6,] -300.1648 193.7878 -231.4551 169.3754
# [7,] -287.5825 186.7556 -222.5503 186.0367
# [8,] -292.2520 192.1551 -224.6602 201.7710
# [9,] -290.5079 189.3965 -233.6349 166.3475
#[10,] -266.7869 166.8228 -203.0646 164.5097
#[11,] -279.1568 179.5897 -218.0020 159.3463
#[12,] -286.9119 186.4577 -225.4129 161.5617

# then add these into the plot_pca function, to call these as limits

# each row is 1 individual, : xlower, xupper, ylower, yupper
#> lims[1,]
#[1] -278.8890  176.4082 -211.3981  173.1496

#+ xlim(limits[i,][1], limits[i,][2]) + ylim(limits[i,][3], limits[i,][4])

# where limits: need to be calculated separately for pc1/2, pc3/4

#lims_12<-check_limits_pca(pca_obj,1,2,random_pca_obj)
#lims_34<-check_limits_pca(pca_obj,3,4,random_pca_obj)


#-----------------------------------------------------------------------------------------------------------------------
# functions to run for all individuals of the first scenario

plot_pca<-function(input,ind,xax,yax,x1,y1,limits) {

pc1pc2<-input[[ind]][[1]][,c(x1:y1,31)]
colnames(pc1pc2)<-c("pc1","pc2","pop")
xaxis<-input[[ind]][[xax]]
yaxis<-input[[ind]][[yax]]
#title<-paste("empirical ",x,"",sep="")
title<-paste("empirical ",ind,"",sep="")

li<-limits[ind,] # each row is 1 individual, : xlower, xupper, ylower, yupper

p<-ggplot(data = pc1pc2, aes(x = pc1, y = pc2, color = pop)) +   geom_point(alpha = 0.8) + xlab(xaxis) + ylab(yaxis) + ggtitle(title) + theme_classic() + xlim(limits[i,][1], limits[i,][2]) + ylim(limits[i,][3], limits[i,][4]) + scale_color_manual(values=c("#3A7E99", "#91C79B", "#EF7C69","#F2BA5B"))
return(p)
}


random_plot_pca<-function(input,ind,xax,yax,x1,y1,limits) {

pc1pc2<-input[[ind]][[1]][,c(x1:y1,31)]
colnames(pc1pc2)<-c("pc1","pc2","pop")
xaxis<-input[[ind]][[xax]]
yaxis<-input[[ind]][[yax]]
#title<-paste("random ",x,"",sep="")
title<-paste("random ",ind,"",sep="")
li<-limits[ind,] # each row is 1 individual, : xlower, xupper, ylower, yupper

p<-ggplot(data = pc1pc2, aes(x = pc1, y = pc2, color = pop)) +   geom_point(alpha = 0.8) + xlab(xaxis) + ylab(yaxis) + ggtitle(title) + theme_classic() + xlim(limits[i,][1], limits[i,][2]) + ylim(limits[i,][3], limits[i,][4]) + scale_color_manual(values=c("#3A7E99", "#91C79B", "#EF7C69","#F2BA5B"))
return(p)
}



#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# for GBB

# 1) read in empirical data
scen=(c(rep("GBB",12)))
# test
SPE="GBB"
spe="gbb"

pca_obj=list()
for (i in 1:length(scen)){
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/pca/",spe,".",i,".pca",sep=""),verbose=T)
pca_obj[[i]]<-tes
}

#-----------------------------------------------------------------------------------------------------------------------

# object tes contains:
#tes<-list(pc1pc2,k,xaxis,yaxis,x3axis,y4axis)

# x[[1]] = pc1pc2 - df of all principal components + last column is pop identifiers
# x[[2]] = k - % of variance explained by each component
# x[[3]] = xaxis for PC1
# x[[4]] = yaxis for PC2
# x[[5]] = xaxis for PC3
# x[[6]] = yaxis for PC4

#-----------------------------------------------------------------------------------------------------------------------

# 2) read in equivalent random data
scen=(c(rep("GBB",12)))
# test
SPE="GBB"
spe="gbb"

random_pca_obj=list()
for (i in 1:length(scen)){
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/pca/",spe,".",i,".random.pca",sep=""),verbose=T)
random_pca_obj[[i]]<-tes
}

# obtain plot limits (xlim, ylim)
lims_12<-check_limits_pca(pca_obj,1,2,random_pca_obj)
lims_34<-check_limits_pca(pca_obj,3,4,random_pca_obj)

hold_plots=list()
for (i in 1:length(scen)){
test_12<-plot_pca(pca_obj,i,3,4,1,2,lims_12)
test_34<-plot_pca(pca_obj,i,5,6,3,4,lims_34)
rand_12<-random_plot_pca(random_pca_obj,i,3,4,1,2,lims_12)
rand_34<-random_plot_pca(random_pca_obj,i,5,6,3,4,lims_34)
hold_plots[[i]]<-ggarrange(test_12,rand_12,test_34,rand_34, ncol = 2, nrow = 2)
}


#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBG/plot_pca
#[hpawar@login1 ~]$ mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/GBB/plot_pca

# replot

SPE="GBB"
spe="gbb"
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/plot_pca/",spe,".allids.privateintrogregions.pca.30sep22.pdf",sep="")) 
hold_plots
dev.off()


# generated for all ids of the first scenario
#-----------------------------------------------------------------------------------------------------------------------

#mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/GBG/plot_pca
#mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/GBB/plot_pca


#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s*.skov/regionsperid/GBB/regi_pca/gbb.allids.allintrogregions.pca.27jul22.pdf  /Users/harvi/Downloads/gorilla_postabc/skov/GBB
# works & looks fine

# perform the same for GBG, & the 40kb reps (these are pcas of all introg regions or filtered by length >=40kb), then for the pcas of regions private to a given individual
#-----------------------------------------------------------------------------------------------------------------------

# 2) read in data for GBG .pca


scen=(c(rep("GBG",9)))
# test
SPE="GBG"
spe="gbg"

pca_obj=list()
for (i in 1:length(scen)){
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/pca/",spe,".",i,".pca",sep=""),verbose=T)
pca_obj[[i]]<-tes
}

random_pca_obj=list()
for (i in 1:length(scen)){
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/random/",SPE,"/pca/",spe,".",i,".random.pca",sep=""),verbose=T)
random_pca_obj[[i]]<-tes
}

# obtain plot limits (xlim, ylim)
lims_12<-check_limits_pca(pca_obj,1,2,random_pca_obj)
lims_34<-check_limits_pca(pca_obj,3,4,random_pca_obj)

hold_plots=list()
for (i in 1:length(scen)){
test_12<-plot_pca(pca_obj,i,3,4,1,2,lims_12)
test_34<-plot_pca(pca_obj,i,5,6,3,4,lims_34)
rand_12<-random_plot_pca(random_pca_obj,i,3,4,1,2,lims_12)
rand_34<-random_plot_pca(random_pca_obj,i,5,6,3,4,lims_34)
hold_plots[[i]]<-ggarrange(test_12,rand_12,test_34,rand_34, ncol = 2, nrow = 2)
}


SPE="GBG"
spe="gbg"
pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/privateregionsperid/",SPE,"/plot_pca/",spe,".allids.privateintrogregions.pca.30sep22.pdf",sep="")) 
hold_plots
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

q()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/'overlap.s*.skov'/privateregionsperid/GBG/plot_pca/gbg.allids.privateintrogregions.pca.30sep22.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/pca
scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/'overlap.s*.skov'/privateregionsperid/GBB/plot_pca/gbb.allids.privateintrogregions.pca.30sep22.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/pca
