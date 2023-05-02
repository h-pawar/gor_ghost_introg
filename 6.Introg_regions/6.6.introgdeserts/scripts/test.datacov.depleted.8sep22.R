# Thu  8 Sep 2022 09:52:14 CEST
# assess data coverage in the regions identified as depleted for fragments of introgression

#-----------------------------------------------------------------------------------------------------------------------

# callable regions
#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)
files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/", pattern="gorilla_"), sep = "")

genomewide=list()
for (i in 1:length(files)){
load(file=files[[i]],verbose=T)
genomewide[[i]]<-aweight3
}
factortonumeric_fun<-function(x){
genomewide[[x]][,c(1:3)] <- sapply(genomewide[[x]][,c(1:3)] ,as.character)
genomewide[[x]][,c(1:3)] <- sapply(genomewide[[x]][,c(1:3)] ,as.numeric)
return(genomewide[[x]])
}

genomewide_num=list()
for (j in 1:22){
genomewide_num[[j]]<-factortonumeric_fun(j)
}

test<-data.frame(do.call(rbind.data.frame, genomewide_num))
colnames(test)<-c('chr','startpos','proportion')

test<-test[order(test$chr),]

tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/callable.proportion.windows.txt")
write.table(test,tmp,sep="\t",row.names=F,col.names=T,quote=F)

#-----------------------------------------------------------------------------------------------------------------------
# intersect this with the deserts, to annotate the depleted regions with the callable proportion


## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1 
require(data.table)
library(GenomicRanges)
library(Biostrings)
library("ggbio")
library('GenomeInfoDb')

load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/reduced.deserts.23aug22",verbose=T)
# unique depleted regions - generated in test.deserts.22aug22.R


#-----------------------------------------------------------------------------------------------------------------------

props.windows<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/callable.proportion.windows.txt")

colnames(props.windows)<-props.windows[1,]

props.windows<-props.windows[-1,]

props.windows[,] <- sapply(props.windows[,],as.numeric)

props.windows$endpos<-props.windows$startpos+40000
props.windows<-props.windows[,c(1,2,4,3)]

props.windows$chr<-sub("^","chr",props.windows$chr)

G_inf<-GRanges(seqnames=props.windows[,1],ranges=IRanges(start=as.numeric(props.windows[,2]),end=as.numeric(props.windows[,3])),strand=rep("*",length(props.windows[,1])))
mcols(G_inf)$proportion<-props.windows$proportion

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
hits<-findOverlaps(G_inf, reduced_deserts, ignore.strand=TRUE)
x1<-unique(as.data.frame(hits)[,1])
test<-G_inf[x1]
overlaps<- as.data.frame(hits)


calc_prop_deserts_fun<-function(it){
start<-min(overlaps[(overlaps$"subjectHits"==it),][,1])
end<-max(overlaps[(overlaps$"subjectHits"==it),][,1])
all_prop<-G_inf[start:end]$proportion
mean_prop<-mean(G_inf[start:end]$proportion)
return(list(all_prop,mean_prop))}


hold_prop<-list()
for (i in (1:length(reduced_deserts))) {
hold_prop[[i]]<-calc_prop_deserts_fun(i) }


hold_means<-list()
for (i in (1:length(reduced_deserts))) {
hold_means[[i]]<-hold_prop[[i]][[2]] }


# filter out by callable sites, then by length of depleted region

index_keep<-which(unlist(hold_means)>0.5)
reduced_deserts[index_keep] # the depleted regions with sufficient callable sites

 
# filter by those >=5mb as the depleted regions
reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=5000000]

# total length of the retained depleted regions
sum(width(reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=5000000]))/1000000

red_deserts_0.5_5mb<-reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=5000000]
red_deserts_0.5_8mb<-reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=8000000]

#-----------------------------------------------------------------------------------------------------------------------
# plotting

library(ggbio)
data(ideoCyto, package = "biovizBase")


all_hg19<-ideoCyto$hg19
sex_chr<-all_hg19[all_hg19@seqnames=="chrX" | all_hg19@seqnames=="chrY"]
autosomes<-setdiff(all_hg19,sex_chr)

seqlengths(red_deserts_0.5_5mb)<-seqlengths(autosomes)[1:22]
seqlengths(red_deserts_0.5_8mb)<-seqlengths(autosomes)[1:22]


pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/plot.deserts.0.5.5mb.9sep22.pdf")
autoplot(red_deserts_0.5_5mb, layout = "karyogram") + labs(title = ">= 5Mb")
autoplot(red_deserts_0.5_8mb, layout = "karyogram") + labs(title = ">= 8Mb")


dev.off()

#-----------------------------------------------------------------------------------------------------------------------

# check coverage of the longest deserts retained

calc_prop_deserts_fun_1<-function(des){
hits<-findOverlaps(G_inf, des, ignore.strand=TRUE)
x1<-unique(as.data.frame(hits)[,1])
test<-G_inf[x1]
overlaps<- as.data.frame(hits)

calc_prop_deserts_fun<-function(it){
start<-min(overlaps[(overlaps$"subjectHits"==it),][,1])
end<-max(overlaps[(overlaps$"subjectHits"==it),][,1])
all_prop<-G_inf[start:end]$proportion
mean_prop<-mean(G_inf[start:end]$proportion)
return(list(all_prop,mean_prop))}


hold_prop<-list()
for (i in (1:length(des))) {
hold_prop[[i]]<-calc_prop_deserts_fun(i) }

return(hold_prop)
}

prop_8<-calc_prop_deserts_fun_1(red_deserts_0.5_8mb)
prop_5<-calc_prop_deserts_fun_1(red_deserts_0.5_5mb)


proc_means_fun<-function(p){
hold_means<-list()
for (i in (1:length(p))) {
hold_means[[i]]<-p[[i]][[2]] }
return(as.data.frame(unlist(hold_means)))}

mcols(red_deserts_0.5_8mb)$proportion<-unlist(proc_means_fun(prop_8))
mcols(red_deserts_0.5_5mb)$proportion<-unlist(proc_means_fun(prop_5))


depleted_0.5_5mb_8mb<-list(red_deserts_0.5_5mb,red_deserts_0.5_8mb)
save(depleted_0.5_5mb_8mb,file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/depleted.0.5.5mb.8mb")



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


calc_prop_deserts_fun<-function(i,t){
x<-paste("chr",i,sep="")
sum(width(t[t@seqnames==x]))/width(autosomes)[i]
}

calc_prop_deserts_perchr<-function(t){
hold_prop<-list()
for (i in (1:22)) {
hold_prop[[i]]<-calc_prop_deserts_fun(i,t) }
return(hold_prop)}

calc_prop_deserts_perchr(red_deserts_0.5_5mb) # for regions >=5mb
calc_prop_deserts_perchr(red_deserts_0.5_8mb) # for regions >=8mb
