# Thu 27 Oct 2022 11:41:33 CEST
# generate haplotype networks for the intersect regions (s*-skov)
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(mgcv)
library(GenomicRanges)


load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/longwindows",sep=""),verbose=T) 

input<-longwindows[[4]]
myGRangesList<-GRangesList(input)
reduced <- reduce(unlist(myGRangesList))
reduced[reduced@seqnames=="chr9",] 
chr9_regions<-which(reduced@seqnames=="chr9")
aut<-reduced[-c(53:56)]


non_overlapping_region_counts<-function(x){
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
return(reducedConsensus)
}

non_overlapping_region_counts(myGRangesList)
rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)]))
length(rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)])))
				
overlaps_counts<-non_overlapping_region_counts(myGRangesList)

extract_regions_fun<-function(){
regions<-as.data.frame(overlaps_GBG)
	
rout <- strsplit(as.character(regions$seqnames),'chr') 
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

for (i in (1:nrow(regions))) {
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun()

chr9_df<-as.data.frame(overlaps_counts)[chr9_regions,]

aut_df<-as.data.frame(overlaps_counts)[-chr9_regions,]

#-----------------------------------------------------------------------------------------------------------------------

# function for the autosomes
extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbg250k")


chr9_extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbg250k")


#-----------------------------------------------------------------------------------------------------------------------

#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


nput="gbg250k"

dir=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/",sep="")

vcffiles <- paste(dir, list.files(path =dir, pattern=".vcf"), sep = "")


chr9_vcffiles<-vcffiles[88:91]
vcffiles<-vcffiles[1:87]

vcffiles1<-c(vcffiles,chr9_vcffiles)


# for all autosomes (minus chr 9)

hhpnet<-list()
hind.hap<-list()
for(i in (1:length(vcffiles))) {
print(i)
locs<-read.vcfR(file=vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[i]]<-try(haploNet(hploc))
hind.hap[[i]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}

# will need to read in the chr 9 regions separately, b/c instead of GBB 12, its 1 GBG, 12 GBB, 8 GBG etc

# for chr 9
for(i in (1:length(chr9_vcffiles))) {
index<-seq(88:91)+87
ptype<-index[[i]]
print(i)
locs<-read.vcfR(file=chr9_vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBG",1),rep("GBB",12),rep("GBG",8),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[ptype]]<-try(haploNet(hploc))
hind.hap[[ptype]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}


t<- strsplit(as.character(vcffiles[[1]]), 'gbg250k_')

t1<-strsplit(as.character(t[[1]][2]), '_')
strsplit(as.character(t1[[1]]), '.overlaps')


format_region_fun<-function(i) {
t<-strsplit(as.character(vcffiles1[[i]]), 'gbg250k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
region<-unlist(t2)[1:3]
return(region)
}


hregions<-list()
for(i in (1:length(hhpnet))) {
hregions[[i]]<-format_region_fun(i)
}



pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/plot/",nput,"_overlaps.hapl.pdf",sep=""))
for(i in (1:length(hhpnet))) {
print(i)
region=hregions[[i]]
plot(hhpnet[[i]], size=attr(hhpnet[[i]], "freq"), scale.ratio = 2, cex = 0.8, pie=hind.hap[[i]],show.mutation=3,main=paste(region,sep="-"))
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20)
}

dev.off()
}


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# next scenario, gbb250k, index 2


input<-longwindows[[2]]
myGRangesList<-GRangesList(input)
reduced <- reduce(unlist(myGRangesList))

reduced[reduced@seqnames=="chr9",] 

chr9_regions<-which(reduced@seqnames=="chr9")

aut<-reduced[-c( chr9_regions)]


non_overlapping_region_counts<-function(x){
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
return(reducedConsensus)
}



rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)]))

length(rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)])))

overlaps_counts<-non_overlapping_region_counts(myGRangesList)


chr9_df<-as.data.frame(overlaps_counts)[chr9_regions,]

aut_df<-as.data.frame(overlaps_counts)[-chr9_regions,]

#-----------------------------------------------------------------------------------------------------------------------

# function for the autosomes
extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbb250k")



chr9_extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 
regions$chr<-as.numeric(do.call(rbind, rout)[,2]) 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbb250k")


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


nput="gbb250k"

dir=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/",sep="")

vcffiles <- paste(dir, list.files(path =dir, pattern=".vcf"), sep = "")


chr9_vcffiles<-vcffiles[100:103]
vcffiles<-vcffiles[1:99]

vcffiles1<-c(vcffiles,chr9_vcffiles)


# for all autosomes (minus chr 9)

hhpnet<-list()
hind.hap<-list()
for(i in (1:length(vcffiles))) {
print(i)
locs<-read.vcfR(file=vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[i]]<-try(haploNet(hploc))
hind.hap[[i]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}

# for chr 9
for(i in (1:length(chr9_vcffiles))) {
index<-seq(100:103)+99
ptype<-index[[i]]
print(i)
locs<-read.vcfR(file=chr9_vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBG",1),rep("GBB",12),rep("GBG",8),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[ptype]]<-try(haploNet(hploc))
hind.hap[[ptype]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}


format_region_fun<-function(i) {
t<-strsplit(as.character(vcffiles1[[i]]), 'gbb250k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
region<-unlist(t2)[1:3]
return(region)
}


hregions<-list()
for(i in (1:length(hhpnet))) {
hregions[[i]]<-format_region_fun(i)
}



pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/plot/",nput,"_overlaps.hapl.pdf",sep=""))
for(i in (1:length(hhpnet))) {
print(i)
region=hregions[[i]]
plot(hhpnet[[i]], size=attr(hhpnet[[i]], "freq"), scale.ratio = 2, cex = 0.8, pie=hind.hap[[i]],show.mutation=3,main=paste(region,sep="-"))
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) 
}

dev.off()


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# next scenario, gbb100k, index 1

input<-longwindows[[1]]
myGRangesList<-GRangesList(input)
reduced <- reduce(unlist(myGRangesList))

reduced[reduced@seqnames=="chr9",] 


chr9_regions<-which(reduced@seqnames=="chr9")

aut<-reduced[-c( chr9_regions)]


non_overlapping_region_counts<-function(x){
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
return(reducedConsensus)
}


overlaps_counts<-non_overlapping_region_counts(myGRangesList)


chr9_df<-as.data.frame(overlaps_counts)[chr9_regions,]

aut_df<-as.data.frame(overlaps_counts)[-chr9_regions,]

#-----------------------------------------------------------------------------------------------------------------------

# function for the autosomes
extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbb100k")

chr9_extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbb100k")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)

nput="gbb100k"

dir=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/",sep="")

vcffiles <- paste(dir, list.files(path =dir, pattern=".vcf"), sep = "")


chr9_vcffiles<-vcffiles[604:628]
vcffiles<-vcffiles[1:603]

vcffiles1<-c(vcffiles,chr9_vcffiles)


# for all autosomes (minus chr 9)

hhpnet<-list()
hind.hap<-list()
for(i in (1:length(vcffiles))) {
print(i)
locs<-read.vcfR(file=vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[i]]<-try(haploNet(hploc))
hind.hap[[i]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}



# for chr 9
for(i in (1:length(chr9_vcffiles))) {
index<-seq(604:628)+603
ptype<-index[[i]]
print(i)
locs<-read.vcfR(file=chr9_vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBG",1),rep("GBB",12),rep("GBG",8),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[ptype]]<-try(haploNet(hploc))
hind.hap[[ptype]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}


format_region_fun<-function(i) {
t<-strsplit(as.character(vcffiles1[[i]]), 'gbb100k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
region<-unlist(t2)[1:3]
return(region)
}


hregions<-list()
for(i in (1:length(hhpnet))) {
hregions[[i]]<-format_region_fun(i)
}



pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/plot/",nput,"_overlaps.hapl.pdf",sep=""))
for(i in (1:length(hhpnet))) {
print(i)
region=hregions[[i]]
plot(hhpnet[[i]], size=attr(hhpnet[[i]], "freq"), scale.ratio = 2, cex = 0.8, pie=hind.hap[[i]],show.mutation=3,main=paste(region,sep="-"))
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20)
}

dev.off()


#-----------------------------------------------------------------------------------------------------------------------

# next scenario, gbg100k, index 1


input<-longwindows[[3]]

myGRangesList<-GRangesList(input)
reduced <- reduce(unlist(myGRangesList))


reduced[reduced@seqnames=="chr9",] 


chr9_regions<-which(reduced@seqnames=="chr9")

aut<-reduced[-c( chr9_regions)]


non_overlapping_region_counts<-function(x){
reduced <- reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
return(reducedConsensus)
}


overlaps_counts<-non_overlapping_region_counts(myGRangesList)


chr9_df<-as.data.frame(overlaps_counts)[chr9_regions,]

aut_df<-as.data.frame(overlaps_counts)[-chr9_regions,]

#-----------------------------------------------------------------------------------------------------------------------

# function for the autosomes
extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbg100k")

chr9_extract_regions_fun<-function(regions,spe){

rout <- strsplit(as.character(regions$seqnames),'chr') 

regions$chr<-as.numeric(do.call(rbind, rout)[,2])


for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbg100k")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


nput="gbg100k"

dir=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/",sep="")

vcffiles <- paste(dir, list.files(path =dir, pattern=".vcf"), sep = "")


chr9_vcffiles<-vcffiles[473:490]
vcffiles<-vcffiles[1:472]

vcffiles1<-c(vcffiles,chr9_vcffiles)


# for all autosomes (minus chr 9)

hhpnet<-list()
hind.hap<-list()
for(i in (1:length(vcffiles))) {
print(i)
locs<-read.vcfR(file=vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[i]]<-try(haploNet(hploc))
hind.hap[[i]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}



# for chr 9
for(i in (1:length(chr9_vcffiles))) {
index<-seq(473:490)+472
ptype<-index[[i]]
print(i)
locs<-read.vcfR(file=chr9_vcffiles[[i]])
blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
rownames(blocs)<-c(rep("GBG",1),rep("GBB",12),rep("GBG",8),rep("GGD",1),rep("GGG",27))
hploc<-haplotype(blocs)
hhpnet[[ptype]]<-try(haploNet(hploc))
hind.hap[[ptype]]<-with(
  stack(setNames(attr(hploc, "index"), rownames(hploc))),
  table(hap=ind, pop=rownames(blocs)[values])
)
}



format_region_fun<-function(i) {
t<-strsplit(as.character(vcffiles1[[i]]), 'gbg100k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')

region<-unlist(t2)[1:3]
return(region)
}


hregions<-list()
for(i in (1:length(hhpnet))) {
hregions[[i]]<-format_region_fun(i)
}



pdf(paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/plot/",nput,"_overlaps.hapl.pdf",sep=""))
for(i in (1:length(hhpnet))) {
print(i)
region=hregions[[i]]
plot(hhpnet[[i]], size=attr(hhpnet[[i]], "freq"), scale.ratio = 2, cex = 0.8, pie=hind.hap[[i]],show.mutation=3,main=paste(region,sep="-"))
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) 
}

dev.off()



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
