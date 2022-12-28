# Thu 27 Oct 2022 11:41:33 CEST
# generate haplotype networks for the intersect regions (s*-skov) (querying the correct vcfs)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# identified in this section, that the old chr 9 was being queried -> go to the next section, where this issue is resolved

# having errors in this step - go through interactively - Wed 27 Jul 2022 10:48:09 CEST

# now read in & plot this - changing R versions **
#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

#require(data.table)
#options(stringsAsFactors=F)
#library(tidyr)
#library('pegas')
#library(vcfR)
#options("scipen"=100)


# interactivley for nput="gbg250k"
# go over all regions

# generate the haplotypes
#generate_hapl<-function(nput) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))

#dir=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",nput,"/",sep="")

#vcffiles <- paste(dir, list.files(path =dir, pattern=".vcf"), sep = "")

#length(vcffiles)
#[1] 91

# write as function & plot for these 8

#hhpnet<-list()
#hind.hap<-list()
#for(i in (1:length(vcffiles))) {
#	print(i)
#locs<-read.vcfR(file=vcffiles[[i]])
#blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
#rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
#hploc<-haplotype(blocs)
#hhpnet[[i]]<-try(haploNet(hploc))
#hind.hap[[i]]<-with(
#  stack(setNames(attr(hploc, "index"), rownames(hploc))),
#  table(hap=ind, pop=rownames(blocs)[values])
#)
#}

# is this getting stuck in one of the regions?
# stuck at 88..(is not proceedign to 89..)

# run specifically for region 89

#> vcffiles[[89]]
#[1] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_134758000_135051000.overlaps.vcf"

# whether this was generated properly, ie perhaps was querying the old tumani rather thna the new tumani **
	# go back over the code & amend *


#agt<-system(paste("/apps/BCFTOOLS/1.9/bin/bcftools view -u -f '%POS %REF;%ALT [%GT,] [%DP,] [%MQ,] [%MQ0,]\n' /scratch/devel/hpawar/admix/sstar/vcf_24jun22/segsite_filt",ft,"_",chrom,".vcf.gz",sep=""),intern=T)

#testsnps<-system(paste("bcftools query -f '%POS [%GT ]\n' -r chr9:134758000-135051000 /scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_134758000_135051000.overlaps.vcf",sep=""),intern=T)

# but htis does have genotypes..

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Gorilla_beringei_beringei-Bwiruka       Gorilla_beringei_bering
#ei-Imfura        Gorilla_beringei_beringei-Kaboko        Gorilla_beringei_beringei-Kahungye      Gorilla_beringei_beringei-Katungi
#       Gorilla_beringei_beringei-Maisha        Gorilla_beringei_beringei-Nyamunwa      Gorilla_beringei_beringei-Semehe        Gorilla_
#beringei_beringei-Tuck  Gorilla_beringei_beringei-Turimaso      Gorilla_beringei_beringei-Umurimo       Gorilla_beringei_beringei-Zirik
#ana      Gorilla_beringei_graueri-9732_Mkubwa    Gorilla_beringei_graueri-A929_Kaisi     Gorilla_beringei_graueri-A967_Victoria  Gorill
#a_beringei_graueri-Dunia  Gorilla_beringei_graueri-Itebero        Gorilla_beringei_graueri-Mukokya        Gorilla_beringei_graueri-Ntabwoba       Gorilla_beringei_graueri-Pinga  Gorilla_beringei_graueri-Tumani Gorilla_gorilla_diehli-B646_Nyango      Gorilla_gorilla_gorilla-9749_Kowali     Gorilla_gorilla_gorilla-9750_Azizi      Gorilla_gorilla_gorilla-9751_Bulera     Gorilla_gorilla_gorilla-9752_Suzie
#      Gorilla_gorilla_gorilla-9753_Kokomo     Gorilla_gorilla_gorilla-A930_Sandra     Gorilla_gorilla_gorilla-A931_Banjo      Gorilla_gorilla_gorilla-A932_Mimi       Gorilla_gorilla_gorilla-A933_Dian       Gorilla_gorilla_gorilla-A934_Delphi     Gorilla_gorilla_gorilla-A936_Coco       Gorilla_gorilla_gorilla-A937_Kolo       Gorilla_gorilla_gorilla-A962_Amani      Gorilla_gorilla_gorilla-B642_Akiba_Beri Gorilla_gorilla_gorilla-B643_Choomba    Gorilla_gorilla_gorilla-B644_Paki       Gorilla_gorilla_gorilla-B647_Anthal     Gorilla_gorilla_gorilla-B650_Katie      Gorilla_gorilla_gorilla-KB3782_Vila     Gorilla_gorilla_gorilla-KB3784_Dolly    Gorilla_gorilla_gorilla-KB4986_Katie    Gorilla_gorilla_gorilla-KB5792_Carolyn  Gorilla_gorilla_gorilla-KB5852_Helen    Gorilla_gorilla_gorilla-KB6039_Oko      Gorilla_gorilla_gorilla-KB7973_Porta    Gorilla_gorilla_gorilla-X00108_Abe      Gorilla_gorilla_gorilla-X00109_Tzambo
#chr9    134928730       .       T       C       148.77  PASS    .       GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0
#     0|0     0|1     0|0     0|0     0|1     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0
#chr9    134928752       .       T       C       1020.77 PASS    .       GT      1|1     1|1     1|1     1|1     1|1     1|1     1|1
#     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1
#chr9    134928786       .       C       T       1031.77 PASS    .       GT      1|1     1|1     1|1     1|1     1|1     1|1     1|1
#     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1
#chr9    134958888       .       C       T       280.78  PASS    .       GT      0|0     0|0     0|0     0|0     0|0     0|0     0|0
#     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     0|0     1|1     0|1     1|1     0|0     0|0     0|0     0|0     0|1     0|1     0|1     0|1     0|1     0|0     0|0     0|0     0|0     0|0     1|1     0|0     0|1     0|0     0|0     0|0     0|1     0|0     0|0     0|0
#chr9    134958903       .       A       G       967.77  PASS    .       GT      1|1     1|1     1|1     1|1     1|1     1|1     1|1
#     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1     1|1

# need to check if this was queried with the right vcf *



#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#vcffiles<- c("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_134758000_135051000.overlaps.vcf")

#hhpnet<-list()
#hind.hap<-list()
#for(i in (1:length(vcffiles))) {
#locs<-read.vcfR(file=vcffiles[[i]])
#blocs<-vcfR2DNAbin(locs,consensus=T,extract.haps=F)
#rownames(blocs)<-c(rep("GBB",12),rep("GBG",9),rep("GGD",1),rep("GGG",27))
#hploc<-haplotype(blocs)
#hhpnet[[i]]<-try(haploNet(hploc))
#hind.hap[[i]]<-with(
#  stack(setNames(attr(hploc, "index"), rownames(hploc))),
#  table(hap=ind, pop=rownames(blocs)[values])
#)
#}

# even when running interactively with only this region, hangs 

#Scanning file to determine attributes.
#File attributes:
#  meta lines: 127
#  header_line: 128
#  variant count: 5
#  column count: 58
#Meta line 127 read in.
#All meta lines processed.
#gt matrix initialized.
#Character matrix gt created.
#  Character matrix gt rows: 5
#  Character matrix gt cols: 58
#  skip: 0
#  nrows: 5
#  row_num: 0
#Processed variant: 5
#All variants processed
#After extracting indels, 5 variants remain.
# this doesnt make sense

# best strategy - go back to region output & rewrite this - yes, executed below
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#Wed 26 Oct 2022 12:51:48 CEST
## module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.12
# have removed previous vcf files & plots from these long windows -> starting again, querying the correct vcf (ie newly generated chr 9 for regions on chr 9)

## module load gcc/6.3.0 openssl/1.0.2q  R/4.0.1  BCFTOOLS/1.12
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
options(scipen=100)
library(adegenet)
library(mgcv)
library(GenomicRanges)

# load this file
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/longwindows",sep=""),verbose=T) 

#longwindows<-list(l_gbb_100k,l_gbb_250k,l_gbg_100k,l_gbg_250k)


#longwindows[[4]] # gbg 250k

# AMEND BELOW FUNCTION, TO ACCOUNT FOR NEWLY GENERATED CHR 9


# B) haplotype networks

# 1) for GBG
# process regions for the haplotypes
#proc_hapl<-function(input,spe) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))

# longestwindowsperid for this scenario

# interactively:
input<-longwindows[[4]]
# check

# check this
myGRangesList<-GRangesList(input) # convert to GRangesList object of length equal to number of individuals
reduced <- reduce(unlist(myGRangesList))
# check per id if fine 

# check which regions are on chr 9


reduced[reduced@seqnames=="chr9",] 


#GRanges object with 4 ranges and 0 metadata columns:
#      seqnames              ranges strand
#         <Rle>           <IRanges>  <Rle>
#  [1]     chr9   31735000-32063000      *
#  [2]     chr9   97658000-97957000      *
#  [3]     chr9 120006000-120288000      *
#  [4]     chr9 134758000-135051000      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths
#which(reduced@seqnames=="chr9")
#[1] 53 54 55 56


chr9_regions<-which(reduced@seqnames=="chr9")


# length(reduced[-c(53:56)])
# [1] 87

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

# for all the granges regions
non_overlapping_region_counts(myGRangesList)
#GRanges object with 91 ranges and 10 metadata columns:
#       seqnames              ranges strand |        V1        V2        V3
#          <Rle>           <IRanges>  <Rle> | <numeric> <numeric> <numeric>
#   [1]     chr1   48152000-48407000      * |         0         0         0
#   [2]     chr1   77372000-77714000      * |         1         0         0
#   [3]     chr1   86933000-87380000      * |         0         0         1




rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)]))

length(rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)])))

#> rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)]))
# [1] 2 1 2 1 1 1 2 3 6 1 1 1 2 7 1 1 1 1 2 1 1 1 1 2 4 1 2 1 1 1 1 1 1 1 1 1 1 1
#[39] 2 1 6 1 3 1 1 2 3 3 4 1 1 2 1 1 1 2 2 1 2 1 4 1 2 4 2 1 4 1 3 1 1 1 2 3 2 1
#[77] 1 1 1 2 1 1 4 1 1 1 1 2 3 1 4
#> 
#> length(rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)])))
#[1] 91

# ie extract the row sums corresponding to 53:56

overlaps_counts<-non_overlapping_region_counts(myGRangesList)

extract_regions_fun<-function(){

regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun()


# head(as.data.frame(overlaps_counts))
#  seqnames     start       end  width strand V1 V2 V3 V4 V5 V6 V7 V8 V9
#1     chr1  48152000  48407000 255001      *  0  0  0  0  0  0  1  0  1
#2     chr1  77372000  77714000 342001      *  1  0  0  0  0  0  0  0  0


#as.data.frame(overlaps_counts)[chr9_regions,]
#   seqnames     start       end  width strand V1 V2 V3 V4 V5 V6 V7 V8 V9
#53     chr9  31735000  32063000 328001      *  0  0  0  0  0  1  0  0  0
#54     chr9  97658000  97957000 299001      *  0  0  1  0  0  0  0  0  0
#55     chr9 120006000 120288000 282001      *  0  0  0  1  0  0  0  0  0
#56     chr9 134758000 135051000 293001      *  1  0  0  0  0  0  1  0  0
#   consensusIDs
#53 consensus_53
#54 consensus_54
#55 consensus_55
#56 consensus_56

# minus these from the df

# & query these separately

chr9_df<-as.data.frame(overlaps_counts)[chr9_regions,]

#nrow(as.data.frame(overlaps_counts)[-chr9_regions,])
#[1] 87
# for rest of the autosomes



aut_df<-as.data.frame(overlaps_counts)[-chr9_regions,]

#-----------------------------------------------------------------------------------------------------------------------

# function for the autosomes
extract_regions_fun<-function(regions,spe){

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbg250k")

# all generated



chr9_extract_regions_fun<-function(regions,spe){

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbg250k")

# have all generated for gbg 250k

#-----------------------------------------------------------------------------------------------------------------------
# go to the next step & plot


#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


# go over all regions

# generate the haplotypes
#generate_hapl<-function(nput) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))


# run interactively : 
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
#[88] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_120006000_120288000.overlaps.vcf" 
#[89] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_134758000_135051000.overlaps.vcf" 
#[90] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_31735000_32063000.overlaps.vcf"   
#[91] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/gbg250k_chr9_97658000_97957000.overlaps.vcf"   

# ie for 88-91 will need to amend 

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

Warning messages:
1: In haplotype.DNAbin(blocs) :
  some sequences of different lengths were assigned to the same haplotype
2: In haplotype.DNAbin(blocs) :
  some sequences were not assigned to the same haplotype because of ambiguities
3: In haplotype.DNAbin(blocs) :
  some sequences of different lengths were assigned to the same haplotype
4: In haplotype.DNAbin(blocs) :
  some sequences were not assigned to the same haplotype because of ambiguities
5: In haplotype.DNAbin(blocs) :
  some sequences of different lengths were assigned to the same haplotype
6: In haplotype.DNAbin(blocs) :
  some sequences were not assigned to the same haplotype because of ambiguities
7: In haplotype.DNAbin(blocs) :
  some sequences of different lengths were assigned to the same haplotype
8: In haplotype.DNAbin(blocs) :
  some sequences were not assigned to the same haplotype because of ambiguities



   strsplit(as.character(vcffiles[[1]]), 'gbg250k_')
[[1]]
[1] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbg250k/"
[2] "chr10_101840000_102104000.overlaps.vcf"                       


t<- strsplit(as.character(vcffiles[[1]]), 'gbg250k_')

t1<-strsplit(as.character(t[[1]][2]), '_')
strsplit(as.character(t1[[1]]), '.overlaps')

#> strsplit(as.character(t1[[1]]), '.overlaps')
#[[1]]
#[1] "chr10"

#[[2]]
#[1] "101840000"

#[[3]]
#[1] "102104000" ".vcf"  


# now obtain the regions
format_region_fun<-function(i) {
#t<-strsplit(as.character(vcffiles[[i]]), 'GGG_')
t<-strsplit(as.character(vcffiles1[[i]]), 'gbg250k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
#region<-unlist(t2)
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
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) ## needs to be changed
}

dev.off()
}


# check how these look like, then run interactively for the next scenario *

# copy to local 

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/hapl/tmp/gbg250k/plot/gbg250k_overlaps.hapl.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots

#-----------------------------------------------------------------------------------------------------------------------


# generate for the second scenario *


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#longwindows<-list(l_gbb_100k,l_gbb_250k,l_gbg_100k,l_gbg_250k)


#longwindows[[4]] # gbg 250k

# AMEND BELOW FUNCTION, TO ACCOUNT FOR NEWLY GENERATED CHR 9


# B) haplotype networks

# 1) for GBG
# process regions for the haplotypes
#proc_hapl<-function(input,spe) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))

# longestwindowsperid for this scenario


# next scenario, gbb250k, index 2

# interactively:
input<-longwindows[[2]]
# check

# check this
myGRangesList<-GRangesList(input) # convert to GRangesList object of length equal to number of individuals
reduced <- reduce(unlist(myGRangesList))
# check per id if fine 

# check which regions are on chr 9


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
 # [1] 1 2 1 1 2 2 1 1 1 4 4 2 1 4 1 1 2 5 3 1 2 1 3 1 1 3 5 1 4 2 2 1 2 3 1 2 1
# [38] 2 3 3 2 1 3 2 1 1 1 2 1 2 1 4 1 2 6 2 1 1 1 2 3 1 1 1 4 1 1 3 1 1 1 2 2 1
# [75] 4 3 1 1 7 4 4 3 1 2 5 1 4 3 6 2 1 1 1 2 3 1 1 5 3 1 2 3 4

length(rowSums(as.matrix(mcols(non_overlapping_region_counts(myGRangesList))[,1:length(input)])))
#[1] 103

overlaps_counts<-non_overlapping_region_counts(myGRangesList)


chr9_df<-as.data.frame(overlaps_counts)[chr9_regions,]

aut_df<-as.data.frame(overlaps_counts)[-chr9_regions,]

#-----------------------------------------------------------------------------------------------------------------------

# function for the autosomes
extract_regions_fun<-function(regions,spe){

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbb250k")

# all generated



chr9_extract_regions_fun<-function(regions,spe){

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbb250k")

# have all generated for gbg 250k
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# go to the next step & plot


#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


# go over all regions

# generate the haplotypes
#generate_hapl<-function(nput) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))


# run interactively : 
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

# will need to read in the chr 9 regions separately, b/c instead of GBB 12, its 1 GBG, 12 GBB, 8 GBG etc
#[88] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbb250k/gbb250k_chr9_120006000_120288000.overlaps.vcf" 
#[89] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbb250k/gbb250k_chr9_134758000_135051000.overlaps.vcf" 
#[90] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbb250k/gbb250k_chr9_31735000_32063000.overlaps.vcf"   
#[91] "/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/gbb250k/gbb250k_chr9_97658000_97957000.overlaps.vcf"   

# ie for 88-91 will need to amend 

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

# same warning messages, only for chr 9
#Warning messages:
#1: In haplotype.DNAbin(blocs) :
#  some sequences of different lengths were assigned to the same haplotype
#2: In haplotype.DNAbin(blocs) :
#  some sequences were not assigned to the same haplotype because of ambiguities
#3: In haplotype.DNAbin(blocs) :
#  some sequences of different lengths were assigned to the same haplotype
#4: In haplotype.DNAbin(blocs) :
#  some sequences were not assigned to the same haplotype because of ambiguities
#5: In haplotype.DNAbin(blocs) :
#  some sequences of different lengths were assigned to the same haplotype
#6: In haplotype.DNAbin(blocs) :
#  some sequences were not assigned to the same haplotype because of ambiguities
#7: In haplotype.DNAbin(blocs) :
#  some sequences of different lengths were assigned to the same haplotype
#8: In haplotype.DNAbin(blocs) :
#  some sequences were not assigned to the same haplotype because of ambiguities


# now obtain the regions
format_region_fun<-function(i) {
#t<-strsplit(as.character(vcffiles[[i]]), 'GGG_')
t<-strsplit(as.character(vcffiles1[[i]]), 'gbb250k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
#region<-unlist(t2)
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
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) ## needs to be changed
}

dev.off()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/hapl/tmp/gbb250k/plot/gbb250k_overlaps.hapl.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#longwindows<-list(l_gbb_100k,l_gbb_250k,l_gbg_100k,l_gbg_250k)


#longwindows[[4]] # gbg 250k


# next scenario, gbb100k, index 1

# interactively:
input<-longwindows[[1]]
# check

# check this
myGRangesList<-GRangesList(input) # convert to GRangesList object of length equal to number of individuals
reduced <- reduce(unlist(myGRangesList))
# check per id if fine 

# check which regions are on chr 9


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

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbb100k")

# all generated



chr9_extract_regions_fun<-function(regions,spe){

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbb100k")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# go to the next step & plot


#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


# go over all regions

# generate the haplotypes
#generate_hapl<-function(nput) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))


# run interactively : 
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


# now obtain the regions
format_region_fun<-function(i) {
#t<-strsplit(as.character(vcffiles[[i]]), 'GGG_')
t<-strsplit(as.character(vcffiles1[[i]]), 'gbb100k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
#region<-unlist(t2)
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
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) ## needs to be changed
}

dev.off()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/hapl/tmp/gbb100k/plot/gbb100k_overlaps.hapl.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots
# with 100k regions, is perhaps too heavy...



#-----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#longwindows<-list(l_gbb_100k,l_gbb_250k,l_gbg_100k,l_gbg_250k)


#longwindows[[4]] # gbg 250k


# next scenario, gbg100k, index 1

# interactively:
input<-longwindows[[3]]
# check

# check this
myGRangesList<-GRangesList(input) # convert to GRangesList object of length equal to number of individuals
reduced <- reduce(unlist(myGRangesList))
# check per id if fine 

# check which regions are on chr 9


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

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",regions$chr[i],".vcf.gz > ",out1," ",sep=""),intern=T)
}

}

extract_regions_fun(aut_df,"gbg100k")

# all generated



chr9_extract_regions_fun<-function(regions,spe){

#regions<-as.data.frame(overlaps_GBG)

rout <- strsplit(as.character(regions$seqnames),'chr') 
#do.call(rbind, rout)
regions$chr<-as.numeric(do.call(rbind, rout)[,2])

# generate vcfs for these regions 

for (i in (1:nrow(regions))) {
print(i)
out1<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/hapl/tmp/",spe,"/",spe,"_",regions$seqnames[i],"_",regions$start[i],"_",regions$end[i],".overlaps.vcf",sep="")
system(paste("bcftools view -r ",regions$seqnames[i],":",regions$start[i],"-",regions$end[i]," /scratch/devel/hpawar/admix/sstar/vcf_24jun22/3_gor.chr9.vcf.gz > ",out1," ",sep=""),intern=T)
}

}

chr9_extract_regions_fun(chr9_df,"gbg100k")

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# go to the next step & plot


#module unload R/4.0.1
#module load gcc/6.3.0 R/3.4.2 BCFTOOLS/1.12

require(data.table)
options(stringsAsFactors=F)
library(tidyr)
library('pegas')
library(vcfR)
options("scipen"=100)


# go over all regions

# generate the haplotypes
#generate_hapl<-function(nput) {

#cn1<-list(c("GBG","GBB","GGG","GGD"))


# run interactively : 
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


# now obtain the regions
format_region_fun<-function(i) {
#t<-strsplit(as.character(vcffiles[[i]]), 'GGG_')
t<-strsplit(as.character(vcffiles1[[i]]), 'gbg100k_')
t1<- strsplit(as.character(t[[1]][2]), '_')
t2<-strsplit(as.character(t1[[1]]), '.overlaps')
#region<-unlist(t2)
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
legend("bottomright", c("GBB","GBG","GGD","GGG"), col=rainbow(ncol(hind.hap[[i]])), pch=20) ## needs to be changed
}

dev.off()

scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/hapl/tmp/gbg100k/plot/gbg100k_overlaps.hapl.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
