# Mon 22 Aug 2022 18:47:30 BST
# recalculate deserts  - excluding short windows at end of chr & excluding centromeres
#-----------------------------------------------------------------------------------------------------------------------

#Q - whether to exclude the short windows at the end of the chr?
#Q - centromere localisation? whether to include this in the counts? (or in the denominators for chr size?)

# MK
#1. Yes, you may exclude windows <1Mbp; yes, centromeres should be
#excluded, do you have coordinates?

# MK re centromeres
#I used the ggbio package to plot the genome. 
#It includes info on the human genome, but since I can't access the cluster anymore, I can't check the details. 
#Not even sure the package still works fine.
#-----------------------------------------------------------------------------------------------------------------------

#library("ggbio")
#data(hg19IdeogramCyto, package = "biovizBase")
#hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X")))
#cols<-getOption("biovizBase")$cytobandColor
#cols[which(names(cols) %in% c("acen"))]<-"black" ## "acen" is the centromeres

#ggplot(hg19) +  layout_karyogram( cytoband = T,aes(linetype="blank",xlab=""))

#Anyhow, this file should contain the same cytobands including centromeres ("acen"):
#http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz

#This is also a rather strict definition of centromeres, I guess...

#(base) harvi@Harvinders-MacBook-Pro ~ % head ~/Downloads/cytoBand.txt
#chr1	0	2300000	p36.33	gneg
#chr1	2300000	5400000	p36.32	gpos25
#chr1	5400000	7200000	p36.31	gneg
#chr1	7200000	9200000	p36.23	gpos25

#(base) harvi@Harvinders-MacBook-Pro ~ % grep 'acen'  ~/Downloads/cytoBand.txt
#chr1	121500000	125000000	p11.1	acen
#chr1	125000000	128900000	q11	acen

#-----------------------------------------------------------------------------------------------------------------------

# https://www.rdocumentation.org/packages/GWASTools/versions/1.18.0/topics/centromeres
#BiocManager::install("GWASTools")
#Centromere base positions from the GRCh36/hg18, GRCh37/hg19 and GRCh38/hg38 genome builds.
#Usage
#data(centromeres.hg18)
#data(centromeres.hg19)

# gives
#A data frame with the following columns.
#chromchromosome (1-22, X, Y) left.basestarting base position of centromere right.baseending base position of centromere

#> (centromeres.hg19)
#   chrom left.base right.base
#1      1 121535434  124535434
#2      2  92326171   95326171
#-----------------------------------------------------------------------------------------------------------------------



## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1 
require(data.table)
library(GenomicRanges)
library(Biostrings)
library("ggbio")
library('GenomeInfoDb')

#-----------------------------------------------------------------------------------------------------------------------
# location of centromeres
data(hg19IdeogramCyto, package = "biovizBase")
# shoudl remove chr x - only considering the autosomes
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22)),pruning.mode="coarse")
centromeres <- hg19[hg19$gieStain=='acen']

#centromeres
#GRanges object with 44 ranges and 2 metadata columns:
#       seqnames              ranges strand |     name gieStain
#          <Rle>           <IRanges>  <Rle> | <factor> <factor>
#   [1]     chr1 121500000-125000000      * |   p11.1      acen
#   [2]     chr1 125000000-128900000      * |   q11        acen
#   [3]    chr10   38000000-40200000      * |   p11.1      acen
#   [4]    chr10   40200000-42300000      * |   q11.1      acen
#   [5]    chr11   51600000-53700000      * |   p11.11     acen
#   ...      ...                 ...    ... .      ...      ...
#  [40]     chr7   59900000-61700000      * |    q11.1     acen
#  [41]     chr8   43100000-45600000      * |    p11.1     acen
#  [42]     chr8   45600000-48100000      * |    q11.1     acen
#  [43]     chr9   47300000-49000000      * |    p11.1     acen
#  [44]     chr9   49000000-50700000      * |    q11       acen
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths

#-----------------------------------------------------------------------------------------------------------------------


# non-overlapping windows across chromosomes

mygenome <- readDNAStringSet("/home/devel/marcmont/scratch/snpCalling_hg19/chimp/assembly/BWA/hg19.fa")
chrSizes <- width(mygenome)
names(chrSizes) <- names(mygenome)
#print(chrSizes)

# or directly read in the lengths of hg19 chr:
	# from overlap.random.regions.20jul22.R
#-----------------------------------------------------------------------------------------------------------------------
# lengths of hg19
#hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

#hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

# hg19 chr in the correct order
#hg19<-hg19[c(1:18,20,19,22,21),]
#colnames(hg19)<-c("chrom","size") 
# or whether need to read in the genome itself?
  # path to hg19 ref

#> hg19
    # checking chrSizes[1:22] & hg19 are equivalent # yes
#-----------------------------------------------------------------------------------------------------------------------


bins   <- tileGenome(chrSizes[1:22], tilewidth=1000000, cut.last.tile.in.chrom=T)
#print(bins)

# 1,000,000 bp length bins across the genome
 	# gives non-overlapping 1mb windows per chr, but at the end of the chr (the window is shorter)
# print(bins)
#GRanges object with 2897 ranges and 0 metadata columns:
#         seqnames            ranges strand
#            <Rle>         <IRanges>  <Rle>
#     [1]     chr1         1-1000000      *
#     [2]     chr1   1000001-2000000      *
#     [3]     chr1   2000001-3000000      *
#     [4]     chr1   3000001-4000000      *
#     [5]     chr1   4000001-5000000      *
#     ...      ...               ...    ...
#  [2893]    chr22 47000001-48000000      *
#  [2894]    chr22 48000001-49000000      *
#  [2895]    chr22 49000001-50000000      *
#  [2896]    chr22 50000001-51000000      *
#  [2897]    chr22 51000001-51304566      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome

#-----------------------------------------------------------------------------------------------------------------------

# find bins of length  <1000000  (ie at end of each chr)
#> bins[width(bins)<1000000]
#GRanges object with 22 ranges and 0 metadata columns:
#       seqnames              ranges strand
#          <Rle>           <IRanges>  <Rle>
#   [1]     chr1 249000001-249250621      *
#   [2]     chr2 243000001-243199373      *
#   [3]     chr3 198000001-198022430      *
#   [4]     chr4 191000001-191154276      *
#   [5]     chr5 180000001-180915260      *
#   ...      ...                 ...    ...
#  [18]    chr18   78000001-78077248      *
#  [19]    chr19   59000001-59128983      *
#  [20]    chr20   63000001-63025520      *
#  [21]    chr21   48000001-48129895      *
#  [22]    chr22   51000001-51304566      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome

#> bins[width(bins)==1000000]
#GRanges object with 2875 ranges and 0 metadata columns:
  
# only retain 1mb regions
bins<-bins[width(bins)==1000000]

#-----------------------------------------------------------------------------------------------------------------------
# remove centromeres - after the overlap? (at end of process) centromeres <- hg19[hg19$gieStain=='acen']
# or remove windows containing a centromere? before overlapping the bins with the regions?

# find overlaps b/n bins & centromeres & then remove windows containing a centromere? - or remove the centromere regions after calc the overlap? - if any of the overlap regions == a centromere
# & then again remove windows < 1mb in length
	# try this approach first - remove bins which overlap w centromere - then remove < 1mb bins, & then follow the below
#-----------------------------------------------------------------------------------------------------------------------


# subtract B from A
# A[!A %over% B]

# subtract centromeres granges from bins
#bins[!bins %over% centromeres]

 bins[!bins %over% centromeres]
#GRanges object with 2779 ranges and 0 metadata columns:
#         seqnames            ranges strand
#            <Rle>         <IRanges>  <Rle>
#     [1]     chr1         1-1000000      *
#     [2]     chr1   1000001-2000000      *
#     [3]     chr1   2000001-3000000      *
#     [4]     chr1   3000001-4000000      *
#     [5]     chr1   4000001-5000000      *
#     ...      ...               ...    ...
#  [2775]    chr22 47000001-48000000      *
#  [2776]    chr22 48000001-49000000      *
#  [2777]    chr22 49000001-50000000      *
#  [2778]    chr22 50000001-51000000      *
#  [2779]    chr22 51000001-51304566      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome

#2875-2779
#[1] 96
# loss of 96 windows

bins<- bins[!bins %over% centromeres]


# & again only retain 1mb regions
bins<-bins[width(bins)==1000000]
#GRanges object with 2757 ranges and 0 metadata columns:

#-----------------------------------------------------------------------------------------------------------------------

# query the bins by the (overlap of s*-skov) windows per eastern individual

 # intersect the non-overlapping windows with the overlap of (s*-skov)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb40_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg40_99",verbose=T)

#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------

# take into account multiple introg fragments in one individual may overlap with any given bin (since bins are of large size)

count_regions_fun<-function(scen,ind){
# overlap the 1mb bins with the putative introg fragments per individual
# will need to go over all ids of ov_gbb_99 & of ov_gbg_99
hits<-findOverlaps(bins, scen[[ind]], ignore.strand=TRUE)
x1<-unique(as.data.frame(hits)[,1])
# where there is a putative introg fragment for the individual
test<-bins[x1]
tmp<-(as.data.frame(hits)[,1]) # indices of the bins
r_freq<-table(tmp) # freq of introg regions overlapping each bin
# add identifer (1 putative introg fragment in this bin, or multiple)
mcols(test)$id<-r_freq
# take non-overlapping ranges, add 0s for id_1
empty<-bins[!(bins %over% test)]
# 0 = no putative introg fragment in this bin for this individual
mcols(empty)$id<-0
# list these 2 grange objects (test & empty), then merge into one object 
grlist <- GenomicRanges::GRangesList(test,empty)
tmp = unlist(grlist)
# sort the ranges numerically per chromosome - 
# verify that seqlevels are sorted:
tmp <- sortSeqlevels(tmp)
# sort your GRanges object:
tmp <- sort(tmp)
# output tmp
return(tmp)
}


count_across_ids_fun<-function(scen){
hold_counts<-list()
for (i in (1:length(scen))) {
# [,1]: overlapping bp  
hold_counts[[i]]<-count_regions_fun(scen,i)}
return(hold_counts)
}


f_gbb<-count_across_ids_fun(ov_gbb_99)
f_gbg<-count_across_ids_fun(ov_gbg_99)

# extract only the metadata (0s and 1s (2s) per individual)

process_fun<-function(scen){
hold_id<-list()
for (i in (1:length(scen))) {
# [,1]: overlapping bp  
hold_id[[i]]<-as.data.frame(scen[[i]]$id) }
out<-do.call(cbind,hold_id)
return(out)
}

g_gbb<-process_fun(f_gbb)

g_gbg<-process_fun(f_gbg)

g_e<-cbind(g_gbb,g_gbg)
g_e1<-cbind(as.data.frame(bins)[,c(1:4)],g_e)

#-----------------------------------------------------------------------------------------------------------------------

deserts<-g_e1[rowSums(g_e)==0,]
# convert deserts back to granges object
fG<-GRanges(seqnames=deserts[,1],ranges=IRanges(start=as.numeric(deserts[,2]),end=as.numeric(deserts[,3]),names=deserts[,2]),strand=rep("*",length(deserts[,1])))
reduce(fG)
sum(width(reduce(fG)))/sum(width(bins)) # deserts / bp of bins
#[1] 0.5437069 # divison by the bp of the genome covered by the bins


#sum(chrSizes[1:22])
#[1] 2881033286

sum(width(reduce(fG)))/ sum(chrSizes[1:22]) # deserts / size of whole genome
#[1] 0.5202994

sum(width(reduce(fG)))/1000000
# [1] 1499 # Mb # still seems rather high

 sum(width(bins))/1000000 # size of the bins covering the genome (excluding centromeres & windows < 1mb)
#[1] 2757


# ie half does not contain any introg signature

#-----------------------------------------------------------------------------------------------------------------------
# unique deserts
reduce(fG)
GRanges object with 580 ranges and 0 metadata columns:
        seqnames            ranges strand
           <Rle>         <IRanges>  <Rle>
    [1]     chr1         1-4000000      *
    [2]     chr1   6000001-9000000      *
    [3]     chr1 11000001-14000000      *
    [4]     chr1 15000001-24000000      *
    [5]     chr1 26000001-31000000      *
    ...      ...               ...    ...
  [576]    chr22 33000001-34000000      *
  [577]    chr22 37000001-38000000      *
  [578]    chr22 39000001-42000000      *
  [579]    chr22 45000001-47000000      *
  [580]    chr22 48000001-51000000      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

# freq for each window (but the number of windows eg for == 0 is inflated, as some are consecutive windows, as seen from the reduce)
#table(rowSums(g_e))
#  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
#1499  338  226  149  121   89   58   58   42   40   29   24   11   19    7   10 
#  16   17   18   19   20   21   25   26   32   33 
#   8    8    4    4    6    2    2    1    1    1 

# could rather break down the regions per chr for g_e1?
	# or generate further granges of the bin positions, widths & the rowsums
		# this may be best for downstream ?

# continue from here

# ie need to assess proportion per chromosome
# & if these depleted regions are correct,  generate random sets of equal length distributions 

# split the unique deserts per chr & calc proportion depleted per chr 


t<-reduce(fG) # unique deserts


reduced_deserts<-reduce(fG) 
# write out this object, in order to then generate random regions of equal dist
#mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/deserts
save(reduced_deserts,file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/reduced.deserts.23aug22")


calc_prop_deserts_fun<-function(i){
x<-paste("chr",i,sep="")
sum(width(t[t@seqnames==x]))/chrSizes[i]
}

#chrSizes[1:22]
#> calc_prop_deserts_fun("chr1",1)
#     chr1 
#0.5987974 
# q - whether am overestimating this??

hold_prop<-list()
for (i in (1:22)) {
hold_prop[[i]]<-calc_prop_deserts_fun(i) }

# still high proportions of deserts per chr...
hold_prop
[[1]]
     chr1 
0.5656957 

[[2]]
     chr2 
0.5427646 

[[3]]
     chr3 
0.4999434 

[[4]]
    chr4 
0.402816 

[[5]]
     chr5 
0.4864156 

[[6]]
     chr6 
0.4266135 

[[7]]
     chr7 
0.5215577 

[[8]]
    chr8 
0.464595 

[[9]]
     chr9 
0.6373331 

[[10]]
 chr10 
0.4058 

[[11]]
    chr11 
0.5259005 

[[12]]
    chr12 
0.4333147 

[[13]]
    chr13 
0.5296524 

[[14]]
    chr14 
0.6520755 

[[15]]
    chr15 
0.6534584 

[[16]]
    chr16 
0.5201719 

[[17]]
    chr17 
0.5788519 

[[18]]
    chr18 
0.4354662 

[[19]]
    chr19 
0.7272237 

[[20]]
    chr20 
0.5553306 

[[21]]
    chr21 
0.6233132 

[[22]]
    chr22 
0.5847433 


# still high proportions inferred

# chr 19 also the most depleted here (same seen in bonobos - kuhlwilm suppl pg 11)

#-----------------------------------------------------------------------------------------------------------------------

# convert g_e1 back to granges
G_e1<-GRanges(seqnames=g_e1[,1],ranges=IRanges(start=as.numeric(g_e1[,2]),end=as.numeric(g_e1[,3]),names=g_e1[,2]),strand=rep("*",length(g_e1[,1])))
mcols(G_e1)$popfreq<-rowSums(g_e)

> G_e1
GRanges object with 2757 ranges and 1 metadata column:
           seqnames            ranges strand |   popfreq
              <Rle>         <IRanges>  <Rle> | <numeric>
         1     chr1         1-1000000      * |         0
   1000001     chr1   1000001-2000000      * |         0
   2000001     chr1   2000001-3000000      * |         0
   3000001     chr1   3000001-4000000      * |         0
   4000001     chr1   4000001-5000000      * |         3
       ...      ...               ...    ... .       ...
  46000001    chr22 46000001-47000000      * |         0
  47000001    chr22 47000001-48000000      * |         1
  48000001    chr22 48000001-49000000      * |         0
  49000001    chr22 49000001-50000000      * |         0
  50000001    chr22 50000001-51000000      * |         0
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

# use this to generate the karyotype plot of fragment density **
	# will need to explore how to do this **


# go through, re help with plotting : 
#https://bioconductor.statistik.tu-dortmund.de/packages/2.11/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
# pg 87 - could use to generate where the deserts are, where multiple fragments across each chr 
# potentially also pg 226

#-----------------------------------------------------------------------------------------------------------------------

nondeserts<-g_e1[rowSums(g_e)!=0,]

G_nondeserts<-GRanges(seqnames=nondeserts[,1],ranges=IRanges(start=as.numeric(nondeserts[,2]),end=as.numeric(nondeserts[,3]),names=nondeserts[,2]),strand=rep("*",length(nondeserts[,1])))

mcols(G_nondeserts)$popfreq<-rowSums(nondeserts[,-c(1:4)])
> G_nondeserts
GRanges object with 1258 ranges and 1 metadata column:
           seqnames            ranges strand |   popfreq
              <Rle>         <IRanges>  <Rle> | <numeric>
   4000001     chr1   4000001-5000000      * |         3
   5000001     chr1   5000001-6000000      * |         3
   9000001     chr1  9000001-10000000      * |         4
  10000001     chr1 10000001-11000000      * |         2
  14000001     chr1 14000001-15000000      * |         1
       ...      ...               ...    ... .       ...
  38000001    chr22 38000001-39000000      * |         2
  42000001    chr22 42000001-43000000      * |         1
  43000001    chr22 43000001-44000000      * |         8
  44000001    chr22 44000001-45000000      * |         8
  47000001    chr22 47000001-48000000      * |         1
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

q()
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# Tue 23 Aug 2022 15:29:10 BST

# MK
#the introgression deserts as defined in the bonobo paper are consecutive regions of at least 8Mbp (Fig. 5a), 
#with some additional information on depleted regions (>5Mbp, Fig. 5b). 
#For Neandertals, they used 10Mbp. So, you should look for such long regions, maybe ask how many you find with 15, 10, 8, 5 Mbp.
#The deserts might be shorter in bonobos due to larger Ne than non-African humans. 
#Then, they might be longer in gorillas, but would be interesting to know.

#-----------------------------------------------------------------------------------------------------------------------


# assess width of deserts
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/reduced.deserts.23aug22",verbose=T)

# frequency table of deserts inferred
> table(width(reduced_deserts))

 1000000  2000000  3000000  4000000  5000000  6000000  7000000  8000000 
     240      154       66       38       29       21        5        8 
 9000000 10000000 11000000 12000000 14000000 15000000 16000000 17000000 
       6        3        1        2        2        1        2        1 
20000000 
       1 

# filter by those >=5mb as the depleted regions

# kuhlwilm 2019: sliding windows of 5 Mbp in 1-Mbp steps # alternate way of calculating
# Putative introgression deserts (>8 Mbp)


reduced_deserts1<-reduced_deserts[width(reduced_deserts)>=5000000]

#sum(width(reduced_deserts1))
#[1] 601000000

sum(width(reduced_deserts1))/1000000
[1] 601 # Mb # for depleted regions

sum(width(reduced_deserts1))/ sum(chrSizes[1:22]) # depleted regions (>=5mb) / size of whole genome
[1] 0.2086057

reduced_deserts2<-reduced_deserts[width(reduced_deserts)>=8000000]
sum(width(reduced_deserts2))/ sum(chrSizes[1:22])
[1] 0.1023938

sum(width(reduced_deserts2))/1000000
[1] 295 # for putative deserts (followign the bonobo definitions)

# depends on definition of a desert



calc_prop_deserts_fun<-function(i,t){
x<-paste("chr",i,sep="")
sum(width(t[t@seqnames==x]))/chrSizes[i]
}

calc_prop_deserts_perchr<-function(t){
hold_prop<-list()
for (i in (1:22)) {
hold_prop[[i]]<-calc_prop_deserts_fun(i,t) }
return(hold_prop)}

# could plot both of these & distribution of deserts across the genome
calc_prop_deserts_perchr(reduced_deserts1) # for regions >=5mb
calc_prop_deserts_perchr(reduced_deserts2) # for regions >=8mb


> calc_prop_deserts_perchr(reduced_deserts1)
[[1]]
     chr1 
0.3329982 

[[2]]
     chr2 
0.1726978 

[[3]]
     chr3 
0.2524966 

[[4]]
      chr4 
0.02615688 

[[5]]
     chr5 
0.1437137 

[[6]]
      chr6 
0.06428423 

[[7]]
     chr7 
0.1633795 

[[8]]
     chr8 
0.1161488 

[[9]]
     chr9 
0.2974221 

[[10]]
    chr10 
0.1106727 

[[11]]
    chr11 
0.2740608 

[[12]]
     chr12 
0.08965133 

[[13]]
    chr13 
0.2257535 

[[14]]
    chr14 
0.3819299 

[[15]]
    chr15 
0.3998775 

[[16]]
    chr16 
0.3098896 

[[17]]
    chr17 
0.2586359 

[[18]]
chr18 
    0 

[[19]]
    chr19 
0.6088385 

[[20]]
    chr20 
0.1586659 

[[21]]
    chr21 
0.4155421 

[[22]]
    chr22 
0.2338973 

> calc_prop_deserts_perchr(reduced_deserts2)
[[1]]
     chr1 
0.1685051 

[[2]]
      chr2 
0.08223705 

[[3]]
     chr3 
0.1666478 

[[4]]
chr4 
   0 

[[5]]
     chr5 
0.0442196 

[[6]]
chr6 
   0 

[[7]]
chr7 
   0 

[[8]]
chr8 
   0 

[[9]]
     chr9 
0.2549333 

[[10]]
chr10 
    0 

[[11]]
    chr11 
0.1703621 

[[12]]
chr12 
    0 

[[13]]
    chr13 
0.1389252 

[[14]]
    chr14 
0.2328841 

[[15]]
    chr15 
0.2243215 

[[16]]
    chr16 
0.1992148 

[[17]]
  chr17 
0.12316 

[[18]]
chr18 
    0 

[[19]]
    chr19 
0.3213314 

[[20]]
chr20 
    0 

[[21]]
    chr21 
0.2077711 

[[22]]
    chr22 
0.2338973 