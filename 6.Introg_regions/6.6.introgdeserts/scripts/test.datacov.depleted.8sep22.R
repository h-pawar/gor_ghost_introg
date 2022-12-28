# Thu  8 Sep 2022 09:52:14 CEST
# assess data coverage in the regions identified as depleted for fragments of introgression
# to do this need to write out bed file of the proportion of callable regions for each 40kb window
# then assess coverage & subsequently filter

# for depleted regions following test.deserts.22aug22.R, test.deserts.2aug22.R 

#Fri  2 Sep 2022 17:12:30 BST
# MK
#What I suggest is to take each 1Mbp bin and ask how much of it is covered with data
# (that is, the cumulative amount of 40kbp windows within the bin). 
# So, for each 1Mbp region, you have two values:
#1) introgressed fragment density, 
#2) proportion of the whole region covered; and you could subset those regions with >0.5 or other values.


#-----------------------------------------------------------------------------------------------------------------------

# amend depleted.in.props.windowsregions.2sep22.R - to extract proportion of callable regions per window

# Fri  2 Sep 2022 08:18:16 BST
# need to filter deserts by sufficiently props.windows windows **

# the window info is in R v 3.4.2, whereas the introgressed regions are in R 4.0.1
# output the props.windows windows as bed file - to be readinto R 4.0.1 as granges object

# filter.het.by.clbl.26aug.R:

#/scratch/devel/hpawar/admix/abc/emp.data/window.info/ : info per chr, of proportion of sites per window, after filtering by repeats & mappability

#module load gcc/6.3.0 R/3.4.2
require(data.table)
options(stringsAsFactors=F)
files <- paste("/scratch/devel/hpawar/admix/abc/emp.data/window.info/", list.files(path = "/scratch/devel/hpawar/admix/abc/emp.data/window.info/", pattern="gorilla_"), sep = "")

#load(file=files[[1]],verbose=T)
#x<-aweight3


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

# loop through
genomewide_num=list()
for (j in 1:22){
genomewide_num[[j]]<-factortonumeric_fun(j)
}

#str(genomewide_num)
#List of 22
# $ :'data.frame': 3389 obs. of  3 variables:
#  ..$ chr                                         : num [1:3389] 10 10 10 10 10 10 10 10 10 10 ...
#  ..$ V2                                          : num [1:3389] 0 40000 80000 120000 160000 200000 240000 280000 320000 360000 ...
#  ..$ format(aweight2, digits = 3, scientific = F): num [1:3389] 0 0.0195 0.0838 0.4378 0.5427 ...
# $ :'data.frame': 3376 obs. of  3 variables:
#  ..$ chr                                         : num [1:3376] 11 11 11 11 11 11 11 11 11 11 ...
#  ..$ V2                                          : num [1:3376] 0 40000 80000 120000 160000 200000 240000 280000 320000 360000 ...
#  ..$ format(aweight2, digits = 3, scientific = F): num [1:3376] 0 0.000725 0.007975 0.029925 0.2069 ...


test<-data.frame(do.call(rbind.data.frame, genomewide_num))
colnames(test)<-c('chr','startpos','proportion')

#head(test)
#  chr startpos proportion
#1  10        0   0.000000
#2  10    40000   0.019475

# order by chr
test<-test[order(test$chr),]

# write this out
tmp<-paste("/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/callable.proportion.windows.txt")
write.table(test,tmp,sep="\t",row.names=F,col.names=T,quote=F)

# MK: apply a threshold of 50%, to keep a large enough part of the genome : ie object 'test' here was filtered by 0.5 for calculation of empirical summary stats

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


#reduced_deserts
#GRanges object with 580 ranges and 0 metadata columns:
#        seqnames            ranges strand
#           <Rle>         <IRanges>  <Rle>
#    [1]     chr1         1-4000000      *
#    [2]     chr1   6000001-9000000      *
#    [3]     chr1 11000001-14000000      *
#    [4]     chr1 15000001-24000000      *
#    [5]     chr1 26000001-31000000      *
#    ...      ...               ...    ...
#  [576]    chr22 33000001-34000000      *
#  [577]    chr22 37000001-38000000      *
#  [578]    chr22 39000001-42000000      *
#  [579]    chr22 45000001-47000000      *
#  [580]    chr22 48000001-51000000      *
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths


# lengths of the depleted regions
#> table(width(reduced_deserts))
# 1000000  2000000  3000000  4000000  5000000  6000000  7000000  8000000 
#     240      154       66       38       29       21        5        8 
# 9000000 10000000 11000000 12000000 14000000 15000000 16000000 17000000 
#       6        3        1        2        2        1        2        1 
#20000000 
#       1 


# then find overlaps

#-----------------------------------------------------------------------------------------------------------------------

props.windows<-read.table("/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/callable.proportion.windows.txt")

colnames(props.windows)<-props.windows[1,]

props.windows<-props.windows[-1,]

# convert cols from character to  numeric
props.windows[,] <- sapply(props.windows[,],as.numeric)


# 40kb windows
props.windows$endpos<-props.windows$startpos+40000
props.windows<-props.windows[,c(1,2,4,3)]

# add chr to col1
props.windows$chr<-sub("^","chr",props.windows$chr)

#> head(props.windows )
#   chr startpos endpos proportion
#2 chr1        0  40000   0.014450
#3 chr1    40000  80000   0.060700
#4 chr1    80000 120000   0.029875
#5 chr1   120000 160000   0.028725
#6 chr1   160000 200000   0.013975
#7 chr1   200000 240000   0.011800
#> str(props.windows)
#'data.frame':   72036 obs. of  4 variables:
# $ chr       : chr  "chr1" "chr1" "chr1" "chr1" ...
# $ startpos  : num  0 40000 80000 120000 160000 200000 240000 280000 320000 360000 ...
# $ endpos    : num  40000 80000 120000 160000 200000 240000 280000 320000 360000 400000 ...
# $ proportion: num  0.0144 0.0607 0.0299 0.0287 0.014 ...


# convert to granges
G_inf<-GRanges(seqnames=props.windows[,1],ranges=IRanges(start=as.numeric(props.windows[,2]),end=as.numeric(props.windows[,3])),strand=rep("*",length(props.windows[,1])))
mcols(G_inf)$proportion<-props.windows$proportion

#-----------------------------------------------------------------------------------------------------------------------

# overlap the depleted regions with the 40kb windows

#hits<-findOverlaps(G_inf, reduced_deserts, ignore.strand=TRUE)

#> hits
#Hits object with 38055 hits and 0 metadata columns:
#          queryHits subjectHits
#          <integer>   <integer>
#      [1]         1           1
#      [2]         2           1
#      [3]         3           1
#      [4]         4           1
#      [5]         5           1
#      ...       ...         ...
#  [38051]     72025         580
#  [38052]     72026         580
#  [38053]     72027         580
#  [38054]     72028         580
#  [38055]     72029         580
#  -------
#  queryLength: 72036 / subjectLength: 580

#x1<-unique(as.data.frame(hits)[,1])

#test<-G_inf[x1]

# test
#GRanges object with 38055 ranges and 1 metadata column:
#          seqnames            ranges strand | proportion
#             <Rle>         <IRanges>  <Rle> |  <numeric>
#      [1]     chr1           0-40000      * |   0.014450
#      [2]     chr1       40000-80000      * |   0.060700
#      [3]     chr1      80000-120000      * |   0.029875
#      [4]     chr1     120000-160000      * |   0.028725
#      [5]     chr1     160000-200000      * |   0.013975
#      ...      ...               ...    ... .        ...
#  [38051]    chr22 50840000-50880000      * |   0.547150
#  [38052]    chr22 50880000-50920000      * |   0.765725
#  [38053]    chr22 50920000-50960000      * |   0.579125
#  [38054]    chr22 50960000-51000000      * |   0.593000
#  [38055]    chr22 51000000-51040000      * |   0.509575
#  -------
#  seqinfo: 22 sequences from an unspecified genome; no seqlengths

# these 40kb windows overlap with a depleted region
# match these proportions to each depleted region


#overlaps<- as.data.frame(hits)
# overlaps[(overlaps$"subjectHits"==1),] # for the first depleted region, the 40kb windows covering expand from window 1-101
#100       100           1
#101       101           1


#min(overlaps[(overlaps$"subjectHits"==1),][,1])
#max(overlaps[(overlaps$"subjectHits"==1),][,1])

# test[1:101]
#GRanges object with 101 ranges and 1 metadata column:
#        seqnames          ranges strand | proportion
#           <Rle>       <IRanges>  <Rle> |  <numeric>
#    [1]     chr1         0-40000      * |   0.014450
#    [2]     chr1     40000-80000      * |   0.060700
##    [3]     chr1    80000-120000      * |   0.029875
#    [4]     chr1   120000-160000      * |   0.028725
#    [5]     chr1   160000-200000      * |   0.013975
#    ...      ...             ...    ... .        ...


# test[1:101]$proportion
#  test[1:101]$proportion
#  [1] 0.014450 0.060700 0.029875 0.028725 0.013975 0.011800 0.021350 0.000000
#  [9] 0.001900 0.001575 0.001675 0.000375 0.000000 0.028975 0.067900 0.008475
# [17] 0.022475 0.055775 0.192300 0.320450 0.324450 0.841050 0.816850 0.510450
# [25] 0.821300 0.590875 0.586475 0.796950 0.633350 0.736500 0.674425 0.860075
# [33] 0.620600 0.464400 0.554675 0.275250 0.457150 0.513375 0.511400 0.356075
# [41] 0.276250 0.195775 0.599075 0.530500 0.362525 0.482200 0.671750 0.559575
# [49] 0.410725 0.629950 0.709350 0.556075 0.668075 0.634825 0.668475 0.801925
# [57] 0.902725 0.603375 0.735825 0.912025 0.971175 0.897025 0.537775 0.532975
# [65] 0.484850 0.229675 0.000000 0.399225 0.612300 0.737125 0.746275 0.672875
# [73] 0.591325 0.696250 0.716500 0.902775 0.876600 0.827550 0.893725 0.832875
# [81] 0.807275 0.790925 0.872150 0.849425 0.699425 0.863500 0.846400 0.780950
# [89] 0.740875 0.652300 0.851075 0.660300 0.589825 0.538325 0.496600 0.619125
# [97] 0.077500 0.000000 0.000000 0.024275 0.493700

 # ie includes regions with 0 callable regions -> will need to split into separate regions
# or whether to take the mean proportion across this single region? # coudl do with discussing this
# perhaps take the mean as a starting point? - go with this, can't think of a good alternative..

#mean(test[1:101]$proportion)
#[1] 0.4965569

# whether shoudl go back to 1mb bins?

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
hits<-findOverlaps(G_inf, reduced_deserts, ignore.strand=TRUE)
x1<-unique(as.data.frame(hits)[,1])
test<-G_inf[x1]
overlaps<- as.data.frame(hits)


# function assess start & end window (of the 40kb windows) which overlap with a given depleted region (for the 580 depleted regions - loose definition >= 1Mb)
# apply to each of the 580 regions


calc_prop_deserts_fun<-function(it){
start<-min(overlaps[(overlaps$"subjectHits"==it),][,1])
end<-max(overlaps[(overlaps$"subjectHits"==it),][,1])
all_prop<-G_inf[start:end]$proportion
mean_prop<-mean(G_inf[start:end]$proportion)
return(list(all_prop,mean_prop))}


hold_prop<-list()
for (i in (1:length(reduced_deserts))) {
hold_prop[[i]]<-calc_prop_deserts_fun(i) }

# now gives output
# shoudl clearly remove 'depleted regions' where all 40kb windows contained within have 0 callable sites
# shoudl establish a threshold? or should be able to filter out from the mean vals? 


#[[572]]
#[[572]][[1]]
#  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [38] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [75] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[112] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[149] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[186] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[223] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[260] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#[297] 0 0 0 0 0

#[[572]][[2]]
#[1] 0

# first extract means across the bins -> assess distribution -> filter out those with v low data coverage
hold_means<-list()
for (i in (1:length(reduced_deserts))) {
hold_means[[i]]<-hold_prop[[i]][[2]] }


# table(unlist(hold_means)) # distribution of mean callable proportions for each 'depleted region'
# filter out by callable sites, then by length of depleted region

#> summary(unlist(hold_means))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.4514  0.4885  0.4767  0.5200  0.6767 

#which(unlist(hold_means)<0.5) # these windows to remove

# assess lengths of depleted regions to keep

index_keep<-which(unlist(hold_means)>0.5)

reduced_deserts[index_keep] # the depleted regions with sufficient callable sites

# freq dist of lengths of these regions
#> table(width(reduced_deserts[index_keep]))

# 1000000  2000000  3000000  4000000  5000000  6000000  7000000  8000000 
#     111       64       17       14       11        6        1        1 
# 9000000 12000000 
#       1        1 
 
# filter by those >=5mb as the depleted regions

reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=5000000] # retains 21 regions

#GRanges object with 21 ranges and 0 metadata columns:
#       seqnames              ranges strand
#          <Rle>           <IRanges>  <Rle>
#   [1]     chr1   78000001-84000000      *
#   [2]     chr2   38000001-45000000      *
#   [3]     chr2   55000001-67000000      *
#   [4]     chr2 206000001-211000000      *
#   [5]     chr2 238000001-243000000      *
#   ...      ...                 ...    ...
#  [17]    chr15   87000001-95000000      *
#  [18]    chr16           1-5000000      *
#  [19]    chr17   76000001-81000000      *
#  [20]    chr20   17000001-22000000      *
#  [21]    chr21   22000001-27000000      *
#  -------
 # seqinfo: 22 sequences from an unspecified genome; no seqlength

# total length of the retained depleted regions
sum(width(reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=5000000]))/1000000
[1] 127 # 127 Mb - fewer than seen in bonobos...

# plot these as a distribution

red_deserts_0.5_5mb<-reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=5000000]
red_deserts_0.5_8mb<-reduced_deserts[index_keep][width(reduced_deserts[index_keep])>=8000000]

#> sum(width(red_deserts_0.5_5mb))/1000000
#[1] 127
#> sum(width(red_deserts_0.5_8mb))/1000000
#[1] 29

#-----------------------------------------------------------------------------------------------------------------------
# plotting- following plot.deserts.25aug22.R

library(ggbio)
data(ideoCyto, package = "biovizBase")

# karyogram lengths hg19
all_hg19<-ideoCyto$hg19
sex_chr<-all_hg19[all_hg19@seqnames=="chrX" | all_hg19@seqnames=="chrY"]
autosomes<-setdiff(all_hg19,sex_chr)

#p <- autoplot(autosomes, layout = "karyogram") + theme_clear()
#autoplot(autosomes, layout = "karyogram") + theme_clear() + layout_karyogram(reduced_deserts1,aes(color = "grey", fill = "grey"))

# set sequence lengths
seqlengths(red_deserts_0.5_5mb)<-seqlengths(autosomes)[1:22]
seqlengths(red_deserts_0.5_8mb)<-seqlengths(autosomes)[1:22]


#http://bioconductor.org/packages/devel/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/plot.deserts.0.5.5mb.9sep22.pdf")
autoplot(red_deserts_0.5_5mb, layout = "karyogram") + labs(title = ">= 5Mb")
autoplot(red_deserts_0.5_8mb, layout = "karyogram") + labs(title = ">= 8Mb")


dev.off()

#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/overlap.s\*.skov/deserts/plot.deserts.0.5.5mb.9sep22.pdf /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov/plots

# at the >= 5mb still looks like some regions at the end of the chr - check coords & proportions of remaining depleted regions
# & how much of the chr is depleted 

#depleted_0.5_5mb_8mb<-list(red_deserts_0.5_5mb,red_deserts_0.5_8mb)
#save(depleted_0.5_5mb_8mb,file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/depleted.0.5.5mb.8mb")

red_deserts_0.5_8mb
GRanges object with 3 ranges and 0 metadata columns:
      seqnames            ranges strand
         <Rle>         <IRanges>  <Rle>
  [1]     chr2 55000001-67000000      *
  [2]    chr14 54000001-63000000      *
  [3]    chr15 87000001-95000000      *
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths
>

> as.data.frame(red_deserts_0.5_5mb)
   seqnames     start       end    width strand
1      chr1  78000001  84000000  6000000      *
2      chr2  38000001  45000000  7000000      *
3      chr2  55000001  67000000 12000000      *
4      chr2 206000001 211000000  5000000      *
5      chr2 238000001 243000000  5000000      *
6      chr4  19000001  24000000  5000000      *
7      chr5 137000001 143000000  6000000      *
8      chr6 157000001 163000000  6000000      *
9      chr7  20000001  26000000  6000000      *
10     chr9   4000001  10000000  6000000      *
11    chr10 103000001 108000000  5000000      *
12    chr13  31000001  36000000  5000000      *
13    chr13  41000001  46000000  5000000      *
14    chr14  54000001  63000000  9000000      *
15    chr14 101000001 106000000  5000000      *
16    chr15  63000001  69000000  6000000      *
17    chr15  87000001  95000000  8000000      *
18    chr16         1   5000000  5000000      *
19    chr17  76000001  81000000  5000000      *
20    chr20  17000001  22000000  5000000      *
21    chr21  22000001  27000000  5000000      *

#-----------------------------------------------------------------------------------------------------------------------
# check data coverage of these regions - ie what was the distribution of callable sites in these regions?

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
# function works, whether should add this as annotation?: callable sites across all the 40kb regions within each depleted window, or just the mean coverage within the depleted window - going for the mean - to have a single point value for each depleted window
prop_5<-calc_prop_deserts_fun_1(red_deserts_0.5_5mb)


proc_means_fun<-function(p){
hold_means<-list()
for (i in (1:length(p))) {
hold_means[[i]]<-p[[i]][[2]] }
return(as.data.frame(unlist(hold_means)))}

mcols(red_deserts_0.5_8mb)$proportion<-unlist(proc_means_fun(prop_8))
mcols(red_deserts_0.5_5mb)$proportion<-unlist(proc_means_fun(prop_5))

> red_deserts_0.5_8mb
GRanges object with 3 ranges and 1 metadata column:
      seqnames            ranges strand | proportion
         <Rle>         <IRanges>  <Rle> |  <numeric>
  [1]     chr2 55000001-67000000      * |   0.506019
  [2]    chr14 54000001-63000000      * |   0.503190
  [3]    chr15 87000001-95000000      * |   0.502111
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

> mcols(red_deserts_0.5_5mb)$proportion<-unlist(proc_means_fun(prop_5))
> red_deserts_0.5_5mb
GRanges object with 21 ranges and 1 metadata column:
       seqnames              ranges strand | proportion
          <Rle>           <IRanges>  <Rle> |  <numeric>
   [1]     chr1   78000001-84000000      * |   0.515784
   [2]     chr2   38000001-45000000      * |   0.516893
   [3]     chr2   55000001-67000000      * |   0.506019
   [4]     chr2 206000001-211000000      * |   0.527189
   [5]     chr2 238000001-243000000      * |   0.564452
   ...      ...                 ...    ... .        ...
  [17]    chr15   87000001-95000000      * |   0.502111
  [18]    chr16           1-5000000      * |   0.558520
  [19]    chr17   76000001-81000000      * |   0.584985
  [20]    chr20   17000001-22000000      * |   0.509677
  [21]    chr21   22000001-27000000      * |   0.501037
  -------
  seqinfo: 22 sequences from an unspecified genome; no seqlengths

# save with the proportion
depleted_0.5_5mb_8mb<-list(red_deserts_0.5_5mb,red_deserts_0.5_8mb)
save(depleted_0.5_5mb_8mb,file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/depleted.0.5.5mb.8mb")



#-----------------------------------------------------------------------------------------------------------------------


# go through whether this has been performed correctly * 

# # Q - the reduced_deserts: removed centromeres & windows overlapping centromeres a priori
# could try without doing this step (ie assess all bins - without removing those overlapping centromeres - as a second step)
    # run this & see if differs

# MK: centromere filtering should not make such a difference

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


calc_prop_deserts_perchr(red_deserts_0.5_5mb) # for regions >=5mb
[[1]]
[1] 0.02407216

[[2]]
[1] 0.1192437

[[3]]
[1] 0

[[4]]
[1] 0.02615688

[[5]]
[1] 0.0331647

[[6]]
[1] 0.03506412

[[7]]
[1] 0.03770297

[[8]]
[1] 0

[[9]]
[1] 0.04248888

[[10]]
[1] 0.03689091

[[11]]
[1] 0

[[12]]
[1] 0

[[13]]
[1] 0.08682826

[[14]]
[1] 0.1304151

[[15]]
[1] 0.1365435

[[16]]
[1] 0.05533743

[[17]]
[1] 0.06157999

[[18]]
[1] 0

[[19]]
[1] 0

[[20]]
[1] 0.07933294

[[21]]
[1] 0.1038855

[[22]]
[1] 0

> calc_prop_deserts_perchr(red_deserts_0.5_8mb) # for regions >=8mb
[[1]]
[1] 0

[[2]]
[1] 0.04934223

[[3]]
[1] 0

[[4]]
[1] 0

[[5]]
[1] 0

[[6]]
[1] 0

[[7]]
[1] 0

[[8]]
[1] 0

[[9]]
[1] 0

[[10]]
[1] 0

[[11]]
[1] 0

[[12]]
[1] 0

[[13]]
[1] 0

[[14]]
[1] 0.08383827

[[15]]
[1] 0.07802488

[[16]]
[1] 0

[[17]]
[1] 0

[[18]]
[1] 0

[[19]]
[1] 0

[[20]]
[1] 0

[[21]]
[1] 0

[[22]]
[1] 0


