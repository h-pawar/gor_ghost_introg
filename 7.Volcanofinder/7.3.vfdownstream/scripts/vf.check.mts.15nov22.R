# Tue 15 Nov 2022 10:50:13 CET
# check frequency of missense mts identified in vf candidate genes

# RUN INTERACTIVELY - function from vf.missense.29jul22.R # to check the volcanofinder mutations identified
#module load  R/4.0.1  BCFTOOLS/1.14 
#-----------------------------------------------------------------------------------------------------------------------

#intersect_introgbed_regions<-function(ind,chrom){ # run through function interactively, as sanity check of gts in westerns in this gene

ind=2
chrom=12

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating
# contingent on the derived allele being the one under selection - ie looking at sites where all westerns are 0/0 not 2/2



#head(test_1)
#     [,1]    [,2]       [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9]  [,10] [,11]
#[1,] "chr12" "11091076" "1|1" "0|1" "0|1" "1|1" "1|1" "1|1" "0|0" "1|1" "1|1"
#[2,] "chr12" "11091164" "0|0" "0|0" "0|0" "0|0" "0|0" "0|0" "0|0" "0|0" "0|0"
#[3,] "chr12" "11091560" "1|1" "0|1" "0|1" "1|1" "1|1" "1|1" "0|0" "1|1" "1|1"



test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1

# then split & take row sums

# cols 3:23 easterns
# 24:51 westerns

westerns<-test_1[,c(24:51)]
# convert columns from character to numeric
westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))


# extract those rows where westerns rowsum == 0

#which(rowSums(westerns)==0)
#[1] 1 2

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))


#easterns[which(rowSums(westerns)==0),]
#> easterns[which(rowSums(westerns)==0),]
#  V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
#1  0  1  0  1  0  0  0  0  0   0   1   1   0   0   0   1   0   0   1   0   2
#2  0  0  0  0  0  0  0  0  0   0   0   0   0   1   0   0   1   1   0   0   0
# these are segregating 

# write as a function, to go over all genes, & to output also the variant position

x<- easterns[which(rowSums(westerns)==0),]


# ie these positions are fixed 0/0 in westerns in this gene
#> which(rowSums(westerns)==0)
#[1]  3 14 16 38 59 83 87 88
#> westerns[which(rowSums(westerns)==0),]
 #  V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
#3   0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0
#14  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0

#   V22 V23 V24 V25 V26 V27 V28
#3    0   0   0   0   0   0   0
#14   0   0   0   0   0   0   0
#16   0   0   0   0   0   0   0


# in the easterns
> x
   V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
3   2  1  1  2  2  2  0  2  2   0   1   1   2   2   2   2   2   2   2   2   2
14  1  1  1  0  1  2  0  0  2   0   1   1   2   2   2   2   2   2   2   2   2
16  1  1  1  0  1  2  0  0  2   0   1   1   2   1   2   2   2   2   2   2   2
38  2  1  1  2  2  2  1  2  1   0   1   1   2   2   2   2   2   2   2   2   1
59  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   1   0   0   0   0   0
83  0  1  1  0  0  1  0  0  1   0   1   1   1   1   1   1   1   0   1   1   0
87  0  1  0  0  0  1  0  0  1   0   1   1   1   1   1   1   1   0   1   1   1
88  2  1  1  2  2  1  1  2  1   0   1   1   1   1   1   1   1   2   1   1   1

# shoudl perhaps exclude pos 59? (ask mk)

test_1[which(rowSums(westerns)==0),]
     [,1]    [,2]       [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,] "chr12" "11091560" "2"  "1"  "1"  "2"  "2"  "2"  "0"  "2"   "2"   "0"  
[2,] "chr12" "11150032" "1"  "1"  "1"  "0"  "1"  "2"  "0"  "0"   "2"   "0"  
[3,] "chr12" "11150177" "1"  "1"  "1"  "0"  "1"  "2"  "0"  "0"   "2"   "0"  
[4,] "chr12" "11174543" "2"  "1"  "1"  "2"  "2"  "2"  "1"  "2"   "1"   "0"  
[5,] "chr12" "11175060" "0"  "0"  "0"  "0"  "0"  "0"  "0"  "0"   "0"   "0"  
[6,] "chr12" "11214095" "0"  "1"  "1"  "0"  "0"  "1"  "0"  "0"   "1"   "0"  
[7,] "chr12" "11214117" "0"  "1"  "0"  "0"  "0"  "1"  "0"  "0"   "1"   "0"  
[8,] "chr12" "11214119" "2"  "1"  "1"  "2"  "2"  "1"  "1"  "2"   "1"   "0"  


mts_withannotn<-out_mts[[ind]][,c(1:6)]

 mts_withannotn
     seqnames    start      end width strand          variant
33      chr12 11090872 11090873     2      * missense_variant
34      chr12 11090905 11090906     2      * missense_variant
36      chr12 11090911 11090912     2      * missense_variant


test_1[which(rowSums(westerns)==0),]

mts_withannotn[which(mts_withannotn$start==(test_1[which(rowSums(westerns)==0),][,2]),]


pos_allw_0<-as.numeric(test_1[which(rowSums(westerns)==0),][,2])


> mts_withannotn[which(mts_withannotn$start== 11091560),]
   seqnames    start      end width strand          variant
86    chr12 11091560 11091561     2      * missense_variant

# this is why was going sequentially

hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==pos_allw_0[[i]]),]
}

 do.call(rbind,hold_out)
     seqnames    start      end width strand          variant
86      chr12 11091560 11091561     2      * missense_variant
1595    chr12 11150032 11150033     2      * missense_variant
1606    chr12 11150177 11150178     2      * missense_variant
2482    chr12 11174543 11174544     2      * missense_variant
2529    chr12 11175060 11175061     2      * missense_variant
3649    chr12 11214095 11214096     2      * missense_variant
3654    chr12 11214119 11214120     2      * missense_variant


 nrow( do.call(rbind,hold_out))
[1] 7
nrow(test_1[which(rowSums(westerns)==0),])
[1] 8
# filter out those which are not missense mts # missing position:11214117
test_1[which(rowSums(westerns)==0),][c(1:6,8),]


missense_tasr214<-cbind(do.call(rbind,hold_out),test_1[which(rowSums(westerns)==0),][c(1:6,8),])

# should output this table **
# ie positions wihtin the gene, where all westerns are 0/0 & easterns are segregating & the site is functionally characterised as a missense variant


# & the equivalent for the other genes 

> missense_tasr214
     seqnames    start      end width strand          variant     1        2 3
86      chr12 11091560 11091561     2      * missense_variant chr12 11091560 2
1595    chr12 11150032 11150033     2      * missense_variant chr12 11150032 1
1606    chr12 11150177 11150178     2      * missense_variant chr12 11150177 1
2482    chr12 11174543 11174544     2      * missense_variant chr12 11174543 2
2529    chr12 11175060 11175061     2      * missense_variant chr12 11175060 0
3649    chr12 11214095 11214096     2      * missense_variant chr12 11214095 0
3654    chr12 11214119 11214120     2      * missense_variant chr12 11214119 2
     4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
86   1 1 2 2 2 0  2  2  0  1  1  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0
1595 1 1 0 1 2 0  0  2  0  1  1  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0
1606 1 1 0 1 2 0  0  2  0  1  1  2  1  2  2  2  2  2  2  2  0  0  0  0  0  0  0
2482 1 1 2 2 2 1  2  1  0  1  1  2  2  2  2  2  2  2  2  1  0  0  0  0  0  0  0
2529 0 0 0 0 0 0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
3649 1 1 0 0 1 0  0  1  0  1  1  1  1  1  1  1  0  1  1  0  0  0  0  0  0  0  0
3654 1 1 2 2 1 1  2  1  0  1  1  1  1  1  1  1  2  1  1  1  0  0  0  0  0  0  0
     31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51
86    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
1595  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
1606  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
2482  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
2529  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
3649  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
3654  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0


 missense_tasr214[,c(9:29)]


#-----------------------------------------------------------------------------------------------------------------------

# also run through for sema5a (chr 5) & tmprss11e (chr 4) 




gene_1<-extract_mts_fun(1,5) # sema5a


gene_4<-extract_mts_fun(4,4) # tmprss11e
#-----------------------------------------------------------------------------------------------------------------------
#gene_1<-extract_mts_fun(1,5) # sema5a
ind=1
chrom=5

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating
# contingent on the derived allele being the one under selection - ie looking at sites where all westerns are 0/0 not 2/2

test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1

# then split & take row sums

# cols 3:23 easterns
# 24:51 westerns

westerns<-test_1[,c(24:51)]
# convert columns from character to numeric
westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))


x<- easterns[which(rowSums(westerns)==0),]


test_1[which(rowSums(westerns)==0),]

test_1[which(rowSums(westerns)==0),]
     [,1]   [,2]      [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
[1,] "chr5" "9054192" "0"  "1"  "0"  "1"  "0"  "0"  "0"  "0"   "0"   "0"  
[2,] "chr5" "9108295" "0"  "0"  "0"  "0"  "0"  "0"  "0"  "0"   "0"   "0"  



pos_allw_0<-as.numeric(test_1[which(rowSums(westerns)==0),][,2])
mts_withannotn<-out_mts[[ind]][,c(1:6)]


hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==pos_allw_0[[i]]),]
}

 do.call(rbind,hold_out)
     seqnames   start     end width strand                              variant
945      chr5 9054192 9054193     2      * splice_region_variant,intron_variant
3146     chr5 9108295 9108296     2      *                     missense_variant



missense_sema5a<-cbind(do.call(rbind,hold_out),test_1[which(rowSums(westerns)==0),])

> missense_sema5a
     seqnames   start     end width strand                              variant
945      chr5 9054192 9054193     2      * splice_region_variant,intron_variant
3146     chr5 9108295 9108296     2      *                     missense_variant
        1       2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
945  chr5 9054192 0 1 0 1 0 0 0  0  0  0  1  1  0  0  0  1  0  0  1  0  2  0  0
3146 chr5 9108295 0 0 0 0 0 0 0  0  0  0  0  0  0  1  0  0  1  1  0  0  0  0  0
     26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
945   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
3146  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
     51
945   0
3146  0


#> missense_sema5a[,c(9:29)]
#     3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
#945  0 1 0 1 0 0 0  0  0  0  1  1  0  0  0  1  0  0  1  0  2
#3146 0 0 0 0 0 0 0  0  0  0  0  0  0  1  0  0  1  1  0  0  0
#> ncol(missense_sema5a[,c(9:29)])
#[1] 21


#sem_num <- as.data.frame(missense_sema5a[,c(9:29)])  %>%
#              mutate(across(everything(), as.numeric))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#gene_4<-extract_mts_fun(4,4) # tmprss11e
ind=4
chrom=4

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating
# contingent on the derived allele being the one under selection - ie looking at sites where all westerns are 0/0 not 2/2

test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1

# then split & take row sums

# cols 3:23 easterns
# 24:51 westerns

westerns<-test_1[,c(24:51)]
# convert columns from character to numeric
westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))


x<- easterns[which(rowSums(westerns)==0),]


test_1[which(rowSums(westerns)==0),]
 [1] "chr4"     "69344621" "2"        "0"        "0"        "1"       
 [7] "1"        "1"        "1"        "1"        "1"        "0"       
[13] "1"        "1"        "0"        "0"        "0"        "0"       
[19] "0"        "0"        "0"        "0"        "0"        "0"       
[25] "0"        "0"        "0"        "0"        "0"        "0"       
[31] "0"        "0"        "0"        "0"        "0"        "0"       
[37] "0"        "0"        "0"        "0"        "0"        "0"       
[43] "0"        "0"        "0"        "0"        "0"        "0"       
[49] "0"        "0"        "0"   


tmprss11e_e<-test_1[which(rowSums(westerns)==0),]


pos_allw_0<-as.numeric(test_1[which(rowSums(westerns)==0),][2])
mts_withannotn<-out_mts[[ind]][,c(1:6)]


hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==pos_allw_0[[i]]),]
}

# only the one entry
> hold_out[[1]]
     seqnames    start      end width strand          variant
1102     chr4 69344621 69344622     2      * missense_variant


tmprss11e_ann<- hold_out[[1]]

unlist(test_1[which(rowSums(westerns)==0),])[3:23]
 [1] "2" "0" "0" "1" "1" "1" "1" "1" "1" "0" "1" "1" "0" "0" "0" "0" "0" "0" "0"
[20] "0" "0"


table(unlist(test_1[which(rowSums(westerns)==0),])[3:23])

 0  1  2 
12  8  1 

 colnames(tmprss11e_out)<-NULL
> tmprss11e_out
                                                                               
1 chr4 69344621 69344622 2 * missense_variant chr4 69344621 2 0 0 1 1 1 1 1 1 0
                                                                               
1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0


tmprss11e_out<-as.data.frame(c(hold_out[[1]],tmprss11e_e))

# output as R objects 
check_missense_15nov22<-list(missense_tasr214,missense_sema5a,tmprss11e_out)

save(check_missense_15nov22,file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/check_missense_15nov22")
