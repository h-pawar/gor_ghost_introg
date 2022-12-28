# Fri 29 Jul 2022 09:24:48 CEST
# module load R/4.0.1 
#MK: Which of the missense/start/stop/splice variants are specific to eastern gorillas, within the overlapping fragments, and not fixed?
	# address this

# follows from 
#-----------------------------------------------------------------------------------------------------------------------
library(dplyr)

load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.21jul22",verbose=T)
#Loading objects:
#  em_mts
load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.meta.21jul22",verbose=T)
#Loading objects:
#  typ


# str(em_mts)
#List of 7
# $ :'data.frame':	21009 obs. of  7 variables:
#  ..$ seqnames: Factor w/ 1 level "chr5": 1 1 1 1 1 1 1 1 1 1 ...
#  ..$ start   : int [1:21009] 9035137 9035146 9035147 9035164 9035172 9035179 9035186 9035187 9035242 9035244 ...
#  ..$ end     : int [1:21009] 9035138 9035147 9035148 9035165 9035173 9035180 9035187 9035188 9035243 9035245 ...
#  ..$ width   : int [1:21009] 2 2 2 2 2 2 2 2 2 2 ...
#  ..$ strand  : Factor w/ 3 levels "+","-","*": 3 3 3 3 3 3 3 3 3 3 ...
#  ..$ variant : chr [1:21009] "downstream_gene_variant" "3_prime_UTR_variant" "3_prime_UTR_variant" "3_prime_UTR_variant" ...
#  ..$ impact  : chr [1:21009] "IMPACT=MODIFIER;DISTANCE=1;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=SEMA5A;SYMBOL_SOURCE=HGNC;HGNC_ID=10736;BIOTYPE=p"| __truncated__ "IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=SEMA5A;SYMBOL_SOURCE=HGNC;HGNC_ID=10736;BIOTYPE=protein_codi"| __truncated__ "IMPACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=SEMA5A;SYMBOL_SOURCE=HGNC;HGNC_ID=10736;B


# str(typ)
#List of 7
# $ :List of 2
#  ..$ : int [1:8, 1] 380 22 1 20492 26 2 13 73
#  .. ..- attr(*, "dimnames")=List of 2
#  .. .. ..$ : chr [1:8] "3_prime_UTR_variant" "5_prime_UTR_variant" "downstream_gene_variant" "intron_variant" ...
#  .. .. ..$ : chr "variant"
#  ..$ :'data.frame':	26 obs. of  7 variables:
#  .. ..$ seqnames: Factor w/ 1 level "chr5": 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..$ start   : int [1:26] 9044606 9052019 9066543 9066620 9108295 9108335 9108353 9122793 9122802 9122807 ...
#  .. ..$ end     : int [1:26] 9044607 9052020 9066544 9066621 9108296 9108336 9108354 9122794 9122803 9122808 ...
#  .. ..$ width   : int [1:26] 2 2 2 2 2 2 2 2 2 2 ...
#  .. ..$ strand  : Factor w/ 3 levels "+","-","*": 3 3 3 3 3 3 3 3 3 3 ...
#  .. ..$ variant : chr [1:26] "missense_variant" "missense_variant" "missense_variant" "missense_variant" ...
#  .. ..$ impact  : chr [1:26] "IMPACT=MODERATE;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=SEMA5A;SYMBOL_SOURCE=HGNC;HGNC_ID=10736;BIOTY


#> typ[[1]][[1]] # frequency of each mutation type within fragments for each of the 7 candidate genes (already added to spreadsheet)
#                                       variant
#3_prime_UTR_variant                        380
#5_prime_UTR_variant                         22


# typ[[1]][[2]]
#      seqnames   start     end width strand          variant
#450       chr5 9044606 9044607     2      * missense_variant
#802       chr5 9052019 9052020     2      * missense_variant
#1515      chr5 9066543 9066544     2      * missense_variant
#1518      chr5 9066620 9066621     2      * missense_variant
#3146      chr5 9108295 9108296     2      * missense_variant

# df of missense variants per candidate gene - write these out as bed files & query the vcf 

#>  typ[[2]][[2]][,c(1:6)]
#     seqnames    start      end width strand          variant
#33      chr12 11090872 11090873     2      * missense_variant
#34      chr12 11090905 11090906     2      * missense_variant
#36      chr12 11090911 11090912     2      * missense_variant
#-----------------------------------------------------------------------------------------------------------------------

#mkdir -p /scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions

# or try assessing all types in one go, the missense, splice, start/stop etc
	# for now don't include the noncoding transcript variants

#missense_variant
#missense_variant,splice_region_variant
#splice_donor_variant
#splice_region_variant,5_prime_UTR_variant
#splice_region_variant,intron_variant
#splice_region_variant,synonymous_variant
#start_lost
#stop_gained
#stop_lost
	
#-----------------------------------------------------------------------------------------------------------------------

#load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.21jul22", verbose=T)
#Loading objects:
#  e_cand_fin

# em_mts #  contains all variant types in these 7 regions

#-----------------------------------------------------------------------------------------------------------------------

assess_types_fun<-function(input_df) {
# 1) read in df with mutations of the candidate gene
mts<-input_df
#count_types <- apply(mts[6], 2, table)
mis_mt<-mts[ which(mts$variant=='missense_variant' | mts$variant=='missense_variant,splice_region_variant' | mts$variant=='splice_donor_variant' | mts$variant=='splice_region_variant,5_prime_UTR_variant' | mts$variant=='splice_region_variant,intron_variant' | mts$variant=='splice_region_variant,synonymous_variant' | mts$variant=='start_lost' | mts$variant=='stop_gained' | mts$variant=='stop_lost'),]
return(mis_mt)
}

out_mts<-list()
for (ind in (1:length(em_mts))) {
out_mts[[ind]]<-assess_types_fun(em_mts[[ind]])
}

# convert these to bed files, then query vcf at these positions

# is already in df format
#> str(out_mts[[1]][,c(1:3)])
#'data.frame':	41 obs. of  3 variables:
# $ seqnames: Factor w/ 1 level "chr5": 1 1 1 1 1 1 1 1 1 1 ...
# $ start   : int  9044606 9052019 9054192 9062992 9063221 9066543 9066620 9108295 9108335 9108353 ...
# $ end     : int  9044607 9052020 9054193 9062993 9063222 9066544 9066621 9108296 9108336 9108354 ...
 # so can directly convert to bed file *
#-----------------------------------------------------------------------------------------------------------------------

df_tobed<-function(ind){

df<-out_mts[[ind]][,c(1:3)]
df$start<-df$start-1

a=ind
tmp<-paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions/gene.",a,".tmp.bed",sep="")

write.table(df,tmp,sep="\t",row.names=F,col.names=F,quote=F) # should write this out only once

}

for (ind in (1:length(em_mts))) {
df_tobed(ind)}

#-----------------------------------------------------------------------------------------------------------------------

# no such mutations in gene 5
#[hpawar@login2 ~]$ ls -lhtr /scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions
#total 3.0K
#-rw-r--r-- 1 hpawar devel  483 Jul 29 09:52 gene.7.tmp.bed
#-rw-r--r-- 1 hpawar devel  676 Jul 29 09:52 gene.6.tmp.bed
#-rw-r--r-- 1 hpawar devel    0 Jul 29 09:52 gene.5.tmp.bed
#-rw-r--r-- 1 hpawar devel  575 Jul 29 09:52 gene.4.tmp.bed
#-rw-r--r-- 1 hpawar devel  850 Jul 29 09:52 gene.3.tmp.bed
#-rw-r--r-- 1 hpawar devel 7.4K Jul 29 09:52 gene.2.tmp.bed
#-rw-r--r-- 1 hpawar devel  861 Jul 29 09:52 gene.1.tmp.bed

#-----------------------------------------------------------------------------------------------------------------------
intersect_introgbed_regions<-function(ind,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions/gene.",a,".tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos & gts
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}
# will need to check the chr 

# gene 1, chr 5
# gene 2, chr 12
# gene 3, chr 3
# gene 4, chr 4
# gene 5, chr 1 (no such mts)
# gene 6, chr 10
# gene 7, chr 7


# check order of individuals in the vcf
#-----------------------------------------------------------------------------------------------------------------------

# samples & their populations
# for ids
samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

# add population identifier
identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
intersect_introgbed_regions<-function(ind,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions/gene.",a,".tmp.bed",sep="")

chrom=chrom
# filter by biallelic snps only, query by the putative introg regions for id 1 & output chr, pos & gts
testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))  ## ie list merged vertically into df

return(testsnps)
}
# will need to check the chr 

# gene 1, chr 5
# gene 2, chr 12
# gene 3, chr 3
# gene 4, chr 4
# gene 5, chr 1 (no such mts)
# gene 6, chr 10
# gene 7, chr 7


extract_mts_fun<-function(ind,chrom){

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


mts_withannotn<-out_mts[[ind]][,c(1:6)]

hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==(test_1[which(rowSums(westerns)==0),2][[i]])),]
}


int_mts<-list(x, hold_out)

return(int_mts)
}


# gene 1, chr 5
# gene 2, chr 12
# gene 3, chr 3
# gene 4, chr 4
# gene 5, chr 1 (no such mts)
# gene 6, chr 10
# gene 7, chr 7

gene_1<-extract_mts_fun(1,5)
gene_2<-extract_mts_fun(2,12)
gene_3<-extract_mts_fun(3,3)
gene_4<-extract_mts_fun(4,4)

gene_6<-extract_mts_fun(6,10)
gene_7<-extract_mts_fun(7,7)



#> gene_3<-extract_mts_fun(3,3)
#Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'x' in selecting a method for function 'which': subscript out of bounds

#> gene_6<-extract_mts_fun(6,10)
#Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'x' in selecting a method for function 'which': subscript out of bounds
#> gene_7<-extract_mts_fun(7,7)
#Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'x' in selecting a method for function 'which': subscript out of bounds


# run through these 3 genes interactively (below) - all had no rows where westerns all carried 0/0
#-----------------------------------------------------------------------------------------------------------------------

# ie so these are the segregating mutations of interest in easterns
> gene_1
[[1]]
  V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
1  0  1  0  1  0  0  0  0  0   0   1   1   0   0   0   1   0   0   1   0   2
2  0  0  0  0  0  0  0  0  0   0   0   0   0   1   0   0   1   1   0   0   0

[[2]]
[[2]][[1]]
    seqnames   start     end width strand                              variant
945     chr5 9054192 9054193     2      * splice_region_variant,intron_variant

[[2]][[2]]
     seqnames   start     end width strand          variant
3146     chr5 9108295 9108296     2      * missense_variant


# gene_1 corresponds to # SEMA5A
> e_cand_fin[[1]]
GRanges object with 1 range and 1 metadata column:
                  seqnames          ranges strand |         gene_id
                     <Rle>       <IRanges>  <Rle> |     <character>
  ENSG00000112902     chr5 9035138-9546187      - | ENSG00000112902
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths



> gene_2
[[1]]
   V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
3   2  1  1  2  2  2  0  2  2   0   1   1   2   2   2   2   2   2   2   2   2
14  1  1  1  0  1  2  0  0  2   0   1   1   2   2   2   2   2   2   2   2   2
16  1  1  1  0  1  2  0  0  2   0   1   1   2   1   2   2   2   2   2   2   2
38  2  1  1  2  2  2  1  2  1   0   1   1   2   2   2   2   2   2   2   2   1
59  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   1   0   0   0   0   0
83  0  1  1  0  0  1  0  0  1   0   1   1   1   1   1   1   1   0   1   1   0
87  0  1  0  0  0  1  0  0  1   0   1   1   1   1   1   1   1   0   1   1   1
88  2  1  1  2  2  1  1  2  1   0   1   1   1   1   1   1   1   2   1   1   1

[[2]]
[[2]][[1]]
   seqnames    start      end width strand          variant
86    chr12 11091560 11091561     2      * missense_variant

[[2]][[2]]
     seqnames    start      end width strand          variant
1595    chr12 11150032 11150033     2      * missense_variant

[[2]][[3]]
     seqnames    start      end width strand          variant
1606    chr12 11150177 11150178     2      * missense_variant

[[2]][[4]]
     seqnames    start      end width strand          variant
2482    chr12 11174543 11174544     2      * missense_variant

[[2]][[5]]
     seqnames    start      end width strand          variant
2529    chr12 11175060 11175061     2      * missense_variant

[[2]][[6]]
     seqnames    start      end width strand          variant
3649    chr12 11214095 11214096     2      * missense_variant

[[2]][[7]]
[1] seqnames start    end      width    strand   variant # check what happened here? - perhaps compare to start pos output by test_1
<0 rows> (or 0-length row.names)

[[2]][[8]]
     seqnames    start      end width strand          variant
3654    chr12 11214119 11214120     2      * missense_variant

# gene_2 corresponds to gene: # TAS2R14

> e_cand_fin[[2]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000212127    chr12 11090005-11324172      - | ENSG00000212127
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths



#-----------------------------------------------------------------------------------------------------------------------
# for [[2]][[7]]: whcih mutation this corresponds to..
which(mts_withannotn$start==11214117,)
# gives integer(0) # which doesnt make sense

> which(mts_withannotn$start==11214116,)
[1] 246

3652    chr12 11214116 11214117     2      * missense_variant
3654    chr12 11214119 11214120     2      * missense_variant
# thse are the possible mts it coudl correspond to - regenerate the tmp bed file && check --
    # or maybe i am doing something wrong? - have regenerated & checked, but can't see where am going wrong here..


#-----------------------------------------------------------------------------------------------------------------------

> gene_4
[[1]]
  V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
8  2  0  0  1  1  1  1  1  1   0   1   1   0   0   0   0   0   0   0   0   0

[[2]]
[[2]][[1]]
     seqnames    start      end width strand          variant
1102     chr4 69344621 69344622     2      * missense_variant


# gene_4 corresponds to gene: # tmprss11e

 e_cand_fin[[4]]
GRanges object with 1 range and 1 metadata column:
                  seqnames            ranges strand |         gene_id
                     <Rle>         <IRanges>  <Rle> |     <character>
  ENSG00000087128     chr4 69313167-69363322      + | ENSG00000087128
  -------
  seqinfo: 265 sequences from an unspecified genome; no seqlengths



missense_mts<-list(gene_1,gene_2,gene_4)
save(missense_mts,file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.meta.21jul22")


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# exploring the genes which did not output - 

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
test_1<-intersect_introgbed_regions(3,3)

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


# extract those rows where westerns rowsum == 0
which(rowSums(westerns)==0)
#integer(0)
#-----------------------------------------------------------------------------------------------------------------------

test_1<-intersect_introgbed_regions(6,10)

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


# extract those rows where westerns rowsum == 0
which(rowSums(westerns)==0)
integer(0)
#-----------------------------------------------------------------------------------------------------------------------
test_1<-intersect_introgbed_regions(7,7)

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


# extract those rows where westerns rowsum == 0
which(rowSums(westerns)==0)
integer(0)
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
