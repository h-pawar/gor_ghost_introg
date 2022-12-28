# Mon 27 Jun 2022 11:20:50 CEST
# follows from vf.overlap.skov.20jun22.R
# determine if/which mutations fall in the intersect of (vf >0.95 outliers with strict skov outliers) 
   # ie if any protein coding mts within the candidate genes identified 

#-----------------------------------------------------------------------------------------------------------------------

# gorilla VEP data -
# Mon 20 Jun 2022 17:31:49 CEST
# MK
#15:49 (1 hour ago)

#the gorilla variant effect prediction can be found here:
#/scratch/devel/shan/gorillas/vep/chr"$chrom".txt
#In principle, done on the same vcf file for all individuals, hg19.
#For example, one can subset to missense and LoF:

#for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
#   echo $chrom
#   egrep "missense|stop_gained|stop_lost|start_lost|splic|synonymous"
#/scratch/devel/shan/gorillas/vep/chr"$chrom".txt >
#/scratch/devel/mkuhlwilm/gori/"$chrom"_vep_subset.txt
#   done

# intersect this with the candidate genes?

# less /scratch/devel/shan/gorillas/vep/chr"$chrom".txt
## MOTIF_POS : The relative position of the variation in the aligned TFBP
## HIGH_INF_POS : A flag indicating if the variant falls in a high information position of the TFBP
## MOTIF_SCORE_CHANGE : The difference in motif score of the reference and variant sequences for the TFBP
#Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Am
#ino_acids     Codons  Existing_variation      Extra
#21_10966571_C/T 21:10966571     T       ENSG00000166157 ENST00000361285 Transcript      intron_variant  -       -       -       -       -       -       IM
#PACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=TPTE;SYMBOL_SOURCE=HGNC;HGNC_ID=12023;BIOTYPE=protein_coding;CANONICAL=YES;CCDS=CCDS13560.2;ENSP=ENSP0000
#0355208;SWISSPROT=TPTE_HUMAN;UNIPARC=UPI000016A18A;INTRON=7/23;HGVSc=ENST00000361285.4:c.173+2504N>A
#21_10966595_T/C 21:10966595     C       ENSG00000166157 ENST00000361285 Transcript      intron_variant  -       -       -       -       -       -       IM
#PACT=MODIFIER;STRAND=-1;VARIANT_CLASS=SNV;SYMBOL=TPTE;SYMBOL_SOURCE=HGNC;HGNC_ID=12023;BIOTYPE=protein_coding;CANONICAL=YES;CCDS=CCDS13560.2;ENSP=ENSP00000355208;SWISSPROT=TPTE_HUMAN;UNIPARC=UPI000016A18A;INTRON=7/23;HGVSc=ENST00000361285.4:c.173+2480N>G
#21_14595380_G/C 21:14595380     C       -       -       -       intergenic_variant      -       -       -       -       -       -       IMPACT=MODIFIER;VARIANT_CLASS=SNV
#21_14595469_C/T 21:14595469     T       -       -       -       intergenic_variant      -       -       -       -       -       -       IMPACT=MODIFIER;VARIANT_CLASS=SNV
#-----------------------------------------------------------------------------------------------------------------------

# for each candidate, read in the equivalent vep file

#-----------------------------------------------------------------------------------------------------------------------
# MG
#em_genes<-unique(unlist(hold_genes2))
# EL
#el_genes<-unique(unlist(o_gbg))
#e_cand<-list(em_genes,el_genes)

# already executed:
#save(e_cand, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.skov.intersect.genes")
#-----------------------------------------------------------------------------------------------------------------------
# module load R/4.0.1 
require(data.table)
options(stringsAsFactors=F)
library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures) # if need to read in the gtf

# read in the candidates (intersection of vf 0.95 outliers - strict skov outliers - genes of gtf)
load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.skov.intersect.genes", verbose=T)


# testing: in_genes=e_cand[[1]], i=11 # ie testing w chr 22

# check data structure of in_genes object
   # may need to split this into multiple functions?

# A) assess variant types
find_mts_fun<-function(in_genes,i) {
# 1) read in vep data for given chr of the candidate gene
x<-(data.frame(in_genes[[i]])[1])
chrom<- as.numeric(gsub('chr', '', x))
vep<-read.table(paste("/scratch/devel/shan/gorillas/vep/chr",chrom,".txt",sep=""))
# convert vep to granges format
vep$chr<-paste('chr',chrom,sep='')
# split column 1 into separate columns
test<-vep[,1]
# need to split V1 * by _ delimiter **
test<-do.call(rbind,strsplit(test,split="_"))
# convert to granges object
vep_ranges<-GRanges(seqnames=vep$chr,ranges=IRanges(start=as.numeric(test[,2]),end=as.numeric(test[,2])+1,names=vep$chr),strand=rep("*",length(vep[,1])),
   variant=(vep[,7]),impact=(vep[,14]))

# 2) overlap the vep data with the candidate genic region (intersect of vf-skov-gtf)
y<- findOverlaps(in_genes[[i]],vep_ranges)
y1<-unique(as.data.frame(y)[,2])

# extract mutations from vep data within range of overlap with the candidate genic region
out_genes<-list()
for (ind in (1:length(y1))) {
out_genes[[ind]]<-data.frame(vep_ranges[y1[[ind]]])
}

mts<-do.call(rbind,out_genes)

return(mts)
}

# try - # works running interactively
#mts<-find_mts_fun(e_cand[[1]], 11) # works

# calculate frequency of variant types of the mutations found
#count_types <- apply(mts[6], 2, table)
#mts[ which(mts$variant=='missense_variant'),]


# run both functions for this one region, then apply to rest of the candidate regions identified


# B) other aspects of metadata

# perhaps this should be a separate function -
   # write as a function *

process_meta_fun<-function(input) {
# assess other aspects of the metadata for this region
# first process the metadata
tes<-as.data.frame(input[,7])
colnames(tes)<-c("impact")
## Split up the values
Split <- strsplit(tes$impact, ";", fixed = TRUE)
## How long is each list element?
Ncol <- vapply(Split, length, 1L)
## Create an empty character matrix to store the results
M <- matrix(NA_character_, nrow = nrow(tes),
            ncol = max(Ncol), 
            dimnames = list(NULL, paste0("V", sequence(max(Ncol)))))
## Use matrix indexing to figure out where to put the results
M[cbind(rep(1:nrow(tes), Ncol), 
        sequence(Ncol))] <- unlist(Split, use.names = FALSE)
M1<-as.data.frame(M)
# not sure how useful this is..
# perhaps could subset by which biotype=protein coding?

# check what output format to generate **

# counts frequencies of entries for each column
counts_M1<-apply(M1, 2, table) # or output all here - ie the dfs of the mts, plus the counts of each type 

#return(list(M1, counts_M1)) # not optimal format

#count_impact <- apply(M1[1], 2, table)
#count_biotype <- apply(M1[7], 2, table)

#return(list(count_impact,count_biotype)) # also not optimal (if all of same biotype, does not specify which biotype)
return(counts_M1)
}


# running interactively -
#process_meta_fun(mts) # this is not optimal output format # perhaps V1 (impact). v7 (biotype), v10 (gene)
   # perhaps assess the impact & biotype categories?


#apply(M1[7], 2, table) # only gives a count, not the category - go back to outputting all
#V7 
#58 

# these would be the counts of interest
#> process_meta_fun(mts) [[1]]
#     IMPACT=LOW IMPACT=MODERATE IMPACT=MODIFIER 
#             11               1              46 
#> process_meta_fun(mts) [[7]]
#BIOTYPE=protein_coding 
#                    58 

# scen=1

mts_per_cand_fun<-function(scen, i){
mts<-find_mts_fun(scen,i)
# calculate frequency of variant types of the mutations found
count_types <- apply(mts[6], 2, table)
mis_mt<-mts[ which(mts$variant=='missense_variant'),]
imp_mt<-process_meta_fun(mts) [[1]]
biot_mt<-process_meta_fun(mts) [[7]]
gene_mt<-process_meta_fun(mts) [[10]]
return(list(mts,count_types,mis_mt,imp_mt,biot_mt,gene_mt))
}

# testing
#mts_per_cand_fun(e_cand[[1]], 11) # works

#$V10

#ENSP=ENSP00000215743 
# could also be good to output process_meta_fun(mts) [[7]] - to have the gene names as sep column

process_mts_cand<-function(scen){
out_mts<-list()
for (ind in (1:length(scen))) {
out_mts[[ind]]<-mts_per_cand_fun(scen, ind)
}
return(out_mts)
}


em.mts<-process_mts_cand(e_cand[[1]]) # see if runs


#> process_mts_cand(e_cand[[1]])
#srun: Job step aborted: Waiting up to 32 seconds for job step to finish.
#srun: Exceeded job memory limit
#srun: error: cnb1: task 0: Killed
#Session closed

# is too heavy, will need to send as a job **

# if fine also run for 
el.mts<-process_mts_cand(e_cand[[2]]) 
# send outputs to vectors & assess which genes in those of interest

e.mts<-list(em.mts,el.mts)

save(e.mts, file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.skov.intersect.genes.mts")
