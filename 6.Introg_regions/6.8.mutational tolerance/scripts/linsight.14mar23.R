# Tue 14 Mar 2023 14:53:27 CET
# 1) assess linsight scores - mutational tolerance in noncoding regions
   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395419/pdf/nihms852361.pdf

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12

#-----------------------------------------------------------------------------------------------------------------------

# not all introg regions are annotated with a linsight score 
# then will need to generate random regions only of those introg fragments which have a corresponding linsight score?


# high > 0.8
# but unclear to me what threshold to take as low??
# from the paper - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395419/pdf/nihms852361.pdf
#, and unannotated intronic and intergenic regions
#exhibit the least constraint (median scores of 0.044–0.048). 

#-----------------------------------------------------------------------------------------------------------------------
# this section has already been run *

# to install & use precomputed linsight scores to hg19 using - https://bioconductor.org/packages/release/bioc/vignettes/GenomicScores/inst/doc/GenomicScores.html
#BiocManager::install("GenomicScores")

#library(GenomicScores)
#linsight <- getGScores("linsight.UCSC.hg19") # needs acess to internet

#Another way to retrieve genomic scores is by using the AnnotationHub, 
#which is a web resource that provides a central location where genomic files 
#(e.g., VCF, bed, wig) and other resources from standard (e.g., UCSC, Ensembl) and distributed sites, can be found
##  [3] "linsight.UCSC.hg19"   


#linsight
#GScores object 
# organism: Homo sapiens (UCSC, hg19)
# provider: UCSC
# provider version: 19Aug2016
# download date: Apr 10, 2018
# loaded sequences: default
# maximum abs. error: 0.005
# use 'citation()' to cite these data in publications


#populations(linsight)
#defaultPopulation(linsight)

# populations(linsight)
#[1] "default"
#> defaultPopulation(linsight)
#[1] "default"
# just the one score

#mkdir -p /scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight

# save as an R object
#save(linsight,file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight/linsight_hg19"))
#-----------------------------------------------------------------------------------------------------------------------

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
#-----------------------------------------------------------------------------------------------------------------------
library(GenomicScores)
require(data.table)
options(stringsAsFactors=F)
options(scipen=100)
library(mgcv)
library(GenomicRanges)
library(tidyr)
library(ggbio)
#library(pheatmap)
#library(ggplot2)
#library(ggpubr)
library(valr)
library(GenomicFeatures)

# read in linsight gscore object
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight/linsight_hg19"),verbose=T)


# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# as a function

lin_perind<-function(scen,id){
gr1sco <- gscores(linsight, scen[[id]])
# retain regions where there is a linsight score
    # some introg regions annotated with NA instead of a linsight score
gr1sco<- gr1sco[!is.na((mcols(gr1sco))),]
# subset this object by high (>0.8), low (<0.0448) scores, & mean val across all introg regions for whcih there is a score 

# output list of 
# 1 & 2) counts of bp with high/low linsight scores,  : raw counts
# 3) mean linsihgt score - across regions where a score is available, 
# 4) proportion of high linsight bp across total introg bp for this ind, : proportions
# 5) proportion of low linsight bp across total introg bp for this ind : proportions
out<-list( sum(width(gr1sco[ gr1sco$default  > 0.8 ])-1),
   sum(width(gr1sco[ gr1sco$default < 0.0448 ])-1),
   mean((mcols(gr1sco))[,1]),
 (sum(width(gr1sco[ gr1sco$default > 0.8 ])-1))/ (sum(width(scen[[id]])-1)),
 (sum(width(gr1sco[ gr1sco$default < 0.0448 ])-1))/ (sum(width(scen[[id]])-1))
   )

return(out)}



proc_lin<-function(scen){
hold<-list()
for (ind in (1:length(scen))) {
print(ind)
hold[[ind]]<-lin_perind(scen,ind)
}
return(hold)}


mg_lin<-proc_lin(ov_gbb_99)


 do.call(rbind,mg_lin)
      [,1] [,2]   [,3]       [,4] [,5]       
 [1,] 0    104000 0.06910786 0    0.003607854
 [2,] 0    330000 0.06935042 0    0.009531512
 [3,] 0    133000 0.07069143 0    0.003887184
 [4,] 0    105000 0.06898641 0    0.003501401
 [5,] 0    65000  0.07160991 0    0.001980258
 [6,] 0    96000  0.07251344 0    0.002920507
 [7,] 0    299000 0.06938208 0    0.009587635
 [8,] 0    419000 0.07048759 0    0.01264638 
 [9,] 0    292000 0.07171551 0    0.008423968
[10,] 0    0      0.07473761 0    0          
[11,] 0    429000 0.07005035 0    0.0159991  
[12,] 0    338000 0.07004959 0    0.01041249 
>

# ok this doesnt look v promising.. no sites with high linsight scores in the introg regions


el_lin<-proc_lin(ov_gbg_99)

 do.call(rbind,el_lin)

>  do.call(rbind,el_lin)
      [,1] [,2]   [,3]       [,4] [,5]       
 [1,] 0    134000 0.07112629 0    0.005030597
 [2,] 0    140000 0.07334326 0    0.005482886
 [3,] 0    65000  0.07127086 0    0.002182746
 [4,] 0    143000 0.07275003 0    0.005527638
 [5,] 0    102000 0.07207721 0    0.003933516
 [6,] 0    356000 0.06972328 0    0.01343802 
 [7,] 0    0      0.07236936 0    0          
 [8,] 0    179000 0.06971464 0    0.00733276 
 [9,] 0    143000 0.07277368 0    0.004175182



# from the paper :    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395419/pdf/nihms852361.pdf
# . The genome-wide average of the LINSIGHT scores is about 0.07, suggesting that about
#7% of noncoding sites are under evolutionary constraint, consistent with numerous previous
#studies3,
#-----------------------------------------------------------------------------------------------------------------------

# population-wise means 
colMeans(data.frame(matrix(unlist(do.call(rbind,mg_lin)), ncol = ncol(do.call(rbind,mg_lin)), byrow = F)))

     X1               X2               X3               X4 
     0.000000000 217500.000000000      0.070723517      0.000000000 
              X5 
     0.006874859 


colMeans(data.frame(matrix(unlist(do.call(rbind,el_lin)), ncol = ncol(do.call(rbind,el_lin)), byrow = F)))
    X1               X2               X3               X4 
     0.000000000 140222.222222222      0.071683179      0.000000000 
              X5 
     0.005233705 

