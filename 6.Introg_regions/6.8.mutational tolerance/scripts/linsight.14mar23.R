# Tue 14 Mar 2023 14:53:27 CET
# 1) assess linsight scores - mutational tolerance in noncoding regions
   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5395419/pdf/nihms852361.pdf


#-----------------------------------------------------------------------------------------------------------------------

# this section has already been run *

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1  BCFTOOLS/1.12
# to use precomputed linsight scores to hg19 using - https://bioconductor.org/packages/release/bioc/vignettes/GenomicScores/inst/doc/GenomicScores.html
#BiocManager::install("GenomicScores")

#library(GenomicScores)
#linsight <- getGScores("linsight.UCSC.hg19") # needs acess to internet

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
library(valr)
library(GenomicFeatures)

# read in linsight gscore object
load(file=paste("/scratch/devel/hpawar/admix/overlap.s*.skov/revisions_8feb23/linsight/linsight_hg19"),verbose=T)


# read in introgressed regions
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbb_99",verbose=T)
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/overlap_objects/ov_gbg_99",verbose=T)


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

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

el_lin<-proc_lin(ov_gbg_99)



#-----------------------------------------------------------------------------------------------------------------------

# population-wise means 
colMeans(data.frame(matrix(unlist(do.call(rbind,mg_lin)), ncol = ncol(do.call(rbind,mg_lin)), byrow = F)))
colMeans(data.frame(matrix(unlist(do.call(rbind,el_lin)), ncol = ncol(do.call(rbind,el_lin)), byrow = F)))


