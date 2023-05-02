# Fri 29 Jul 2022 09:24:48 CEST
# module load R/4.0.1 
# Which of the missense/start/stop/splice variants are specific to eastern gorillas, within the overlapping fragments, and not fixed?

#-----------------------------------------------------------------------------------------------------------------------
library(dplyr)

load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.21jul22",verbose=T)
#Loading objects:
#  em_mts
load(file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.meta.21jul22",verbose=T)
#Loading objects:
#  typ

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

assess_types_fun<-function(input_df) {
mts<-input_df
mis_mt<-mts[ which(mts$variant=='missense_variant' | mts$variant=='missense_variant,splice_region_variant' | mts$variant=='splice_donor_variant' | mts$variant=='splice_region_variant,5_prime_UTR_variant' | mts$variant=='splice_region_variant,intron_variant' | mts$variant=='splice_region_variant,synonymous_variant' | mts$variant=='start_lost' | mts$variant=='stop_gained' | mts$variant=='stop_lost'),]
return(mis_mt)
}

out_mts<-list()
for (ind in (1:length(em_mts))) {
out_mts[[ind]]<-assess_types_fun(em_mts[[ind]])
}

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
#-----------------------------------------------------------------------------------------------------------------------
intersect_introgbed_regions<-function(ind,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions/gene.",a,".tmp.bed",sep="")

chrom=chrom

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" "))

return(testsnps)
}
#-----------------------------------------------------------------------------------------------------------------------

samples.in.vcf<-system(paste("bcftools query -l /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_22.vcf.gz",sep=""),intern=T) # -l, --list-samples: list sample names and exit

identifiers<-cbind(samples.in.vcf, c(rep("GBB",12), rep("GBG",9), rep("GGD",1),rep("GGG",27)))
colnames(identifiers)<-c("ids","pop")
#-----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
intersect_introgbed_regions<-function(ind,chrom){
a=ind
tmp<-paste("/scratch/devel/hpawar/admix/volcanofinder/output/merged/mt_regions/gene.",a,".tmp.bed",sep="")

chrom=chrom

testsnps<-system(paste("bcftools view -m2 -M2 -v snps -R ",tmp," /scratch/devel/mkuhlwilm/arch/subsets/3_gorilla_",chrom,".vcf.gz | bcftools query -f '%CHROM %POS [%GT ]\n' ",sep=""),intern=T) 
testsnps<-do.call(rbind,strsplit(testsnps,split=" ")) 

return(testsnps)
}


extract_mts_fun<-function(ind,chrom){

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating

test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1


westerns<-test_1[,c(24:51)]
westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))


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



# run through these 3 genes interactively (below) - all had no rows where westerns all carried 0/0 -> all gave  #which(rowSums(westerns)==0) = #integer(0) when running through interactively
#-----------------------------------------------------------------------------------------------------------------------

# ie so these are the segregating mutations of interest in easterns


missense_mts<-list(gene_1,gene_2,gene_4)
save(missense_mts,file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/em.el.vf.s*.skov.intersect.genes.mts.meta.21jul22")


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

