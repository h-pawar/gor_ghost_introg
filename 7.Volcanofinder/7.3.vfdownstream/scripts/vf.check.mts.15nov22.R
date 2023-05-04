# Tue 15 Nov 2022 10:50:13 CET
# check frequency of missense mts identified in vf candidate genes

#module load  R/4.0.1  BCFTOOLS/1.14 
#-----------------------------------------------------------------------------------------------------------------------

#intersect_introgbed_regions<-function(ind,chrom){ # run through function interactively, as sanity check 

ind=2
chrom=12

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating

test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1

# cols 3:23 easterns
# 24:51 westerns

westerns<-test_1[,c(24:51)]

westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))

x<- easterns[which(rowSums(westerns)==0),]

mts_withannotn<-out_mts[[ind]][,c(1:6)]

test_1[which(rowSums(westerns)==0),]

mts_withannotn[which(mts_withannotn$start==(test_1[which(rowSums(westerns)==0),][,2]),]


pos_allw_0<-as.numeric(test_1[which(rowSums(westerns)==0),][,2])


hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==pos_allw_0[[i]]),]
}

 do.call(rbind,hold_out)

 nrow( do.call(rbind,hold_out))

nrow(test_1[which(rowSums(westerns)==0),])

test_1[which(rowSums(westerns)==0),][c(1:6,8),]


missense_tasr214<-cbind(do.call(rbind,hold_out),test_1[which(rowSums(westerns)==0),][c(1:6,8),])
# ie positions wihtin the gene, where all westerns are 0/0 & easterns are segregating & the site is functionally characterised as a missense variant

missense_tasr214[,c(9:29)]


#-----------------------------------------------------------------------------------------------------------------------

gene_1<-extract_mts_fun(1,5) # sema5a


gene_4<-extract_mts_fun(4,4) # tmprss11e
#-----------------------------------------------------------------------------------------------------------------------
#gene_1<-extract_mts_fun(1,5) # sema5a
ind=1
chrom=5

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating

test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1

# then split & take row sums

# cols 3:23 easterns
# 24:51 westerns

westerns<-test_1[,c(24:51)]
westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))


x<- easterns[which(rowSums(westerns)==0),]


test_1[which(rowSums(westerns)==0),]


pos_allw_0<-as.numeric(test_1[which(rowSums(westerns)==0),][,2])
mts_withannotn<-out_mts[[ind]][,c(1:6)]


hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==pos_allw_0[[i]]),]
}


missense_sema5a<-cbind(do.call(rbind,hold_out),test_1[which(rowSums(westerns)==0),])


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#gene_4<-extract_mts_fun(4,4) # tmprss11e
ind=4
chrom=4

test_1<-intersect_introgbed_regions(ind,chrom)

# rows where last 28 (all westerns) are 0/0 & first 21 are segregating

test_1[which(test_1=="0|0")]<-0
test_1[which(test_1=="1|1")]<-2
test_1[which(test_1=="0|1")]<-1


# cols 3:23 easterns
# 24:51 westerns

westerns<-test_1[,c(24:51)]

westerns <- as.data.frame(westerns)  %>%
              mutate(across(everything(), as.numeric))

easterns<-test_1[,c(3:23)]
easterns <- as.data.frame(easterns)  %>%
              mutate(across(everything(), as.numeric))


x<- easterns[which(rowSums(westerns)==0),]


tmprss11e_e<-test_1[which(rowSums(westerns)==0),]


pos_allw_0<-as.numeric(test_1[which(rowSums(westerns)==0),][2])
mts_withannotn<-out_mts[[ind]][,c(1:6)]


hold_out<-list()
for (i in (1:length(which(rowSums(westerns)==0)))) {  
hold_out[[i]]<-mts_withannotn[which(mts_withannotn$start==pos_allw_0[[i]]),]
}


tmprss11e_ann<- hold_out[[1]]

table(unlist(test_1[which(rowSums(westerns)==0),])[3:23])

 colnames(tmprss11e_out)<-NULL

tmprss11e_out<-as.data.frame(c(hold_out[[1]],tmprss11e_e))

check_missense_15nov22<-list(missense_tasr214,missense_sema5a,tmprss11e_out)

save(check_missense_15nov22,file="/scratch/devel/hpawar/admix/volcanofinder/output/merged/check_missense_15nov22")
