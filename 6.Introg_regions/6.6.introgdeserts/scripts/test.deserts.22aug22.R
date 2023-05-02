# Mon 22 Aug 2022 18:47:30 BST
# recalculate deserts  - excluding short windows at end of chr & excluding centromeres
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

hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22)),pruning.mode="coarse")
centromeres <- hg19[hg19$gieStain=='acen']

#-----------------------------------------------------------------------------------------------------------------------

# non-overlapping windows across chromosomes
mygenome <- readDNAStringSet("/home/devel/marcmont/scratch/snpCalling_hg19/chimp/assembly/BWA/hg19.fa")
chrSizes <- width(mygenome)
names(chrSizes) <- names(mygenome)

#-----------------------------------------------------------------------------------------------------------------------

bins   <- tileGenome(chrSizes[1:22], tilewidth=1000000, cut.last.tile.in.chrom=T)

#-----------------------------------------------------------------------------------------------------------------------
  
# only retain 1mb regions
bins<-bins[width(bins)==1000000]

#-----------------------------------------------------------------------------------------------------------------------

# subtract centromeres granges from bins
bins<- bins[!bins %over% centromeres]
# & again only retain 1mb regions
bins<-bins[width(bins)==1000000]

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
hits<-findOverlaps(bins, scen[[ind]], ignore.strand=TRUE)
x1<-unique(as.data.frame(hits)[,1])
test<-bins[x1]
tmp<-(as.data.frame(hits)[,1])
r_freq<-table(tmp)
mcols(test)$id<-r_freq
empty<-bins[!(bins %over% test)]
mcols(empty)$id<-0
grlist <- GenomicRanges::GRangesList(test,empty)
tmp = unlist(grlist)
tmp <- sortSeqlevels(tmp)
tmp <- sort(tmp)
return(tmp)
}


count_across_ids_fun<-function(scen){
hold_counts<-list()
for (i in (1:length(scen))) {
hold_counts[[i]]<-count_regions_fun(scen,i)}
return(hold_counts)
}


f_gbb<-count_across_ids_fun(ov_gbb_99)
f_gbg<-count_across_ids_fun(ov_gbg_99)

process_fun<-function(scen){
hold_id<-list()
for (i in (1:length(scen))) {  
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
fG<-GRanges(seqnames=deserts[,1],ranges=IRanges(start=as.numeric(deserts[,2]),end=as.numeric(deserts[,3]),names=deserts[,2]),strand=rep("*",length(deserts[,1])))
reduce(fG)
sum(width(reduce(fG)))/sum(width(bins)) # deserts / bp of bins

sum(width(reduce(fG)))/ sum(chrSizes[1:22]) # deserts / size of whole genome

sum(width(reduce(fG)))/1000000

sum(width(bins))/1000000 # size of the bins covering the genome (excluding centromeres & windows < 1mb)


#-----------------------------------------------------------------------------------------------------------------------

t<-reduce(fG) # unique deserts

reduced_deserts<-reduce(fG) 

save(reduced_deserts,file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/reduced.deserts.23aug22")


calc_prop_deserts_fun<-function(i){
x<-paste("chr",i,sep="")
sum(width(t[t@seqnames==x]))/chrSizes[i]
}


hold_prop<-list()
for (i in (1:22)) {
hold_prop[[i]]<-calc_prop_deserts_fun(i) }


#-----------------------------------------------------------------------------------------------------------------------

G_e1<-GRanges(seqnames=g_e1[,1],ranges=IRanges(start=as.numeric(g_e1[,2]),end=as.numeric(g_e1[,3]),names=g_e1[,2]),strand=rep("*",length(g_e1[,1])))
mcols(G_e1)$popfreq<-rowSums(g_e)

#-----------------------------------------------------------------------------------------------------------------------

nondeserts<-g_e1[rowSums(g_e)!=0,]

G_nondeserts<-GRanges(seqnames=nondeserts[,1],ranges=IRanges(start=as.numeric(nondeserts[,2]),end=as.numeric(nondeserts[,3]),names=nondeserts[,2]),strand=rep("*",length(nondeserts[,1])))

mcols(G_nondeserts)$popfreq<-rowSums(nondeserts[,-c(1:4)])


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# Tue 23 Aug 2022 15:29:10 BST


# assess width of deserts
load(file="/scratch/devel/hpawar/admix/overlap.s*.skov/deserts/reduced.deserts.23aug22",verbose=T)

# frequency table of deserts inferred
table(width(reduced_deserts))


# filter by those >=5mb as the depleted regions
reduced_deserts1<-reduced_deserts[width(reduced_deserts)>=5000000]

sum(width(reduced_deserts1))/1000000
sum(width(reduced_deserts1))/ sum(chrSizes[1:22]) # depleted regions (>=5mb) / size of whole genome


reduced_deserts2<-reduced_deserts[width(reduced_deserts)>=8000000]
sum(width(reduced_deserts2))/ sum(chrSizes[1:22])

sum(width(reduced_deserts2))/1000000



calc_prop_deserts_fun<-function(i,t){
x<-paste("chr",i,sep="")
sum(width(t[t@seqnames==x]))/chrSizes[i]
}

calc_prop_deserts_perchr<-function(t){
hold_prop<-list()
for (i in (1:22)) {
hold_prop[[i]]<-calc_prop_deserts_fun(i,t) }
return(hold_prop)}


calc_prop_deserts_perchr(reduced_deserts1) # for regions >=5mb
calc_prop_deserts_perchr(reduced_deserts2) # for regions >=8mb
