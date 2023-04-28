#Fri  7 Oct 2022 09:38:36 CEST
# calculate how much of each chr is covered by skov fragments (prop of skov introg regions per chr, incl x chr)
# how mcuh of the chr is usable = denominator for calc 
#& calc x chr:autosome ratio

#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1 
require(data.table)
library(GenomicRanges)

# 1) read in weight files for skov hmm

# callable sites in 1kb non-overlapping windows
cov_windows<-read.table("/scratch/devel/mkuhlwilm/arch/N2_gorilla_weights_float.txt.gz")


cov_windows[,1]<-sapply("chr",paste, cov_windows[,1], sep="")
# filter eg 1/2 of each chr should be covered
cov_windows<-cov_windows[which(cov_windows$V3 > 0.5),]
informranges<-GRanges(seqnames=cov_windows[,1],ranges=IRanges(start=as.numeric(cov_windows[,2]),end=as.numeric(cov_windows[,2]+1000),names=cov_windows[,1]),strand=rep("*",length(cov_windows[,1])),prop=as.numeric(cov_windows[,3]))

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# 2) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects
# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")

# individuals for GBB & GBG
sk_ids_gbb<-unique(sk$V4)[1:12]
sk_ids_gbg<-unique(sk$V4)[13:21]

#-----------------------------------------------------------------------------------------------------------------------

# convert skov outliers to granges objects per individual
skov_proc_proportion<-function(lids) {

# split skov data -> skov fragments per individual
sk_per_id<-list()
for (ind in (1:length(lids))) {
sk_per_id[[ind]]<-subset(sk, V4 == lids[ind])
}

# convert skov per id also to granges 
skov_converttoranges_function<-function(nput) {
sk<-sk_per_id[[nput]]  
skranges<-GRanges(seqnames=sk[,1],ranges=IRanges(start=as.numeric(sk[,2]),end=as.numeric(sk[,3]),names=sk[,1]),strand=rep("*",length(sk[,1])))
return(skranges)
}

sk_allids<-list()
for (ind in (1:length(lids))) {
sk_allids[[ind]]<-skov_converttoranges_function(ind)
}

reduce_sk_allids<-list()
for (ind in (1:length(lids))) {
reduce_sk_allids[[ind]]<-reduce(sk_allids[[ind]])
}

return(reduce_sk_allids)
}

#-----------------------------------------------------------------------------------------------------------------------

sk_regions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
sk_regions_gbg<-skov_proc_proportion(sk_ids_gbg)

#-----------------------------------------------------------------------------------------------------------------------

# intersect skov regions per id with the weight file to annotate these regions with the callable proportion

# 1) calculate the outliers as they are / callable proportion of each chr -> start with this 

chrs<-unique(informranges@seqnames)
chrs<-as.character(chrs)


# calculate this per chromosome per individual

check_prop_perid<-function(scen){
hold_id<-list()
for (i in (1:length(chrs))) {
chrom<-chrs[[i]]
# outlier windows for given chr / callable proportion for given chr 
hold_id[[i]]<-sum(width(scen[scen@seqnames==chrom])-1)/sum(width(informranges[informranges@seqnames==chrom])-1)}

df_id<-as.data.frame(unlist(hold_id))
return(df_id)}


allids_gbb<-list()
for (i in (1:length(sk_regions_gbb))) {
allids_gbb[[i]]<-check_prop_perid(sk_regions_gbb[[i]])}

sk_gbb_prop<-do.call(cbind,allids_gbb)


allids_gbg<-list()
for (i in (1:length(sk_regions_gbg))) {
allids_gbg[[i]]<-check_prop_perid(sk_regions_gbg[[i]])}

sk_gbg_prop<-do.call(cbind,allids_gbg)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Tue 18 Oct 2022 11:13:09 CEST
# calculate autosome: x chr ratio per individual

# ie calc per autosome:x per
calc_autx_perid<-function(scen, id){
outc<-list()
# for the 22 autosomes / xchr
for (i in (1:22)) {
   # autosome proportion for id i / x chr proportion for id i 
outc[[i]]<-scen[[id]][i]/scen[[id]][23]}
vals_1id<-unlist(outc)
# directly output mean across autosomes (1 val per id, here)
return(mean(vals_1id))
}


calc_autx_onescen<-function(scen){
all_onescen<-list()
for (i in (1:length(scen))) {
all_onescen[[i]]<-calc_autx_perid(scen,i)}
return(all_onescen)
}

m_gbb<-calc_autx_onescen(sk_gbb_prop)
m_gbg<-calc_autx_onescen(sk_gbg_prop)


#-----------------------------------------------------------------------------------------------------------------------


# plot violinplot of A:X for all individuals. - using proportions
x<-list(unlist(m_gbb),unlist(m_gbg))

prop1_df<-cbind( as.data.frame(unlist(x)), as.data.frame(unlist(list(rep("MG", 12),rep("EL",9)))))
colnames(prop1_df)<-c("x_aut","pop")

# last entry is infinte (b/c division by 0) -> remove this row for plotting
prop1_df<-prop1_df[-21,]

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/plots/xaut_ratio_skovonly.17oct22.v3.1.pdf")


ggplot(prop1_df, aes(x=pop, y=x_aut,fill=pop)) +
geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 5.0, ymax = 8.9),
            alpha = 0.01,
            fill = "grey") + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#91C79B","#3A7E99")) + 
geom_boxplot(width = 0.07) + theme_classic()   + 
theme(legend.position="none") + scale_x_discrete(limits=c("MG", "EL")) + ylab("Autosome:X ratio") + xlab("Population")  +
 geom_hline(yintercept=7.54, color = "red") + 
  theme(axis.title.x = element_blank()) +  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15))  +  theme(axis.text.y = element_text(size = 15)) +
 theme(axis.line.y = element_line(), 
   axis.line.x = element_blank(),
   axis.ticks.x=element_blank())
dev.off()


