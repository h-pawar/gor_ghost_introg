#Fri  7 Oct 2022 09:38:36 CEST
#skov fragments how much of each chr is covered  -> for me to calculate (prop of skov introg regions per chr, incl x chr)
#id coverage for each id for each chr -> how mcuh of the chr is usable = denominator for calc 
#& filter eg 1/2 of each piece should be covered
#& calc x chr:autosome ratio


 # outlier windows / total windows


# Sat 15 Oct 2022 15:31:03 CEST
# coverage for skov hmm windows
#-----------------------------------------------------------------------------------------------------------------------

# MK
#the weight files for the SkovHMM are here:

#/scratch/devel/mkuhlwilm/arch/N2_gorilla_weights_float.txt.gz

#They are filtered for RepeatMask and Mapability35 - i.e. any base for
#which any individual has genotypes is counted. Im not sure how the file
#looks like because Skov changed the format at some point, it either
#contains some sort of position and the weight, or just the weight.

#Please tell if you cannot interpret it, then we may reconstruct from
#other files.

#zcat /scratch/devel/mkuhlwilm/arch/N2_gorilla_weights_float.txt.gz | head
#1  0  0.000
#1  1000  0.000
#1  2000  0.000
#1  3000  0.000

#zcat /scratch/devel/mkuhlwilm/arch/N2_gorilla_weights_float.txt.gz | tail
#X  155265000   0.000
#X  155266000   0.000


# MK
#these are non-overlapping windows - if I remember correctly in 1kbp steps 
#(so column 2 is the start coordinate, but please check the tail of the file).


#-----------------------------------------------------------------------------------------------------------------------
## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1 
require(data.table)
library(GenomicRanges)

# 1) read in weight files for skov hmm

# callable sites in 1kb non-overlapping windows
cov_windows<-read.table("/scratch/devel/mkuhlwilm/arch/N2_gorilla_weights_float.txt.gz")

#> head(cov_windows)
#  V1   V2 V3
#1  1    0  0
#2  1 1000  0

# to match the granges object will need to add 'chr'

cov_windows[,1]<-sapply("chr",paste, cov_windows[,1], sep="")
# apply cutoff to the proportion
   # whether to implement this here? or only later?
cov_windows<-cov_windows[which(cov_windows$V3 > 0.5),]
informranges<-GRanges(seqnames=cov_windows[,1],ranges=IRanges(start=as.numeric(cov_windows[,2]),end=as.numeric(cov_windows[,2]+1000),names=cov_windows[,1]),strand=rep("*",length(cov_windows[,1])),prop=as.numeric(cov_windows[,3]))

# pre-filter by proportion
# informranges
#GRanges object with 3036269 ranges and 1 metadata column:
#       seqnames              ranges strand |      prop
#          <Rle>           <IRanges>  <Rle> | <numeric>
#  chr1     chr1              0-1000      * |         0
#  chr1     chr1           1000-2000      * |         0
#  chr1     chr1           2000-3000      * |         0
#  chr1     chr1           3000-4000      * |         0
#  chr1     chr1           4000-5000      * |         0
#   ...      ...                 ...    ... .       ...
#  chrX     chrX 155264000-155265000      * |         0
#  chrX     chrX 155265000-155266000      * |         0
#  chrX     chrX 155266000-155267000      * |         0
#  chrX     chrX 155267000-155268000      * |         0
#  chrX     chrX 155268000-155269000      * |         0
#  -------
#  seqinfo: 23 sequences from an unspecified genome; no seqlengths

# post-filter by proportion > 0.5
#> informranges
#GRanges object with 1414961 ranges and 1 metadata column:
#       seqnames              ranges strand |      prop
#          <Rle>           <IRanges>  <Rle> | <numeric>
#  chr1     chr1       755000-756000      * |     0.569
#  chr1     chr1       756000-757000      * |     0.923
#  chr1     chr1       757000-758000      * |     0.917
#  chr1     chr1       758000-759000      * |     0.803
#  chr1     chr1       761000-762000      * |     0.547
#   ...      ...                 ...    ... .       ...
#  chrX     chrX 155230000-155231000      * |     0.755
#  chrX     chrX 155231000-155232000      * |     0.705
#  chrX     chrX 155232000-155233000      * |     0.676
#  chrX     chrX 155233000-155234000      * |     0.897
#  chrX     chrX 155234000-155235000      * |     0.669
#  -------
#  seqinfo: 23 sequences from an unspecified genome; no seqlengths
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# 2) read in skov strict outliers & using function below split per inidivdual & generate reduced granges objects
# skov HMM bed file with the putative introgressed fragments ( ~2% with strict cutoff)
sk<-read.table("/scratch/devel/mkuhlwilm/arch/skov/gor_fragments_strict.bed")
#head(sk)
#    V1       V2       V3                                V4
#1 chr1  5199000  5219000 Gorilla_beringei_beringei-Bwiruka
#2 chr1  5402000  5460000 Gorilla_beringei_beringei-Bwiruka

#> nrow(sk)
#[1] 15691

# convert startpos & endpos from int to numeric values
#sk[,c(2:3)]<-as.data.frame(lapply(sk[,c(2:3)], as.numeric)) # not necessary - b/c doing this step in the convert to granges function below

#head(sk)
#    V1       V2       V3                                V4
#1 chr1  5199000  5219000 Gorilla_beringei_beringei-Bwiruka
#2 chr1  5402000  5460000 Gorilla_beringei_beringei-Bwiruka

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
# perhaps already reduced? in the sk_allids step?
return(reduce_sk_allids)
}

#-----------------------------------------------------------------------------------------------------------------------

sk_regions_gbb<-skov_proc_proportion(sk_ids_gbb) # works
sk_regions_gbg<-skov_proc_proportion(sk_ids_gbg)

#-----------------------------------------------------------------------------------------------------------------------

# intersect skov regions per id with the weight file to annotate these regions with the callable proportion

# 1) calculate the outliers as they are / callable proportion of each chr -> start with this 

#scen = sk_regions_gbb[[1]]
#sum(width(scen[scen@seqnames=="chr1"])-1)/sum(width(informranges[informranges@seqnames=="chr1"])-1)
#[1] 0.0287833


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

# Q - whether shoudl filter the outlier regions by proportion covered/callable sites?

#-----------------------------------------------------------------------------------------------------------------------

# 2) x chr:autosome ratio per individual


#> sum(sk_gbb_prop[c(1:22),1])
#[1] 0.9254712
#> sk_gbb_prop[23,1]/sum(sk_gbb_prop[c(1:22),1])
#[1] 0.00791827

prop_xaut_perid<-function(scen){
out<-list()
for (i in (1:length(scen))) {
   # x chr proportion for id i / autosomes proportion for id i
out[[i]]<-scen[23,i]/sum(scen[c(1:22),i])}
return(out)}

# perhaps instead of proportions should be counting the windows? ie sum of x chr windows / sum of autosomal windows per id


prop_xaut_gbb<-prop_xaut_perid(sk_gbb_prop)
prop_xaut_gbg<-prop_xaut_perid(sk_gbg_prop)


> unlist(prop_xaut_gbb)
 [1] 0.007918270 0.003589299 0.003886604 0.005378198 0.004956640 0.001997716
 [7] 0.009110938 0.009021666 0.006407656 0.005385476 0.004534494 0.004277122
> unlist(prop_xaut_gbg)
[1] 0.001457314 0.005379361 0.013243692 0.004117436 0.011673187 0.008382243
[7] 0.001442747 0.009979615 0.000000000

#-----------------------------------------------------------------------------------------------------------------------

# this may be a more correct way of calculating
win_xaut_perid<-function(scen){
rat<-list()
for (i in (1:length(scen))) {
# sum of x chr outlier windows / sum of autosomal outlier windows for given individual
rat[[i]]<-sum(width(scen[[i]][scen[[i]]@seqnames=="chrX"])-1)/sum(width(scen[[i]][scen[[i]]@seqnames!="chrX"])-1)}
return(unlist(rat))}

win_xaut_gbb<-win_xaut_perid(sk_regions_gbb)
win_xaut_gbg<-win_xaut_perid(sk_regions_gbg)


> win_xaut_gbb
 [1] 0.007671646 0.003491524 0.003784540 0.005064092 0.004592661 0.001947010
 [7] 0.008588692 0.008290171 0.006620677 0.005372252 0.004584726 0.004153585
> win_xaut_gbg
[1] 0.001354463 0.005124672 0.011862057 0.003761018 0.011244087 0.008200899
[7] 0.001356637 0.009284945 0.000000000
#-----------------------------------------------------------------------------------------------------------------------

# plot of iolinplot of X:A for all individuals. - using proportions
x<-list(prop_xaut_gbb,prop_xaut_gbg)

# as.data.frame(unlist(x))
#     unlist(x)
#1  0.007918270
#2  0.003589299
#3  0.003886604


prop_df<-cbind( as.data.frame(unlist(x)), as.data.frame(unlist(list(rep("MG", 12),rep("EL",9)))))
colnames(prop_df)<-c("x_aut","pop")


# plot as violins 
library(ggplot2)
o<-ggplot(prop_df, aes(x=pop, y=x_aut,,fill=id)) + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#3A7E99","#91C79B")) + 
geom_boxplot(width = 0.07) + theme_classic() 

pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/plots/xaut_ratio_skovonly.17oct22.pdf")
ggplot(prop_df, aes(x=pop, y=x_aut,,fill=pop)) + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#3A7E99","#91C79B")) + 
geom_boxplot(width = 0.07) + theme_classic()  


dev.off()
scp -r hpawar@172.16.10.21:/scratch/devel/hpawar/admix/'overlap.s*.skov'/plots/xaut_ratio_skovonly.17oct22.pdf  /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov

# colours are the wrong way around **
# & should perhaps also plots violins for the proportion per individual **
   # try this when the lcuster is back up

#-----------------------------------------------------------------------------------------------------------------------

# replotted -

pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/plots/xaut_ratio_skovonly.17oct22.pdf")

ggplot(prop_df, aes(x=pop, y=x_aut,,fill=pop)) + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#91C79B","#3A7E99")) + 
geom_boxplot(width = 0.07) + theme_classic()   + 
theme(legend.position="none") + scale_x_discrete(limits=c("MG", "EL")) + ylab("X:Autosome ratio") + xlab("Population")
dev.off()
#-----------------------------------------------------------------------------------------------------------------------

# Tue 18 Oct 2022 11:13:09 CEST
# MK
#could you calculate the other way around (Aut:X)? 
#And add a grey bar behind that [5.0,8.9] for humans, as well as a red line at 7.54 for bonobos?

#-----------------------------------------------------------------------------------------------------------------------

q()
#-----------------------------------------------------------------------------------------------------------------------

# results from outliers / callable proportion of each chr 
sk_gbg_prop
   unlist(hold_id) unlist(hold_id) unlist(hold_id) unlist(hold_id)
1      0.038607102     0.035791246      0.05311193     0.030883871
2      0.042929774     0.043636797      0.05372577     0.049944391
3      0.038642938     0.037289798      0.04823258     0.050811394
4      0.045490221     0.037118429      0.04357517     0.043177513
5      0.036738675     0.048115389      0.04674536     0.026786791
6      0.051353297     0.043028776      0.04651608     0.039428983
7      0.077721895     0.051688643      0.06749362     0.059102517
8      0.076179470     0.064106050      0.06020394     0.062855201
9      0.054586512     0.048363294      0.05912057     0.039117370
10     0.061408781     0.054020377      0.07186103     0.060499890
11     0.056475539     0.046515500      0.03941632     0.044305522
12     0.044981892     0.021511727      0.03551611     0.040303402
13     0.045270724     0.043354097      0.07586009     0.045462386
14     0.027570362     0.034131578      0.05266646     0.027216896
15     0.032709345     0.025727013      0.02948096     0.033234897
16     0.060356390     0.067829086      0.08020146     0.057290669
17     0.031248362     0.038142925      0.02309548     0.044618046
18     0.065152392     0.047333140      0.05876254     0.054374035
19     0.051492573     0.033980352      0.05206208     0.029993830
20     0.042331433     0.045569535      0.03619253     0.041521908
21     0.032024169     0.023564955      0.04949190     0.059489151
22     0.072019968     0.105564349      0.01802021     0.028430537
23     0.001581611     0.005359904      0.01458597     0.003989175
   unlist(hold_id) unlist(hold_id) unlist(hold_id) unlist(hold_id)
1       0.03270376     0.040073882     0.033192691     0.044899770
2       0.03524786     0.035263743     0.034802987     0.035430569
3       0.03403442     0.033553954     0.040790312     0.039643085
4       0.03679402     0.059000199     0.036030097     0.044496070
5       0.03909512     0.038415590     0.036344107     0.033308125
6       0.04850721     0.053929399     0.035919185     0.045278646
7       0.03691301     0.032209308     0.048131547     0.050711419
8       0.05866757     0.070958532     0.082991162     0.063657376
9       0.04350918     0.053750822     0.030547110     0.049376789
10      0.05295023     0.076830609     0.056835007     0.078487136
11      0.03122729     0.046197627     0.045410511     0.045228869
12      0.03194118     0.041158276     0.022475403     0.039495158
13      0.07555343     0.064647820     0.063402012     0.056712985
14      0.03574427     0.033999028     0.035191976     0.019020899
15      0.02590220     0.038265178     0.021973072     0.030907453
16      0.06443489     0.045000411     0.062573564     0.069252457
17      0.04052850     0.025743197     0.026005348     0.027892833
18      0.04822531     0.058955440     0.046513310     0.043885031
19      0.04456362     0.034217645     0.051824783     0.032319301
20      0.04334334     0.054811617     0.041758019     0.036226262
21      0.03674815     0.065586377     0.025102994     0.081406207
22      0.06234019     0.110678193     0.060087666     0.034335809
23      0.01119429     0.009331506     0.001353156     0.009999297
   unlist(hold_id)
1       0.06232005
2       0.06278201
3       0.05361573
4       0.06256867
5       0.04672344
6       0.06993723
7       0.08634753
8       0.08040789
9       0.00341388
10      0.08641794
11      0.07833313
12      0.05421453
13      0.05997125
14      0.05558256
15      0.06306622
16      0.10710864
17      0.09280134
18      0.08830054
19      0.05400788
20      0.08176207
21      0.13710519
22      0.05728723
23      0.00000000


>  sk_gbb_prop
   unlist(hold_id) unlist(hold_id) unlist(hold_id) unlist(hold_id)
1      0.028783297     0.034750014     0.040128207     0.025080130
2      0.034389895     0.043279314     0.051048618     0.037646965
3      0.031671324     0.037387851     0.038603716     0.055527774
4      0.044548394     0.046055317     0.061333836     0.050952814
5      0.044969805     0.057902870     0.048203071     0.061158057
6      0.035469211     0.045807366     0.028719599     0.036897878
7      0.043610257     0.044665659     0.043766613     0.027036535
8      0.048484024     0.051638341     0.067763426     0.034575119
9      0.040913213     0.040006401     0.042495688     0.042549030
10     0.044814190     0.071303965     0.063299861     0.054577439
11     0.039900702     0.050662993     0.034542262     0.045910027
12     0.032314220     0.041919891     0.051991855     0.029298848
13     0.055563009     0.054336368     0.061734547     0.046631529
14     0.038417355     0.045729687     0.040913710     0.029978350
15     0.033034686     0.031758346     0.031558136     0.025727013
16     0.046752251     0.093531875     0.060164783     0.083376673
17     0.046846327     0.045168563     0.038457505     0.029518167
18     0.041618441     0.051070602     0.048249421     0.029851466
19     0.034360021     0.039390632     0.040861848     0.040766931
20     0.036563565     0.046716363     0.060410834     0.046648902
21     0.067838506     0.035100247     0.021367756     0.047679209
22     0.054608547     0.108121271     0.150249604     0.059661512
23     0.007328132     0.004006748     0.004375791     0.005061156
   unlist(hold_id) unlist(hold_id) unlist(hold_id) unlist(hold_id)
1      0.031952266     0.038914945     0.026755156      0.03867048
2      0.044359708     0.041476009     0.032904353      0.04554337
3      0.044741874     0.040339266     0.035397362      0.04383978
4      0.058947875     0.048493601     0.055714271      0.05983738
5      0.050997929     0.050296474     0.050691042      0.06703274
6      0.057461696     0.038146557     0.050329606      0.04080140
7      0.056874446     0.049017564     0.031258144      0.05312191
8      0.064038069     0.049381373     0.046172672      0.03827328
9      0.050852581     0.041606657     0.032538539      0.04100212
10     0.051249725     0.072476728     0.071670454      0.05991351
11     0.046152216     0.041141923     0.053917413      0.03915900
12     0.046256431     0.038733544     0.033806363      0.05724544
13     0.063133685     0.063536176     0.074614279      0.02936272
14     0.035677992     0.029359784     0.044183272      0.04144391
15     0.021497572     0.054857600     0.036137945      0.02167276
16     0.075028057     0.085730709     0.047573427      0.09257384
17     0.044827767     0.034499030     0.041996540      0.03837886
18     0.061511381     0.042920525     0.050516011      0.05652006
19     0.015186750     0.031892174     0.050828152      0.02159366
20     0.064391001     0.062771950     0.051573515      0.07444261
21     0.057511673     0.058225762     0.037187586      0.06009338
22     0.056434920     0.085778644     0.039510532      0.03330086
23     0.005447772     0.002196682     0.009067904      0.00950724
   unlist(hold_id) unlist(hold_id) unlist(hold_id) unlist(hold_id)
1      0.036189631      0.03052170     0.030376836     0.030403998
2      0.030036543      0.03860820     0.028630442     0.037909120
3      0.033965779      0.02813159     0.032200814     0.025670442
4      0.058811834      0.02791992     0.033026716     0.051716740
5      0.049748463      0.02653471     0.054285996     0.058264558
6      0.042207573      0.02216123     0.034209283     0.024163611
7      0.051063220      0.03338198     0.054372752     0.049239068
8      0.045030591      0.03623385     0.043766145     0.060543848
9      0.054302022      0.03300084     0.042655714     0.043402500
10     0.062009822      0.03878912     0.055251777     0.054137653
11     0.052479414      0.02767014     0.030909421     0.063438484
12     0.033091378      0.02005067     0.036977167     0.046629467
13     0.067714423      0.04149497     0.039329181     0.049487302
14     0.040493969      0.01868952     0.031193390     0.044580922
15     0.057510386      0.04186896     0.025101356     0.040492517
16     0.048257740      0.05816659     0.054827143     0.066706813
17     0.042101400      0.03560006     0.056467257     0.048261941
18     0.058328511      0.03354070     0.045042438     0.035469715
19     0.064780979      0.03398035     0.034122728     0.032461677
20     0.036900867      0.05406955     0.045535805     0.053293757
21     0.082779456      0.05817083     0.062070860     0.066520187
22     0.128759284      0.04456350     0.102398636     0.060818215
23     0.007539013      0.00421763     0.004410938     0.004463658

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# Tue 18 Oct 2022 11:13:09 CEST
# MK
#could you calculate the other way around (Aut:X)? 
#And add a grey bar behind that [5.0,8.9] for humans, as well as a red line at 7.54 for bonobos?

# me
#Sorry, were these values also calculated using the proportion covered on the autosomes / proportion covered on the x?

# MK
#yes, they were. Im not sure what you are calculating there. Just taking the first example, chr1/chrX:
#0.028783297/0.007328132
#[1] 3.927781
#I suggest to use the average across autosomes, to get one value per individual.


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


> unlist(m_gbb)
 [1]  5.740464 12.663907 11.695184  8.451630  9.170435 22.753261  4.989008
 [8]  5.038376  7.093787  8.440210 10.024171 10.627367

> unlist(m_gbg)
[1] 31.190629  8.449803  3.432166 11.039526  3.893927  5.422719 31.505563
[8]  4.554739       Inf


#-----------------------------------------------------------------------------------------------------------------------


# plot violinplot of A:X for all individuals. - using proportions
x<-list(unlist(m_gbb),unlist(m_gbg))

# as.data.frame(unlist(x))
#     unlist(x)
#1  0.007918270
#2  0.003589299
#3  0.003886604


prop1_df<-cbind( as.data.frame(unlist(x)), as.data.frame(unlist(list(rep("MG", 12),rep("EL",9)))))
colnames(prop1_df)<-c("x_aut","pop")


# last entry is infinte (b/c division by 0) -> remove this row for plotting
prop1_df<-prop1_df[-21,]

#-----------------------------------------------------------------------------------------------------------------------


pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/plots/xaut_ratio_skovonly.17oct22.v1.pdf")

ggplot(prop1_df, aes(x=pop, y=x_aut,fill=pop)) + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#91C79B","#3A7E99")) + 
geom_boxplot(width = 0.07) + theme_classic()   + 
theme(legend.position="none") + scale_x_discrete(limits=c("MG", "EL")) + ylab("Autosome:X ratio") + xlab("Population")  +
 geom_hline(yintercept=7.54, color = "red") + 
  geom_hline(yintercept=5.0, linetype="dashed", color = "grey") +
  geom_hline(yintercept=8.9, linetype="dashed", color = "grey")
dev.off()



pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/plots/xaut_ratio_skovonly.17oct22.v2.pdf")

ggplot(prop1_df, aes(x=pop, y=x_aut,fill=pop)) + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#91C79B","#3A7E99")) + 
geom_boxplot(width = 0.07) + theme_classic()   + 
theme(legend.position="none") + scale_x_discrete(limits=c("MG", "EL")) + ylab("Autosome:X ratio") + xlab("Population")  +
 geom_hline(yintercept=7.54, color = "red") +
geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 5.0, ymax = 8.9),
            alpha = 0.01,
            fill = "grey")

dev.off()
# need to find a way to send the rectangle behind the violin plots *


# i think this version is better 
pdf("/scratch/devel/hpawar/admix/overlap.s*.skov/plots/xaut_ratio_skovonly.17oct22.v3.pdf")

ggplot(prop1_df, aes(x=pop, y=x_aut,fill=pop)) +
geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 5.0, ymax = 8.9),
            alpha = 0.01,
            fill = "grey") + 
geom_violin(alpha = 0.5) +  scale_fill_manual(values=c("#91C79B","#3A7E99")) + 
geom_boxplot(width = 0.07) + theme_classic()   + 
theme(legend.position="none") + scale_x_discrete(limits=c("MG", "EL")) + ylab("Autosome:X ratio") + xlab("Population")  +
 geom_hline(yintercept=7.54, color = "red") 
dev.off()


scp -r hpawar@172.16.10.21:/scratch/devel/hpawar/admix/'overlap.s*.skov'/plots/xaut_ratio_skovonly.17oct22.v2.pdf  /Users/harvi/Downloads/gorilla_postabc/overlap.sstar.skov


# https://stackoverflow.com/questions/20980622/change-background-color-on-just-one-region-of-the-graph


#-----------------------------------------------------------------------------------------------------------------------

# will need to increase text size & remove x axis (to match the popwise bp plots)
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

# have added this into the panel figure for the introgressed regions *


o1_2<-ggplot(ran_el, aes(x=id, y=protein_coding_bp,,fill=id)) + 
geom_violin(color="#91C79B",alpha = 0.5) +  scale_fill_manual(values=c("#91C79B")) + 
 geom_boxplot(width = 0.07,color="#91C79B") + theme_classic() + ylim(0.0186,0.037)  +
  geom_hline(yintercept=emp_pt_sk_gbg[[1]], color = "red")  + ylab("Protein coding bp")  + 
  theme(axis.title.x = element_blank()) +  theme(axis.text.x = element_text(size = 15)) + 
  theme(axis.title = element_text(size = 15))  +  theme(axis.text.y = element_text(size = 15))  +
     theme(legend.position = "none")  +   
theme(legend.position = "none") +
 theme(axis.line.y = element_blank(), 
   axis.line.x = element_line(),
   axis.title.y = element_blank(),
   axis.ticks.y=element_blank(),
   axis.text.y=element_blank()
   )
#-----------------------------------------------------------------------------------------------------------------------


prop1_df
       x_aut pop
1   5.740464  MG
2  12.663907  MG
3  11.695184  MG
4   8.451630  MG
5   9.170435  MG
6  22.753261  MG
7   4.989008  MG
8   5.038376  MG
9   7.093787  MG
10  8.440210  MG
11 10.024171  MG
12 10.627367  MG
13 31.190629  EL
14  8.449803  EL
15  3.432166  EL
16 11.039526  EL
17  3.893927  EL
18  5.422719  EL
19 31.505563  EL
20  4.554739  EL
#-----------------------------------------------------------------------------------------------------------------------