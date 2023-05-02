#module load R/4.0.1
#!/usr/bin/r

#-----------------------------------------------------------------------------------------------------------------------
require(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
a=args[1]
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# number of individuals
nind <- 21

#-----------------------------------------------------------------------------------------------------------------------

# number of SNPs per GT file
nrows<-c(8649941,10590950,8801570,8064032,7794661,7592642,6306813,6322633,1325944,5550879,5110818,5372361,4328829,3658387,3259675,2671532,2406836,3586343,946585,2275219,1367710,782934)

input<-paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_",a,"_Egor.biallelic.filt.GT")
searchString <- ' '
replacementString <- ''
input = gsub(searchString,replacementString,input)
# read in the GT matrix for this chr
md_geno <- matrix(scan(file=input), ncol=nind, nrow=nrows[[a]], byrow=T)


obs.sfs<-table(rowSums(md_geno))

obs.sfs # should be 0 -42

# then only normalise once have calculated site categories for all chromosomes

# only once
savelist<-list()
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))


# in every job
load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))
savelist[[a]]<-obs.sfs 
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))


q()
#-----------------------------------------------------------------------------------------------------------------------

# Run the following interactively - 

countsites<-function(a){

nind <- 21
# number of SNPs per GT file
nrows<-c(8649941,10590950,8801570,8064032,7794661,7592642,6306813,6322633,1325944,5550879,5110818,5372361,4328829,3658387,3259675,2671532,2406836,3586343,946585,2275219,1367710,782934)

input<-paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/4_",a,"_Egor.biallelic.filt.GT")
searchString <- ' '
replacementString <- ''
input = gsub(searchString,replacementString,input)
# read in the GT matrix for this chr
md_geno <- matrix(scan(file=input), ncol=nind, nrow=nrows[[a]], byrow=T)

obs.sfs<-table(rowSums(md_geno))
obs.sfs # should be 0 -42

# in every job
load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))
savelist[[a]]<-obs.sfs 
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))

}

countsites(3)

load(file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""),verbose=T)
savelist # check chr 3 also added - yes worked

for (i in 4:22){
countsites(i)
}    

savelist
save(savelist,file=paste("/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfscounts",sep=""))


#-----------------------------------------------------------------------------------------------------------------------

all<-data.frame(matrix(unlist(savelist), ncol=length(savelist[[1]]), byrow=TRUE),stringsAsFactors=FALSE)

aut.obs.sfs <-colSums(all)

aut.obs.sfs

sum(aut.obs.sfs)


# normalise 
test<-aut.obs.sfs/sum(aut.obs.sfs)

names(test)<-c(0, seq(1:42))


# drop first category
test1<-test[-1]

x<-as.data.table(test1)

x$V1<-as.numeric(names(test1))

write.table(x,file="/scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfs.input.txt",sep="\t",row.names=FALSE,col.names=FALSE)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# choosing how many blocks to split each chr into **

require(data.table)
options(stringsAsFactors=F)
hg19<-fread('/home/devel/marcmont/scratch/Harvi/geneflow/test/test.reconst.ids/ora_1_1/proportion/hg19.autosomes.chrom.sizes.txt')

hg19.lengths<-apply(hg19[,2],2,function(x) as.numeric(as.character(x)))

hg19[c(1:18,20,19,22,21),]

hg19<-hg19[c(1:18,20,19,22,21),]
hg19$V3<-hg19[,2]/51304566

#hg19
#       V1        V2        V3
# 1:  chr1 249250621 4.8582542   # 50 blocks
# 2:  chr2 243199373 4.7403066   # 50 blocks
# 3:  chr3 198022430 3.8597428   # 40 blocks
# 4:  chr4 191154276 3.7258726   # 40 blocks
# 5:  chr5 180915260 3.5262994   # 35 blocks
# 6:  chr6 171115067 3.3352795   # 35 blocks
# 7:  chr7 159138663 3.1018421   # 30 blocks
# 8:  chr8 146364022 2.8528459   # 30 blocks
# 9:  chr9 141213431 2.7524535   # 30 blocks
#10: chr10 135534747 2.6417677   # 25 blocks
#11: chr11 135006516 2.6314717   # 25 blocks
#12: chr12 133851895 2.6089665   # 25 blocks
#13: chr13 115169878 2.2448271   # 20 blocks
#14: chr14 107349540 2.0923974   # 20 blocks
#15: chr15 102531392 1.9984847   # 20 blocks
#16: chr16  90354753 1.7611445   # 20 blocks
#17: chr17  81195210 1.5826118   # 15 blocks
#18: chr18  78077248 1.5218382   # 15 blocks
#19: chr19  59128983 1.1525092   # 10 blocks
#20: chr20  63025520 1.2284583   # 10 blocks
#21: chr21  48129895 0.9381211   # 10 blocks
#22: chr22  51304566 1.0000000   # 10 blocks


# ie most should be 10 blocks, but could increase to 50/40 for the larger chr? # think that makes sense
#-----------------------------------------------------------------------------------------------------------------------
