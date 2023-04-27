#!/usr/bin/r
options(scipen=100)
library(magrittr)
library(admixtools)
library(tidyverse)
require(gridExtra) 
library("gridGraphics")

# run once to calculate f2 stats
extract_f2(pref="/scratch/admixlab/gorilla/genos/gori",outdir="/scratch/admixlab/gorilla/genos/goris/",blgsize=500000,overwrite=T)

leaves<-c("MG","ELG","WLG","CRG","Orang")
f2_blocks = f2_from_precomp("/scratch/admixlab/gorilla/genos/goris/",pops=leaves)

## get certain number of admixture nodes
ev<-c(0:4)
allres<-list()

igr<-list()
for (nm in ev) {
  res_NUMADMIX = find_graphs(f2_blocks, outpop = 'Orang',numadmix = nm)
  
  #Find best graph
  winner_NUMADMIX = res_NUMADMIX %>% slice_min(score, with_ties = FALSE)
  winner_NUMADMIX$score[[1]]
  
  allres[[nm+1]]<-res_NUMADMIX
  igr[[nm+1]]<-winner_NUMADMIX$edges[[1]]
  save(allres,file="/scratch/admixlab/gorilla/genos/res/ft2stat_results2")
}

## plot all of them
apl=list()
#apl[[1]]<-plot_graph(gph)
for (nm in 0:5) {
  nn=nm+1
  res_NUMADMIX<-allres[[nn]]
  winner_NUMADMIX = res_NUMADMIX %>% slice_min(score, with_ties = FALSE)
  winner_NUMADMIX$score[[1]]
  apl[[nn]]<-plot_graph(winner_NUMADMIX$edges[[1]],title=paste("edges:",nm))
}

pdf(paste("/scratch/admixlab/gorilla/genos/res/f2_5nice.pdf",sep=""),10,15)
grid.arrange(grobs=apl,ncol=2,nrow=3)
dev.off()

# compare fit for increasing number of edges
pops = dimnames(f2_blocks)[[1]]
afit<-list()
for (j in c(1:4)) {
  print(j)
  graph1 = allres[[j]] %>% slice_min(score, with_ties = FALSE)
  graph2 = allres[[j+1]] %>% slice_min(score, with_ties = FALSE)
  fits = qpgraph_resample_multi(f2_blocks, list(graph1$graph[[1]], graph2$graph[[1]]), nboot = 100)
  afit[[j+1]]<-compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)
}

afits<-do.call(rbind,afit)

## explicit test of f4-statistics
f4t1<-qpdstat("/scratch/admixlab/gorilla/genos/goris",pop1=c("MG"), pop2=c("ELG"),pop3=c("WLG","CRG"), pop4=c("Orang"))
f4t2<-qpdstat("/scratch/admixlab/gorilla/genos/goris",pop1=c("WLG"), pop2=c("CRG"),pop3=c("MG","ELG"), pop4=c("Orang"))
write.table(rbind(f4t1,f4t2),file="/scratch/admixlab/gorilla/genos/res/f4tab.txt",sep="\t",col.names=T,row.names=F,quote=F)


## explicit test of graph with ghost admixture edge vs. 0 and 1 edge best graph
g = matrix(c('R', 'Orang', 'R', 'GAnc', 'GAnc', 'Ggh', 'Ggh', 'EGAnc', 'GAnc','GAnc2','GAnc2', 'EGAnc', 'EGAnc', 'ELG','EGAnc','MG','GAnc2','WAnc','WAnc','WLG','WAnc','CRG' ), , 2, byrow = T) %>% edges_to_igraph()

graph1 = allres[[1]] %>% slice_min(score, with_ties = FALSE)
graph2 = allres[[2]] %>% slice_min(score, with_ties = FALSE)

fits = qpgraph_resample_multi(f2_blocks, list(graph1$graph[[1]], g), nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test
             
fits2 = qpgraph_resample_multi(f2_blocks, list(graph2$graph[[1]], g), nboot = 100)
compare_fits(fits2[[1]]$score_test, fits2[[2]]$score_test
             
app<-plot_graph(g,title=paste("ghost model"))

pdf(paste("/scratch/admixlab/gorilla/genos/res/f2_ghost.pdf",sep=""),6,6)
app
dev.off()

  
