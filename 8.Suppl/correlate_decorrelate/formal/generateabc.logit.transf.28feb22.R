#!/usr/bin/env
#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
# Mon 28 Feb 2022 12:03:07 CET
# perform ABC Parameter Inference introducing logit transform in the ABC - to force posteriors to be within the priors

#-----------------------------------------------------------------------------------------------------------------------

library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/param_6oct", verbose=T)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat1_6oct", verbose=T)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct", verbose=T)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

lowerp<-c(3,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.19,5,0,0,0.1,10,1.5,10)
upperp<-c(100,20,30,20,20,25,5,20,0.25,30,0.6,50,100,100,6,100,15,100)
priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)


myabc_28feb22<-abc(target=target1,param=param,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)


myabc_28feb22[[17]][[1]]<-colnames(param)
myabc_28feb22[[17]][[2]]<-names(target1)

myabc_28feb22
summary(myabc_28feb22)

#-----------------------------------------------------------------------------------------------------------------------

save(myabc_28feb22,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/logit.transf_myabc_28feb22"))

q()

#-----------------------------------------------------------------------------------------------------------------------

# then plot on local

setwd("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_28feb22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/param_6oct",verbose=T)
library(abc)
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_posterior.28feb22.pdf") 
hist(myabc_28feb22)
dev.off() 
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_diagnostic.28feb22.pdf") 
plot(myabc_28feb22, param)
dev.off() 


#-----------------------------------------------------------------------------------------------------------------------
#Mon 28 Feb 2022 13:40:09 CET
# generate plots comparing the null models: null window ABC (no PLS, no logit transf) vs ABC_PLS (reduced dimensionality ABC) vs logit-transformed ABC (no PLS)

library(abc)
library(ggplot2)
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------
# read in the window abc model posteriors (no PLS, no logit transformation)
load(file=paste("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/myabc_6oct",sep=""), verbose=T)
#-----------------------------------------------------------------------------------------------------------------------
# read in PLS_ABC posteriors
load(file=paste("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_18feb22",sep=""), verbose=T)
#-----------------------------------------------------------------------------------------------------------------------
# read in logit transformed ABC posteriors (no PLS)
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_28feb22",verbose=T)
#-----------------------------------------------------------------------------------------------------------------------

PLS_ABC<-summary(myabc2)
window_ABC<-summary(myabc_6oct)
logit_ABC<-summary(myabc_28feb22)


process_stats<-function(a){
w1<-as.data.frame(window_ABC[,a])
w1$ID<-c("window_ABC")
colnames(w1)[1]<-'posteriors'
p1<-as.data.frame(PLS_ABC[,a])
p1$ID<-c("PLS_ABC")
colnames(p1)[1]<-'posteriors'
rownames(w1)<-NULL
rownames(p1)<-NULL
l1<-as.data.frame(logit_ABC[,a])
l1$ID<-c("logit_ABC")
colnames(l1)[1]<-'posteriors'
rownames(l1)<-NULL
vtes<-rbind(w1,p1,l1)
return(vtes)
}

#-----------------------------------------------------------------------------------------------------------------------

h1<-process_stats(1)
max(h1$posteriors)
v1<-ggplot(h1, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[1]) + ylim(3,120) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 

fin_v1<- v1 + geom_segment(aes(x='logit_ABC',y=3,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))


#w_lowl_t0   runif(1000,min=3,max=100)
#w_cros_t0   runif(1000,min=0.1,max=20)
h2<-process_stats(2)
cbind(min(h2$posteriors), max(h2$posteriors))

v2<-ggplot(h2, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[2]) + ylim(min(h2$posteriors), max(h2$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v2<- v2 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#e_lowl_t0   runif(1000,min=0.1,max=30)
h3<-process_stats(3)
cbind(min(h3$posteriors), max(h3$posteriors))

v3<-ggplot(h3, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[3]) + ylim(min(h3$posteriors), max(h3$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v3<- v3 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=30), position = position_nudge(x = -0.5))


#e_moun_t0   runif(1000,min=0.1,max=20)
h4<-process_stats(4)
cbind(min(h4$posteriors), max(h4$posteriors))

v4<-ggplot(h4, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[4]) + ylim(min(h4$posteriors), max(h4$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v4<- v4 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))


#e_lowl_t1   runif(1000,min=0.01,max=20) 
h5<-process_stats(5)
cbind(min(h5$posteriors), max(h5$posteriors))

v5<-ggplot(h5, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[5]) + ylim(0.01, max(h5$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v5<- v5 + geom_segment(aes(x='logit_ABC',y=0.01,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#e_lowl_t2   runif(1000,min=1,max=25)
h6<-process_stats(6)
cbind(min(h6$posteriors), max(h6$posteriors))

v6<-ggplot(h6, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[6]) + ylim(1, max(h6$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v6<- v6 + geom_segment(aes(x='logit_ABC',y=1,xend='logit_ABC',yend=25), position = position_nudge(x = -0.5))


#e_moun_t3   runif(1000,min=0.01,max=5)
h7<-process_stats(7)
cbind(min(h7$posteriors), max(h7$posteriors))

v7<-ggplot(h7, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[7]) + ylim(min(h7$posteriors), max(h7$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v7<- v7 + geom_segment(aes(x='logit_ABC',y=0.01,xend='logit_ABC',yend=5), position = position_nudge(x = -0.5))

#e_moun_t3.1 runif(1000,min=0.1,max=20)
h8<-process_stats(8)
cbind(min(h8$posteriors), max(h8$posteriors))

v8<-ggplot(h8, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[8]) + ylim(min(h8$posteriors), max(h8$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v8<- v8 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#t4  runif(1000,min=0.14,max=0.25) 
h9<-process_stats(9)
cbind(min(h9$posteriors), max(h9$posteriors))

v9<-ggplot(h9, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[9]) + ylim(0.14, max(h9$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v9<- v9 + geom_segment(aes(x='logit_ABC',y=0.14,xend='logit_ABC',yend=0.25), position = position_nudge(x = -0.5))

#e_anc_t4    runif(1000,min=0.1,max=30)
h10<-process_stats(10)
cbind(min(h10$posteriors), max(h10$posteriors))

v10<-ggplot(h10, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[10]) + ylim(min(h10$posteriors), max(h10$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v10<- v10 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=30), position = position_nudge(x = -0.5))

#t5  runif(1000,min=0.19,max=0.6)
h11<-process_stats(11)
cbind(min(h11$posteriors), max(h11$posteriors))

v11<-ggplot(h11, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[11]) + ylim(0.19, max(h11$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v11<- v11 + geom_segment(aes(x='logit_ABC',y=0.19,xend='logit_ABC',yend=0.6), position = position_nudge(x = -0.5))

#w_lowl_t5   runif(1000,min=5,max=50)
h12<-process_stats(12)
cbind(min(h12$posteriors), max(h12$posteriors))

v12<-ggplot(h12, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[12]) + ylim(min(h12$posteriors), max(h12$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v12<- v12 + geom_segment(aes(x='logit_ABC',y=5,xend='logit_ABC',yend=50), position = position_nudge(x = -0.5))

#admix_w_e_t6    runif(1000,min=0,max=100)
h13<-process_stats(13)
cbind(min(h13$posteriors), max(h13$posteriors))

v13<-ggplot(h13, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[13]) + ylim(min(h13$posteriors), max(h13$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v13<- v13 + geom_segment(aes(x='logit_ABC',y=0,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#admix_e_w_t6    runif(1000,min=0,max=100)
h14<-process_stats(14)
cbind(min(h14$posteriors), max(h14$posteriors))

v14<-ggplot(h14, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[14]) + ylim(min(h14$posteriors), max(h14$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v14<- v14 + geom_segment(aes(x='logit_ABC',y=0,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#t7  runif(1000,min=0.1,max=6) 
h15<-process_stats(15)
cbind(min(h15$posteriors), max(h15$posteriors))

v15<-ggplot(h15, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[15]) + ylim(0.1, max(h15$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v15<- v15 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=6), position = position_nudge(x = -0.5))

#w_anc_t7    runif(1000,min=10,max=100)
h16<-process_stats(16)
cbind(min(h16$posteriors), max(h16$posteriors))

v16<-ggplot(h16, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[16]) + ylim(10, max(h16$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v16<- v16 + geom_segment(aes(x='logit_ABC',y=10,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#t8  runif(1000,min=1.5,max=15)
h17<-process_stats(17)
cbind(min(h17$posteriors), max(h17$posteriors))

v17<-ggplot(h17, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[17]) + ylim(1.5, max(h17$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v17<- v17 + geom_segment(aes(x='logit_ABC',y=1.5,xend='logit_ABC',yend=15), position = position_nudge(x = -0.5))

#gor_anc runif(1000,min=10,max=100)
h18<-process_stats(18)
cbind(min(h18$posteriors), max(h18$posteriors))

v18<-ggplot(h18, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[18]) + ylim(min(h18$posteriors), max(h18$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v18<- v18 + geom_segment(aes(x='logit_ABC',y=10,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))


library(ggpubr)

ggarrange(fin_v1,fin_v2,fin_v3,fin_v4,fin_v5,fin_v6,fin_v7,fin_v8,fin_v9,fin_v10,fin_v11,fin_v12,fin_v13,fin_v14,fin_v15,fin_v16,fin_v17,fin_v18, ncol=6, nrow=3, common.legend = TRUE, legend="bottom")

