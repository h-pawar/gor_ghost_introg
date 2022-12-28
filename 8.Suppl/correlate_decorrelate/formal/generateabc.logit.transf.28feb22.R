#!/usr/bin/env
#module load R/4.0.1
#-----------------------------------------------------------------------------------------------------------------------
# Mon 28 Feb 2022 12:03:07 CET
# reperform ABC Parameter Inference introducing logit transform in the ABC - to force posteriors to be within the priors

# read in these files & try to run the logit transf
	# files generated from - /Volumes/"Ultra USB 3.0"/IBE/further.analysis.feb.2020/gorillas/abc/generateabc.6oct.R
		# ie window-based ABC, null model
#save(param,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/param_6oct"))
#save(sumstat1,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat1_6oct"))
#save(target1,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct"))
#-----------------------------------------------------------------------------------------------------------------------

library(abc)
options(scipen=100)
require(data.table)
options(stringsAsFactors=F)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/param_6oct", verbose=T)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/sumstat1_6oct", verbose=T)
load("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/target1_6oct", verbose=T)


#> str(param)
# num [1:35543, 1:18] 97.6 32.3 77.7 50.3 52.7 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : NULL
#  ..$ : chr [1:18] "w_lowl_t0" "w_cros_t0" "e_lowl_t0" "e_moun_t0" ...

#   str( sumstat1)
# num [1:35543, 1:40] 3.11 2.79 2.08 3.32 1.14 ...
# - attr(*, "dimnames")=List of 2
#  ..$ : NULL
#  ..$ : chr [1:40] "X1" "X2" "X3" "X4" ...

#  > str(target1)
# Named num [1:40] 0.911 0.701 0.412 0.43 0.062 ...
# - attr(*, "names")= chr [1:40] "het_mu_WL" "het_mu_WC" "het_mu_EL" "het_mu_EM" ...

#-----------------------------------------------------------------------------------------------------------------------

# referring back to the manual
#https://cran.r-project.org/web/packages/abc/abc.pdf

# abc(target, param, sumstat, tol, method, hcorr = TRUE, transf = "none",
#logit.bounds, subset = NULL, kernel = "epanechnikov", numnet =
#10, sizenet = 5, lambda = c(0.0001,0.001,0.01), trace = FALSE, maxit =
#500, ...)

#transf = a vector of character strings indicating the kind of transformation to be applied
#to the parameter values. The possible values are "log", "logit", and "none"
#(default), when no is transformation applied. See also Details

#logit.bounds = a matrix of bounds if transf is "logit". 
#The matrix has as many lines as parameters (including the ones that are not "logit" transformed) and 2 columns.
#First column is the minimum bound and second column is the maximum bound.
#-----------------------------------------------------------------------------------------------------------------------

#str(param)
# num [1:35543, 1:18]
 # ie matrix of the priors - need to refer back tp

# refering back to the priors used to generate this simulated data
# /scratch/devel/hpawar/admix/abc/simul/scripts/test.abc.model.v5.R # new version supercedes test.abc.model.v4.R (used to generate prev 1-2000 simns

lowerp<-c(3,0.1,0.1,0.1,0.01,1,0.01,0.1,0.14,0.1,0.19,5,0,0,0.1,10,1.5,10)
upperp<-c(100,20,30,20,20,25,5,20,0.25,30,0.6,50,100,100,6,100,15,100)
priordf<-cbind(lowerp,upperp)
prior_ranges<-data.matrix(priordf)


myabc_28feb22<-abc(target=target1,param=param,tol=0.005,sumstat=sumstat1,method="neuralnet",numnet=100,transf="logit",logit.bounds=prior_ranges)
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
123456789101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100
Warning message:
All parameters are "logit" transformed. 

myabc_28feb22[[17]][[1]]<-colnames(param)
myabc_28feb22[[17]][[2]]<-names(target1)

myabc_28feb22
summary(myabc_28feb22)

q()
#-----------------------------------------------------------------------------------------------------------------------


> myabc_28feb22
Call:
abc(target = target1, param = param, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Method:
Non-linear regression via neural networks
with correction for heteroscedasticity

Parameters:
w_lowl_t0, w_cros_t0, e_lowl_t0, e_moun_t0, e_lowl_t1, e_lowl_t2, e_moun_t3, e_moun_t3.1, t4, e_anc_t4, t5, w_lowl_t5, admix_w_e_t6, admix_e_w_t6, t7, w_anc_t7, t8, gor_anc

Statistics:
het_mu_WL, het_mu_WC, het_mu_EL, het_mu_EM, het_sd_WL, het_sd_EL, het_sd_EM, fixedsites_WL, fixedsites_WC, fixedsites_EL, fixedsites_EM, segsites_WL, segsites_WC, segsites_EL, segsites_EM, fixedperid_mu_WL, fixedperid_mu_WC, fixedperid_mu_EL, fixedperid_mu_EM, fixedperid_sd_WL, fixedperid_sd_EL, fixedperid_sd_EM, fst_WL.WC, fst_WL.EL, fst_WL.EM, fst_WC.EL, fst_WC.EM, fst_EL.EM, pi_mu_WL, pi_mu_EL, pi_mu_EM, pi_sd_WL, pi_sd_EL, pi_sd_EM, tajima_mu_WL, tajima_mu_EL, tajima_mu_EM, tajima_sd_WL, tajima_sd_EL, tajima_sd_EM

Total number of simulations 35543 

Number of accepted simulations:  178 

> summary(myabc_28feb22)
Call: 
abc(target = target1, param = param, sumstat = sumstat1, tol = 0.005, 
    method = "neuralnet", transf = "logit", logit.bounds = prior_ranges, 
    numnet = 100)
Data:
 abc.out$adj.values (178 posterior samples)
Weights:
 abc.out$weights

                       w_lowl_t0 w_cros_t0 e_lowl_t0 e_moun_t0 e_lowl_t1
Min.:                    25.7410    1.4003    0.1087    0.1290    5.1680
Weighted 2.5 % Perc.:    31.6450    1.8195    0.3933    0.1613    6.7757
Weighted Median:         71.3598    8.5047   14.3944    0.3063   12.4350
Weighted Mean:           70.4555    9.7955   14.1856    0.6191   12.4301
Weighted Mode:           88.4138    7.5839    6.6926    0.2627   13.2578
Weighted 97.5 % Perc.:   97.2996   19.5923   28.4326    2.9315   18.7555
Max.:                    98.5464   19.8937   29.8618   13.0511   19.9133
                       e_lowl_t2 e_moun_t3 e_moun_t3.1      t4 e_anc_t4      t5
Min.:                     4.1169    0.0510      0.7070  0.1468   2.8892  0.1910
Weighted 2.5 % Perc.:    11.2885    0.1192      1.9213  0.1543   6.2215  0.2474
Weighted Median:         22.4913    0.9068      8.7514  0.2108  14.4836  0.4605
Weighted Mean:           21.1732    1.1892      9.5399  0.2098  15.4194  0.4434
Weighted Mode:           23.5678    0.5809      6.0940  0.2298  11.5730  0.4905
Weighted 97.5 % Perc.:   24.8410    3.5420     19.1953  0.2462  28.9200  0.5867
Max.:                    24.9976    4.5624     19.8338  0.2498  29.8302  0.5974
                       w_lowl_t5 admix_w_e_t6 admix_e_w_t6      t7 w_anc_t7
Min.:                     5.0723       2.0913       0.1441  1.1182  14.5561
Weighted 2.5 % Perc.:     6.9223      11.9353       2.4502  2.0377  20.0460
Weighted Median:         26.3827      73.4137      64.9857  3.9768  43.2303
Weighted Mean:           27.8733      64.8076      58.1413  4.0190  46.5845
Weighted Mode:           42.6831      87.5613      88.3961  4.7855  41.5012
Weighted 97.5 % Perc.:   47.6586      98.1674      98.9707  5.5838  87.2841
Max.:                    49.9164      99.9281      99.9649  5.9772  98.6036
                            t8 gor_anc
Min.:                   5.0783 10.3031
Weighted 2.5 % Perc.:   6.0933 10.8609
Weighted Median:       11.6288 16.6371
Weighted Mean:         11.2930 17.0940
Weighted Mode:         11.8299 12.9181
Weighted 97.5 % Perc.: 14.5678 27.1546
Max.:                  14.7361 73.7401

save(myabc_28feb22,file=paste("/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/logit.transf_myabc_28feb22"))


# transfer from cluster->local
#scp -r hpawar@172.16.10.20:/scratch/devel/hpawar/admix/abc/results/test/11sep21simns/segrecalc/logit.transf_myabc_28feb22 /Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc
#-----------------------------------------------------------------------------------------------------------------------

# then plot on local

setwd("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc")
require(data.table)

# A) tol 0.005
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_28feb22",verbose=T)
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/param_6oct",verbose=T)
library(abc)
# histogram of the weighted posterior sample  : https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_posterior.28feb22.pdf") 
hist(myabc_28feb22)
dev.off() 
# generate the diagnostic plots of the estimation of the posterior distribution 
# a density plot of the prior distribution, a density plot of the posterior distribution estimated with and without regression- based correction, a scatter plot of the Euclidean distances as a function of the parameter values, and a normal Q-Q plot of the residuals from the regression
pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_diagnostic.28feb22.pdf") 
plot(myabc_28feb22, param)
dev.off() 



# 2) plot on local - violin plots of posterior distributions
#-----------------------------------------------------------------------------------------------------------------------
#Mon 28 Feb 2022 13:40:09 CET
# read in logit transformed ABC posteriors (no PLS)
# & add as further violin plot comparison

#   # null models: null window ABC vs ABC_PLS (reduced dimensionality ABC) vs logit-transformed ABC (no PLS)

library(abc)
library(ggplot2)
require(data.table)
options(stringsAsFactors=F)

#-----------------------------------------------------------------------------------------------------------------------
# read in the window abc model posteriors
# /Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/myabc_6oct
load(file=paste("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/myabc_6oct",sep=""), verbose=T)
#Loading objects:
#  myabc_6oct
#-----------------------------------------------------------------------------------------------------------------------

# read in PLS_ABC posteriors
load(file=paste("/Users/harvi/Downloads/gorilla_abc/modelchoice/pls/abc_pls_18feb22",sep=""), verbose=T)
#Loading objects:
#  myabc2
#-----------------------------------------------------------------------------------------------------------------------
# read in logit transformed ABC posteriors (no PLS)
load("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/logit.transf_myabc_28feb22",verbose=T)
#Loading objects:
#  myabc_28feb22
#-----------------------------------------------------------------------------------------------------------------------

PLS_ABC<-summary(myabc2)
window_ABC<-summary(myabc_6oct)
logit_ABC<-summary(myabc_28feb22)


# make the first part into a function
#-----------------------------------------------------------------------------------------------------------------------

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

h1<-process_stats(1)
max(h1$posteriors)
v1<-ggplot(h1, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[1]) + ylim(3,120) + theme_classic() + geom_boxplot(width=0.1)

#https://stackoverflow.com/questions/68940577/line-segment-across-bar-plot-with-multiple-series-ggplot2

# adds line to show the prior ranges
#fin_v1<- v1 + geom_segment(aes(x='PLS_ABC',y=3,xend='PLS_ABC',yend=100), position = position_nudge(x = -0.5))
fin_v1<- v1 + geom_segment(aes(x='logit_ABC',y=3,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))


#w_lowl_t0   runif(1000,min=3,max=100)
#w_cros_t0   runif(1000,min=0.1,max=20)
h2<-process_stats(2)
cbind(min(h2$posteriors), max(h2$posteriors))
#   [,1]     [,2]
#[1,] -0.199427 32.96607
v2<-ggplot(h2, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[2]) + ylim(min(h2$posteriors), max(h2$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v2<- v2 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#e_lowl_t0   runif(1000,min=0.1,max=30)
h3<-process_stats(3)
cbind(min(h3$posteriors), max(h3$posteriors))
#  [,1]     [,2]
#[1,] -4.823436 33.79758
v3<-ggplot(h3, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[3]) + ylim(min(h3$posteriors), max(h3$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v3<- v3 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=30), position = position_nudge(x = -0.5))


#e_moun_t0   runif(1000,min=0.1,max=20)
h4<-process_stats(4)
cbind(min(h4$posteriors), max(h4$posteriors))
#    [,1]     [,2]
#[1,] -2.239405 29.27488
v4<-ggplot(h4, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[4]) + ylim(min(h4$posteriors), max(h4$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v4<- v4 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))


#e_lowl_t1   runif(1000,min=0.01,max=20) 
h5<-process_stats(5)
cbind(min(h5$posteriors), max(h5$posteriors))
#      [,1]     [,2]
#[1,] 2.365969 40.33518
v5<-ggplot(h5, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[5]) + ylim(0.01, max(h5$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v5<- v5 + geom_segment(aes(x='logit_ABC',y=0.01,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#e_lowl_t2   runif(1000,min=1,max=25)
h6<-process_stats(6)
cbind(min(h6$posteriors), max(h6$posteriors))
#      [,1]    [,2]
#[1,] 5.844928 31.2664
v6<-ggplot(h6, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[6]) + ylim(1, max(h6$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v6<- v6 + geom_segment(aes(x='logit_ABC',y=1,xend='logit_ABC',yend=25), position = position_nudge(x = -0.5))


#e_moun_t3   runif(1000,min=0.01,max=5)
h7<-process_stats(7)
cbind(min(h7$posteriors), max(h7$posteriors))
#        [,1]     [,2]
#[1,] -1.582033 5.580239
v7<-ggplot(h7, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[7]) + ylim(min(h7$posteriors), max(h7$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v7<- v7 + geom_segment(aes(x='logit_ABC',y=0.01,xend='logit_ABC',yend=5), position = position_nudge(x = -0.5))

#e_moun_t3.1 runif(1000,min=0.1,max=20)
h8<-process_stats(8)
cbind(min(h8$posteriors), max(h8$posteriors))
# [,1]     [,2]
#[1,] -7.143765 23.97649
v8<-ggplot(h8, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[8]) + ylim(min(h8$posteriors), max(h8$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v8<- v8 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#t4  runif(1000,min=0.14,max=0.25) 
h9<-process_stats(9)
cbind(min(h9$posteriors), max(h9$posteriors))
#       [,1]      [,2]
#[1,] 0.1536823 0.2906091
v9<-ggplot(h9, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[9]) + ylim(0.14, max(h9$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v9<- v9 + geom_segment(aes(x='logit_ABC',y=0.14,xend='logit_ABC',yend=0.25), position = position_nudge(x = -0.5))

#e_anc_t4    runif(1000,min=0.1,max=30)
h10<-process_stats(10)
cbind(min(h10$posteriors), max(h10$posteriors))
#  [,1]     [,2]
#[1,] -5.325377 43.71968
v10<-ggplot(h10, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[10]) + ylim(min(h10$posteriors), max(h10$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v10<- v10 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=30), position = position_nudge(x = -0.5))

#t5  runif(1000,min=0.19,max=0.6)
h11<-process_stats(11)
cbind(min(h11$posteriors), max(h11$posteriors))
#      [,1]      [,2]
#[1,] 0.2135755 0.6972605
v11<-ggplot(h11, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[11]) + ylim(0.19, max(h11$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v11<- v11 + geom_segment(aes(x='logit_ABC',y=0.19,xend='logit_ABC',yend=0.6), position = position_nudge(x = -0.5))

#w_lowl_t5   runif(1000,min=5,max=50)
h12<-process_stats(12)
cbind(min(h12$posteriors), max(h12$posteriors))
#  [,1]     [,2]
#[1,] -4.332538 70.06947
v12<-ggplot(h12, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[12]) + ylim(min(h12$posteriors), max(h12$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v12<- v12 + geom_segment(aes(x='logit_ABC',y=5,xend='logit_ABC',yend=50), position = position_nudge(x = -0.5))

#admix_w_e_t6    runif(1000,min=0,max=100)
h13<-process_stats(13)
cbind(min(h13$posteriors), max(h13$posteriors))
#   [,1]     [,2]
#[1,] -28.19599 165.7739
v13<-ggplot(h13, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[13]) + ylim(min(h13$posteriors), max(h13$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v13<- v13 + geom_segment(aes(x='logit_ABC',y=0,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#admix_e_w_t6    runif(1000,min=0,max=100)
h14<-process_stats(14)
cbind(min(h14$posteriors), max(h14$posteriors))
#    [,1]     [,2]
#[1,] -49.94324 146.6777
v14<-ggplot(h14, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[14]) + ylim(min(h14$posteriors), max(h14$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v14<- v14 + geom_segment(aes(x='logit_ABC',y=0,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#t7  runif(1000,min=0.1,max=6) 
h15<-process_stats(15)
cbind(min(h15$posteriors), max(h15$posteriors))
#    [,1]     [,2]
#[1,] 2.335279 9.976224
v15<-ggplot(h15, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[15]) + ylim(0.1, max(h15$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v15<- v15 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=6), position = position_nudge(x = -0.5))

#w_anc_t7    runif(1000,min=10,max=100)
h16<-process_stats(16)
cbind(min(h16$posteriors), max(h16$posteriors))
#   [,1]     [,2]
#[1,] 19.03338 102.1845
v16<-ggplot(h16, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[16]) + ylim(10, max(h16$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v16<- v16 + geom_segment(aes(x='logit_ABC',y=10,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#t8  runif(1000,min=1.5,max=15)
h17<-process_stats(17)
cbind(min(h17$posteriors), max(h17$posteriors))
#     [,1]     [,2]
#[1,] 4.080417 18.00849
v17<-ggplot(h17, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[17]) + ylim(1.5, max(h17$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v17<- v17 + geom_segment(aes(x='logit_ABC',y=1.5,xend='logit_ABC',yend=15), position = position_nudge(x = -0.5))

#gor_anc runif(1000,min=10,max=100)
h18<-process_stats(18)
cbind(min(h18$posteriors), max(h18$posteriors))
#     [,1]     [,2]
#[1,] -1.279905 123.2393
v18<-ggplot(h18, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[18]) + ylim(min(h18$posteriors), max(h18$posteriors)) + theme_classic() + geom_boxplot(width=0.1)
fin_v18<- v18 + geom_segment(aes(x='logit_ABC',y=10,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

# arrange all in one figure

library(ggpubr)

#for i in {1..18}
#do
#echo "fin_v$i"
#done

ggarrange(fin_v1,fin_v2,fin_v3,fin_v4,fin_v5,fin_v6,fin_v7,fin_v8,fin_v9,fin_v10,fin_v11,fin_v12,fin_v13,fin_v14,fin_v15,fin_v16,fin_v17,fin_v18, 
          ncol = 5, nrow = 4)


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


# Tue 27 Sep 2022 14:11:44 CEST
# replot to separate the x axis labels (rmd & one legend), also changed colour scheme
# replotted to add into the paper

h1<-process_stats(1)
max(h1$posteriors)
v1<-ggplot(h1, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[1]) + ylim(3,120) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 

# removing the axis ticks may be best (for clarity)

#https://stackoverflow.com/questions/68940577/line-segment-across-bar-plot-with-multiple-series-ggplot2

# adds line to show the prior ranges
#fin_v1<- v1 + geom_segment(aes(x='PLS_ABC',y=3,xend='PLS_ABC',yend=100), position = position_nudge(x = -0.5))
fin_v1<- v1 + geom_segment(aes(x='logit_ABC',y=3,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))


#w_lowl_t0   runif(1000,min=3,max=100)
#w_cros_t0   runif(1000,min=0.1,max=20)
h2<-process_stats(2)
cbind(min(h2$posteriors), max(h2$posteriors))
#   [,1]     [,2]
#[1,] -0.199427 32.96607
v2<-ggplot(h2, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[2]) + ylim(min(h2$posteriors), max(h2$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v2<- v2 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#e_lowl_t0   runif(1000,min=0.1,max=30)
h3<-process_stats(3)
cbind(min(h3$posteriors), max(h3$posteriors))
#  [,1]     [,2]
#[1,] -4.823436 33.79758
v3<-ggplot(h3, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[3]) + ylim(min(h3$posteriors), max(h3$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v3<- v3 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=30), position = position_nudge(x = -0.5))


#e_moun_t0   runif(1000,min=0.1,max=20)
h4<-process_stats(4)
cbind(min(h4$posteriors), max(h4$posteriors))
#    [,1]     [,2]
#[1,] -2.239405 29.27488
v4<-ggplot(h4, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[4]) + ylim(min(h4$posteriors), max(h4$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v4<- v4 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))


#e_lowl_t1   runif(1000,min=0.01,max=20) 
h5<-process_stats(5)
cbind(min(h5$posteriors), max(h5$posteriors))
#      [,1]     [,2]
#[1,] 2.365969 40.33518
v5<-ggplot(h5, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[5]) + ylim(0.01, max(h5$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v5<- v5 + geom_segment(aes(x='logit_ABC',y=0.01,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#e_lowl_t2   runif(1000,min=1,max=25)
h6<-process_stats(6)
cbind(min(h6$posteriors), max(h6$posteriors))
#      [,1]    [,2]
#[1,] 5.844928 31.2664
v6<-ggplot(h6, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[6]) + ylim(1, max(h6$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v6<- v6 + geom_segment(aes(x='logit_ABC',y=1,xend='logit_ABC',yend=25), position = position_nudge(x = -0.5))


#e_moun_t3   runif(1000,min=0.01,max=5)
h7<-process_stats(7)
cbind(min(h7$posteriors), max(h7$posteriors))
#        [,1]     [,2]
#[1,] -1.582033 5.580239
v7<-ggplot(h7, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[7]) + ylim(min(h7$posteriors), max(h7$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v7<- v7 + geom_segment(aes(x='logit_ABC',y=0.01,xend='logit_ABC',yend=5), position = position_nudge(x = -0.5))

#e_moun_t3.1 runif(1000,min=0.1,max=20)
h8<-process_stats(8)
cbind(min(h8$posteriors), max(h8$posteriors))
# [,1]     [,2]
#[1,] -7.143765 23.97649
v8<-ggplot(h8, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[8]) + ylim(min(h8$posteriors), max(h8$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v8<- v8 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=20), position = position_nudge(x = -0.5))

#t4  runif(1000,min=0.14,max=0.25) 
h9<-process_stats(9)
cbind(min(h9$posteriors), max(h9$posteriors))
#       [,1]      [,2]
#[1,] 0.1536823 0.2906091
v9<-ggplot(h9, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[9]) + ylim(0.14, max(h9$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v9<- v9 + geom_segment(aes(x='logit_ABC',y=0.14,xend='logit_ABC',yend=0.25), position = position_nudge(x = -0.5))

#e_anc_t4    runif(1000,min=0.1,max=30)
h10<-process_stats(10)
cbind(min(h10$posteriors), max(h10$posteriors))
#  [,1]     [,2]
#[1,] -5.325377 43.71968
v10<-ggplot(h10, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[10]) + ylim(min(h10$posteriors), max(h10$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v10<- v10 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=30), position = position_nudge(x = -0.5))

#t5  runif(1000,min=0.19,max=0.6)
h11<-process_stats(11)
cbind(min(h11$posteriors), max(h11$posteriors))
#      [,1]      [,2]
#[1,] 0.2135755 0.6972605
v11<-ggplot(h11, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[11]) + ylim(0.19, max(h11$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v11<- v11 + geom_segment(aes(x='logit_ABC',y=0.19,xend='logit_ABC',yend=0.6), position = position_nudge(x = -0.5))

#w_lowl_t5   runif(1000,min=5,max=50)
h12<-process_stats(12)
cbind(min(h12$posteriors), max(h12$posteriors))
#  [,1]     [,2]
#[1,] -4.332538 70.06947
v12<-ggplot(h12, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[12]) + ylim(min(h12$posteriors), max(h12$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v12<- v12 + geom_segment(aes(x='logit_ABC',y=5,xend='logit_ABC',yend=50), position = position_nudge(x = -0.5))

#admix_w_e_t6    runif(1000,min=0,max=100)
h13<-process_stats(13)
cbind(min(h13$posteriors), max(h13$posteriors))
#   [,1]     [,2]
#[1,] -28.19599 165.7739
v13<-ggplot(h13, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[13]) + ylim(min(h13$posteriors), max(h13$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v13<- v13 + geom_segment(aes(x='logit_ABC',y=0,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#admix_e_w_t6    runif(1000,min=0,max=100)
h14<-process_stats(14)
cbind(min(h14$posteriors), max(h14$posteriors))
#    [,1]     [,2]
#[1,] -49.94324 146.6777
v14<-ggplot(h14, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[14]) + ylim(min(h14$posteriors), max(h14$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v14<- v14 + geom_segment(aes(x='logit_ABC',y=0,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#t7  runif(1000,min=0.1,max=6) 
h15<-process_stats(15)
cbind(min(h15$posteriors), max(h15$posteriors))
#    [,1]     [,2]
#[1,] 2.335279 9.976224
v15<-ggplot(h15, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[15]) + ylim(0.1, max(h15$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v15<- v15 + geom_segment(aes(x='logit_ABC',y=0.1,xend='logit_ABC',yend=6), position = position_nudge(x = -0.5))

#w_anc_t7    runif(1000,min=10,max=100)
h16<-process_stats(16)
cbind(min(h16$posteriors), max(h16$posteriors))
#   [,1]     [,2]
#[1,] 19.03338 102.1845
v16<-ggplot(h16, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[16]) + ylim(10, max(h16$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v16<- v16 + geom_segment(aes(x='logit_ABC',y=10,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

#t8  runif(1000,min=1.5,max=15)
h17<-process_stats(17)
cbind(min(h17$posteriors), max(h17$posteriors))
#     [,1]     [,2]
#[1,] 4.080417 18.00849
v17<-ggplot(h17, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[17]) + ylim(1.5, max(h17$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v17<- v17 + geom_segment(aes(x='logit_ABC',y=1.5,xend='logit_ABC',yend=15), position = position_nudge(x = -0.5))

#gor_anc runif(1000,min=10,max=100)
h18<-process_stats(18)
cbind(min(h18$posteriors), max(h18$posteriors))
#     [,1]     [,2]
#[1,] -1.279905 123.2393
v18<-ggplot(h18, aes(x=ID, y=posteriors, fill=ID)) +  geom_violin() +  xlab(colnames(window_ABC)[18]) + ylim(min(h18$posteriors), max(h18$posteriors)) + theme_classic() + geom_boxplot(width=0.1)  + scale_fill_brewer(palette = "Set2") + theme(axis.text.x=element_text(size=10)) + theme(axis.text.x = element_blank()) 
fin_v18<- v18 + geom_segment(aes(x='logit_ABC',y=10,xend='logit_ABC',yend=100), position = position_nudge(x = -0.5))

# arrange all in one figure

library(ggpubr)

#for i in {1..18}
#do
#echo "fin_v$i"
#done



#pdf("/Users/harvi/Downloads/gorilla_abc/11sep21simns/segrecalc/posterior.dist.28feb22.replot.pdf") 

#ggarrange(fin_v1,fin_v2,fin_v3,fin_v4,fin_v5,fin_v6,fin_v7,fin_v8,fin_v9,fin_v10,fin_v11,fin_v12,fin_v13,fin_v14,fin_v15,fin_v16,fin_v17,fin_v18, 
 #         ncol = 5, nrow = 4)

#dev.off()


ggarrange(fin_v1,fin_v2,fin_v3,fin_v4,fin_v5,fin_v6,fin_v7,fin_v8,fin_v9,fin_v10,fin_v11,fin_v12,fin_v13,fin_v14,fin_v15,fin_v16,fin_v17,fin_v18, ncol=5, nrow=4, common.legend = TRUE, legend="bottom")

ggarrange(fin_v1,fin_v2,fin_v3,fin_v4,fin_v5,fin_v6,fin_v7,fin_v8,fin_v9,fin_v10,fin_v11,fin_v12,fin_v13,fin_v14,fin_v15,fin_v16,fin_v17,fin_v18, ncol=6, nrow=3, common.legend = TRUE, legend="bottom")


#https://ggplot2-book.org/scale-colour.html

# exported as svg & pdfs by changign the aspect ratio: posterior.comp.pdf

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
