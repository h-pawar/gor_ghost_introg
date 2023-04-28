# Mon 25 Jul 2022 11:01:46 CEST
# check pvals for proportion of overlapping bp: comparison b/n empirical & random reps per scenario

## module load gcc/6.3.0 openssl/1.0.2q R/4.0.1
load("/scratch/devel/hpawar/admix/overlap.s*.skov/overlap.random.regions.pval.20jul22",verbose=T) # check this **
# generated using scripts -   /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/overlap.random.regions.20jul22.R
# called by - /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/overlap.random.regions.4may22.arr


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# run through relevant functions of R script interactively -
    # /scratch/devel/hpawar/admix/sstar/scripts/simul_postabc/downstream/overlap.random.regions.20jul22.R

# empirical vals for the mean overlap

recalc_empirical_mean<-function(scen,skov_version){
overlap40kb_gbb<-list()
pt_estimate<-list()
for (ind in (1:length(scen))) {
overlap40kb_gbb[[ind]]<-intersect_function(scen[[ind]],skov_version)[[2]]
pt_estimate[[ind]]<-mean(overlap40kb_gbb[[ind]][,1])}
return(pt_estimate)
}


emp_mean_all_comp<-cbind(recalc_empirical_mean(scen_GBB,autosomes_sk_gbb),
recalc_empirical_mean(scen_GBB,autosomes_40sk_gbb),
recalc_empirical_mean(scen_GBG,autosomes_sk_gbg),
recalc_empirical_mean(scen_GBG,autosomes_40sk_gbg))

#-----------------------------------------------------------------------------------------------------------------------

calc_p_fun<-function(scen,sstar_ver){
out <- length( random_p[[scen]][[sstar_ver]][[1]][ which(random_p[[scen]][[sstar_ver]][[1]] >= emp_mean_all_comp[,scen][sstar_ver]), ])/100
return(out)
}


# for all scernarios (gbb, gbg, gbb.40, gbg.40)

hold_j<-c()
for (j in (1:ncol(emp_mean_all_comp))) { 
per_ci<-c()
for (i in (1:nrow(emp_mean_all_comp))) { 
per_ci[[i]]<-calc_p_fun(j,i)
}
ids<-do.call(rbind,per_ci)
hold_j[[j]]<-ids
}

