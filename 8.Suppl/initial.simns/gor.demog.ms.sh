## Tue 26 Jan 2021 11:39:19 CET
## converting msprime gor demog to ms (coalescent simulator with functionality to sample parameter vals from a probability distribution)

#------------------------------------------------------------------------------------------------------------------------
# taking demog from msprime.gordemog.wcem1.py - besenbacher et al. 2019 mutation rate parameters, xue et al, & mcmanus et al. demog parameters
#------------------------------------------------------------------------------------------------------------------------

ms 88 20000 \  ## 88 chr (isu), 20,000 replicates 
\  ## reference N is 1000, generation time is 19years. 40kb locus.
-t tbs \  ## scaled mutation rate 4 * N * u * length - pass from tetaroAr_corr.txt (generated for gor)
-r tbs 40000.0 \  ## scaled recombination rate. R * 4 * N * length of fragment (rho tbs, nsites=40000.0)
-s val2 \  # number of segregating sites - also pass step sizes of number of mutations per region
-I 4 44 2 18 24 \  # sample configuration: 4 populations : 44 wes_lowl (pop1), 2 wes_cros (pop2), 18 eas_lowl (pop3), 24 eas_moun (pop4) initial migration set to 0 (ie sample 44 chr from pop1, etc)
\  # parameters at sampling time - in the present 
-n 1 25.161 \  # wes_lowl pop size = 25.161 *10^3 - (parameters mcmanus 2015, using 12my divergence for humans & gorillas)
-n 2 3.080 \  # wes_cros pop size = 3.080 *10^3 - mcmanus
-n 3 4.280 \  # eas_lowl pop size = 4.280 *10^3 - mcmanus
-n 4 0.800 \  # eas_moun pop size = 800 - xue et al. 2015 
\  # 15kya: 15,000/(4 * 1000 * 19) =  0.1973684
\  # 1. E subspecies split # merge eas_moun into eas_lowl 
# ej t i j : move all lineages in i to j at time t
-ej 0.1973684 4 3 \
\  # 34kya: 34,000/(4 * 1000 * 19) = 0.4473684
\  # 2. migration pulse W > E 3% - for admixture need to give as a proportion not percent** 
\  # # western lowland -> eastern lowland - (30% (0.3) inferred from mcmanus et al 2015 gphocs modelling 
-em 0.4473684 3 1 1200 \  # migration to 3 (eas_lowl) from 1 (wes_lowl) # 0.3*4000 = 1200
-em 0.4473684 3 2 16 \   # migration to 3 (eas_lowl) from 2 (wes_cros) # 0.004*4000 = 16 # western cross river -> eastern lowland 0.4% (0.004) inferred from mcmanus et al 2015 gphocs
\  #The migration pulse should happen for only one generation and then stop:
\  # ie set migrations back to 0
-em 0.4521184 3 1 0 \
-em 0.4521184 3 2 0 \
\  # 68kya: 68,000/(4 * 1000 * 19) = 0.8947368
\  # 3. W subspecies split ie merge wes_cros into wes_lowl
-ej 0.8947368 2 1 \
\  # 4. pop size change for (wes_lowl) - to size of anc_wes (common anc of wes_lowl & wes_cros)
-en 0.8947368 1 30.693 \  # 30693/N 
\  # 5. gor split (W - E) # merge eastern into western at time 261,000/(4 * 1000 * 19) = 3.434211
-ej 3.434211 3 1 \
\  # 6. pop size change to size of anc_gor at time in years/ (4 * N * gen_time)
-en 3.434211 1 39.751

#------------------------------------------------------------------------------------------------------------------------

# 'high divergence' ms simulations - scripts: /scratch/devel/hpawar/admix/sstar/scripts/high.divergence.simul.arr, /scratch/devel/hpawar/admix/sstar/scripts/high.divergence.simul.R
	# replacing lines 37-40 with: 
#\  # 5. gor split (W - E) # merge eastern into western at time 429,000/(4 * 1000 * 19) =  5.644737
-ej 5.644737 3 1 \
#\  # 6. pop size change to size of anc_gor at time in years/ (4 * N * gen_time)
-en 5.644737 1 39.751

#------------------------------------------------------------------------------------------------------------------------

# Mon  1 Mar 2021 15:22:55 CET
# run a set of simulations explicitly stating the E pop size after merging - does not make a difference to the results
# ie after line 22 add:
# specify pop size of E anc to 4.280 *10^3 at 15kya: 15,000/(4 * 1000 * 19) =  0.1973684
-en 0.1973684 3 4.280
