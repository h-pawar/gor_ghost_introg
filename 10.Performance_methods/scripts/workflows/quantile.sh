#!/bin/bash
#SBATCH -J sstar_quantile
#SBATCH -N 16

sstar quantile --model config/simulation/models/null_model.yaml--ms-dir ext/msdir --N0 1000 --nsamp 22 --nreps 20000 --ref-index 4 --ref-size 20 --tgt-index 3 --tgt-size 2 --mut-rate 1.29e-8 --rec-rate 1e-8 --seq-len 40000 --snp-num-range 25 705 5 --output-dir results/inference/sstar/GorillaGhost/quantiles --thread 16
