#!/bin/bash

module load Anaconda3/5.3.0
source activate /project2/nobrega/grace/conda/HiC-env
module load python/2.7.12

snakemake \
    --snakefile /project2/nobrega/grace/metabolic_traits_GWAS/scripts/snake/Snakefile \
    -kp \
    -j 500 \
    --use-conda \
    --rerun-incomplete \
    --resources load=100 \
    --cluster-config /project2/nobrega/grace/metabolic_traits_GWAS/scripts/snake/cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --time={cluster.time} \
        --tasks-per-node=1 \
        --partition=broadwl \
        --job-name={cluster.name} \
        --output={cluster.logfile}" \
    $*
