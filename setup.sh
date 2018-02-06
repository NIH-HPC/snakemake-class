#! /bin/bash

function fail() {
    echo "$@" >&2
    exit 1
}

module load singularity || fail "Please run the setup script in an sinteractive session"
module load snakemake || fail "Could not load snakemake module"

snakemake --use-singularity \          
    --singularity-args '-B $PWD:/data --pwd /data' \
    --singularity-prefix=00container setup          

