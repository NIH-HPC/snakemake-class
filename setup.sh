#! /bin/bash

function fail() {
    echo "$@" >&2
    exit 1
}

module load singularity || fail "*** Please run the setup script in an sinteractive session ***"
module load snakemake || fail "Could not load snakemake module"

snakemake --use-singularity \
    --singularity-args '-B $PWD:/data' \
    --singularity-prefix=00container setup
if [[ $? -eq 0 ]]; then
    cat <<EOF
+------------------------------------------------------------------------------+
|                                                                              |
|                Class materials have been set up successfully                 |
|                                                                              |
+------------------------------------------------------------------------------+
EOF
else
    cat <<EOF
+------------------------------------------------------------------------------+
|                                                                              |
|                        An error occured during setup                         |
|                                                                              |
+------------------------------------------------------------------------------+
EOF
fi
