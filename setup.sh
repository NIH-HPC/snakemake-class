#! /bin/bash

function fail() {
    echo "FAIL: $@" >&2
    exit 1
}

function info() {
    echo "INFO: $@"
}

module load singularity || fail "*** Please run the setup script in an sinteractive session ***"
module load snakemake/7 || fail "Could not load snakemake 7 module"
module load git || fail "Could not load git module"

## fetch the biowulf profile

if [[ ! -d biowulf_profile ]]
then
    info "Fetching snakemake profile for Biowulf from https://github.com/NIH-HPC/snakemake_profile"
    git clone https://github.com/NIH-HPC/snakemake_profile.git bwprofile &> /dev/null \
        || fail "unable to clone profile repo"
    echo "use-singularity: true" >> bwprofile/config.yaml
else
    info "Snakemake profile for biowulf already downloaded"
fi
info "The profile has been configured to use singularity"

# running the setup workflow will also cache the latest rnaseq container
export SNAKEMAKE_PROFILE=${PWD}/bwprofile
if snakemake -s setup.smk setup
then
    cat <<-EOF
	+------------------------------------------------------------------------------+
	|                                                                              |
	|                Class materials have been set up successfully                 |
	|                                                                              |
	+------------------------------------------------------------------------------+
	EOF
else
    cat <<-EOF
	+------------------------------------------------------------------------------+
	|                                                                              |
	|                        An error occured during setup                         |
	|                                                                              |
	+------------------------------------------------------------------------------+
	EOF
fi
