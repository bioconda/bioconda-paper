# Data analysis for the Bioconda paper

This Snakemake workflow automatically generates all results and figures from the Bioconda paper.

## Requirements

Any 64-bit Linux installation with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).


## Usage

This workflow can be used to recreate all results found in the Bioconda paper.

### Step 1: Setup Bioconda

Setup Bioconda as shown [here](https://bioconda.github.io).

### Step 2: Install Snakemake

Install Snakemake into an isolated environment

    conda create -n snakemake snakemake

### Step 3: Download the workflow

First, create a working directory:

    mkdir bioconda-workflow
    cd bioconda-workflow

Then, download the workflow archive from https://doi.org/10.5281/zenodo.1068297 and unpack it with

    tar -xf bioconda-paper-workflow.tar.gz

### Step 4: Run the workflow

Execute the analysis workflow with Snakemake

    snakemake --use-conda -j --resources api_requests=1

This tells Snakemake to run at most as many jobs in parallel as you have CPU cores.
In addition the `api_requests` resource tells Snakemake to run only one job that
queries the anaconda.org API at a time (we don't want to perform a denial of
service attack :-)).

Results can be found in the folder `plots/`.
