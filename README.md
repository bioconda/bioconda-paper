[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1068297.svg)](https://doi.org/10.5281/zenodo.1068297)

# Data analysis for the Bioconda paper

This Snakemake workflow automatically generates all results and figures from the Bioconda paper.

## Requirements

Any 64-bit Linux installation with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).


## Usage

This workflow can be used to recreate all results found in the Bioconda paper.

### Step 1: Install Miniconda3

First, install the Miniconda3 distribution, a version of the conda package manager, together with a minimal Python 3. To reproduce the results of this workflow, do this on a 64bit Linux system with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).
Download Miniconda3 with

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Make sure to answer `yes` to the question whether your PATH variable shall be modified.
Afterwards, open a new shell/terminal.

If you already have a Miniconda installation that defaults to Python 3, you can skip this step.

### Step 2: Setup Bioconda channel

Setup Bioconda with

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

### Step 3: Install bioconda-utils and Snakmake

Install bioconda-utils and Snakemake with

    conda install bioconda-utils snakemake

### Step 4: Download the workflow

First, create a working directory:

    mkdir bioconda-workflow
    cd bioconda-workflow

Then, download the workflow archive from https://doi.org/10.5281/zenodo.1068297 and unpack it with

    tar -xf bioconda-paper-workflow.tar.gz

### Step 5: Run the workflow

Execute the analysis workflow with Snakemake

    snakemake --use-conda

This tells Snakemake to run at most as many jobs in parallel as you have CPU cores.
In addition the `api_requests` resource tells Snakemake to run only one job that
queries the anaconda.org API at a time (we don't want to perform a denial of
service attack :-)).

Results can be found in the folder `plots/`.
