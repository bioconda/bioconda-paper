[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1068297.svg)](https://doi.org/10.5281/zenodo.1068297)

# Data analysis for the Bioconda paper

This Snakemake workflow automatically generates all results and figures from the Bioconda paper.

## Requirements

Any 64-bit Linux installation with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).


## Usage

This workflow can be used to recreate all results found in the Bioconda paper.

### Step 1: Setup system

* If you are on a Linux system with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6), you can simply install Miniconda3 with

      curl -o /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/miniconda.sh

   Make sure to answer `yes` to the question whether your PATH variable shall be modified.
   Afterwards, open a new shell/terminal.

* Otherwise, e.g., on MacOS or if you don't want to modify your system setup, install [Docker](https://www.docker.com/), run

      docker run -it continuumio/miniconda3 /bin/bash
  
  and execute all the following steps within that container.

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

Results can be found in the folder `plots/`.
