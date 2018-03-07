[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1068297.svg)](https://doi.org/10.5281/zenodo.1068297)
[![Snakemake](https://img.shields.io/badge/snakemake-≥4.6.0-brightgreen.svg)](https://snakemake.bitbucket.io)

# Data analysis for the Bioconda paper

This Snakemake workflow automatically generates all results and figures from the Bioconda paper.

## Requirements

Any 64-bit Linux installation with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6).
Note that the restriction of this workflow to Linux is purely a design decision (to save space and ensure reproducibility) and not related to Conda/Bioconda. Bioconda packages are available for both Linux and MacOS in general.


## Usage

This workflow can be used to recreate all results found in the Bioconda paper.

### Step 1: Setup system

#### Variant a: Installing Miniconda on your system

If you are on a Linux system with [GLIBC 2.5](http://unix.stackexchange.com/a/120381) or newer (i.e. any Linux distribution that is newer than CentOS 6), you can simply install Miniconda3 with

    curl -o /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/miniconda.sh

Make sure to answer `yes` to the question whether your PATH variable shall be modified.
Afterwards, open a new shell/terminal.

#### Variant b: Use a Docker container

Otherwise, e.g., on MacOS or if you don't want to modify your system setup, install [Docker](https://www.docker.com/), run

    docker run -it continuumio/miniconda3 /bin/bash
  
and execute all the following steps within that container.

#### Variant c: Use an existing Miniconda installation

If you want to use an existing Miniconda installation, please be aware that this is only possible if it uses Python 3 by default. You can check this via
  
    python --version

Further, ensure it is up to date with

    conda update --all

### Step 2: Setup Bioconda channel

Setup Bioconda with

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

### Step 3: Install bioconda-utils and Snakmake

Install bioconda-utils and Snakemake >=4.6.0 with

    conda install bioconda-utils snakemake

If you already have an older version of Snakemake, please make sure it is updated to >=4.6.0.

### Step 4: Download the workflow

First, create a working directory:

    mkdir bioconda-workflow
    cd bioconda-workflow

Then, download the workflow archive from https://doi.org/10.5281/zenodo.1068297 and unpack it with

    tar -xf bioconda-paper-workflow.tar.gz

### Step 5: Run the workflow

Execute the analysis workflow with Snakemake

    snakemake --use-conda

Please wait a few minutes for the analysis to finish.
Results can be found in the folder `figs/`.
If you have been running the workflow in the docker container (see above), 
you can obtain the results with

    docker cp <container-id>:/bioconda-workflow/figs .

whith `<container-id>` being the ID of the container.


## Known errors

* If you see an error like
  ```
  ImportError: No module named 'appdirs'
  ```
  when starting Snakemake, you are likely suffering from a bug in an older conda version. Make sure to update your conda installation with 

      conda update --all

  and then reinstall the `appdirs` and `snakemake` package with

      conda install -f appdirs snakemake
* If you see an error like
  ```
  ImportError: Missing required dependencies ['numpy']
  ```
  you are likely suffering from a bug in an older conda version. Make sure to update your conda installation with
  
      conda update --all
  
  and then reinstall the `snakemake` package with

      conda install -f snakemake
