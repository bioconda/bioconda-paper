# Data analysis related to the Bioconda paper

## Usage

This workflow can be used to recreate all results found in the Bioconda paper.

### Step 1:

Setup Bioconda as shown [here](https://bioconda.github.io).

### Step 2:

Install Snakemake into an isolated environment

    conda create -n snakemake snakemake

### Step 3:

Download the workflow archive from https://doi.org/10.5281/zenodo.106829 and unpack it. Then, enter the resulting directory.

### Step 4:

Execute the analysis workflow with Snakemake

    snakemake --use-conda -j --resources api_requests=1

This tells Snakemake to run at most as many jobs in parallel as you have CPU cores.
In addition the `api_requests` resource tells Snakemake to run only one job that
queries the anaconda.org API at a time (we don't want to perform a denial of
service attack :-)).

Results can be found in the folder `plots/`.
