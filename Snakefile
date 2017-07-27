

import requests


repodata = requests.get("https://conda.anaconda.org/bioconda/linux-64/repodata.json").json()
packages = set(p["name"] for p in repodata["packages"].values())


rule all:
    input:
        expand("plots/{plot}.pdf",
               plot=["downloads", "ecosystems"])



rule get_package_data:
    output:
        "package-data/{package}.json"
    resources:
        api_requests=1
    shell:
        "curl -X GET --header 'Accept: application/json' "
        "https://api.anaconda.org/package/bioconda/{wildcards.package} "
        "> {output} && sleep 1"


rule collect_package_data:
    input:
        expand("package-data/{package}.json", package=packages)
    output:
        "package-data/all.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/collect-pkg-data.py"


rule plot_downloads:
    input:
        "package-data/all.tsv"
    output:
        "plots/downloads.pdf"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-downloads.py"


rule plot_ecosystems:
    input:
        "package-data/all.tsv"
    output:
        "plots/ecosystems.pdf"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-ecosystems.py"
