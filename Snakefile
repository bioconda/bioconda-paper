import requests


repodata = requests.get("https://conda.anaconda.org/bioconda/linux-64/repodata.json").json()
packages = set(p["name"] for p in repodata["packages"].values())


rule all:
    input:
        expand("plots/{plot}.pdf",
               plot=["downloads", "ecosystems", "downloads_violin", "age-vs-downloads"])



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


rule clone:
    output:
        "bioconda-recipes/.git/index"
    shell:
        "rm -rf bioconda-recipes && "
        "git clone https://github.com/bioconda/bioconda-recipes.git bioconda-recipes"


rule git_log:
    input:
        "bioconda-recipes/.git/index"
    output:
        "git-log/bioconda-recipes.log"
    shell:
        '(cd bioconda-recipes && git pull); '
        '(cd bioconda-recipes && '
        'git log '
        '--pretty=format:'
        '"%h\t%aN\t%aI" '
        '--name-only '

        # parsing the full log can take a while; for debugging try this arg:
        # '--max-count 10 '

        'recipes/*) '
        '> {output}'


rule plot_downloads:
    input:
        "package-data/all.tsv"
    output:
        "plots/downloads.pdf",
        "plots/downloads_violin.pdf",
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


rule parse_git_log:
    input:
        "git-log/bioconda-recipes.log"
    output:
        "git-log/parsed-log.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/parse-log.py"


rule contributions_plot:
    input:
        "git-log/parsed-log.tsv"
    output:
        "plots/contributions.pdf"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-contributions.py"

rule plot_age_vs_downloads:
    input:
        log="git-log/parsed-log.tsv",
        pkg="package-data/all.tsv"
    output:
        "plots/age-vs-downloads.pdf"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-age-vs-downloads.py"
