import requests


repodata = requests.get("https://conda.anaconda.org/bioconda/linux-64/repodata.json").json()
packages = set(p["name"] for p in repodata["packages"].values())


rule all:
    input:
        expand("plots/{plot}.svg",
               plot=["downloads",
                     "ecosystems",
                     "contributions",
                     "downloads_violin",
                     "age-vs-downloads",
                     "dag"])


############# Collect data ##############


rule get_package_data:
    output:
        "package-data/{package}.json"
    conda:
        "envs/analysis.yaml"
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


rule collect_contributions:
    input:
        "git/bioconda-recipes.log"
    output:
        "git/contributions.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/collect-contributions.py"


rule get_dag:
    input:
        "bioconda-recipes/.git/index"
    output:
        "dag/dag.dot"
    shell:
        "cd bioconda-recipes; "
        "git pull; "
        "bioconda-utils dag --hide-singletons --format dot "
        "recipes config.yml > ../{output}"


rule parse_git_log:
    input:
        "git-log/bioconda-recipes.log"
    output:
        "git-log/parsed-log.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/parse-log.py"


################# Plots #################


rule plot_dag:
    input:
        "dag/dag.dot"
    output:
        "plots/dag.svg"
    conda:
        "envs/analysis.yaml"
    shell:
        "twopi -Tsvg -o {output} "
        '-Nlabel="" -Nstyle=filled -Nfillcolor="#3333335f" '
        '-Ecolor="#3333335f" '
        "-Nshape=circle -Npenwidth=0 {input}"


rule plot_downloads:
    input:
        "package-data/all.tsv"
    output:
        "plots/downloads.svg",
        "plots/downloads_violin.svg",
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-downloads.py"


rule plot_ecosystems:
    input:
        "package-data/all.tsv"
    output:
        "plots/ecosystems.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-ecosystems.py"


rule plot_contributions:
    input:
        "git-log/parsed-log.tsv"
    output:
        "plots/contributions.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-contributions.py"


rule plot_age_vs_downloads:
    input:
        log="git-log/parsed-log.tsv",
        pkg="package-data/all.tsv"
    output:
        "plots/age-vs-downloads.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-age-vs-downloads.py"


########### Figures #############
