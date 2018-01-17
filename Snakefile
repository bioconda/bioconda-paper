import requests


def get_packages(arch):
    repodata = requests.get("https://conda.anaconda.org/bioconda/{arch}/repodata.json".format(arch=arch)).json()
    packages = set(p["name"] for p in repodata["packages"].values())
    return packages

packages = get_packages("linux-64") | get_packages("osx-64") | get_packages("noarch")
print(len(packages))


rule all:
    input:
        "plots/add+del.pdf",
        expand("figs/fig{f}.pdf", f=[1, 2]),
        expand("plots/{plot}.svg",
               plot=["downloads",
                     "ecosystems",
                     "contributions",
                     "downloads_violin",
                     "age-vs-downloads",
                     "dag"]),
        expand("plots/{pkg}.dag.colored.svg",
               # identified via scripts/find-most-deps.py
               pkg=[
                   'multigps',
                   'cnvkit',
                   'jaffa',
                   'gimmemotifs',
                   'mageck-vispr',
                   'cap-mirseq',
                   'qcumber'])


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
        '(cd bioconda-recipes && '
        'git log '
        '--pretty=format:'
        '"%h\t%aN\t%aI" '
        '--name-only '

        # parsing the full log can take a while; for debugging try this arg:
        # '--max-count 10 '

        'recipes/*) '
        '> {output}'


rule get_dag:
    input:
        "bioconda-recipes/.git/index"
    output:
        "dag/dag.dot"
    shell:
        "cd bioconda-recipes; "
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


rule collect_pr_data:
    output:
        "pr/all.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/collect-pr-data.py"


# This rule purposefully left off the default DAG because this needs
# bioconda_utils + conda-build to parse recipes, which in turn need to be in
# the root environment.
rule collect_summaries:
    input:
        "bioconda-recipes/.git/index"
    output:
        "summary/summaries.tsv"
    script:
        "scripts/collect-summaries-and-urls.py"


################# Plots #################


rule plot_adddel:
    output:
        "plots/add+del.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-add-del.py"


rule plot_package_degree:
    input:
        "package-data/all.tsv"
    output:
        "plots/package_degrees.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-package-degrees.py"


rule plot_dag:
    input:
        "dag/dag.dot"
    output:
        "plots/dag.svg"
    conda:
        "envs/analysis.yaml"
    shell:
        "set +o pipefail; ccomps -zX#0 {input} | neato -Tsvg -o {output} "
        '-Nlabel="" -Nstyle=filled -Nfillcolor="#1f77b4" '
        '-Ecolor="#3333335f" -Nwidth=0.2 -LC10 -Gsize="12,12" '
        "-Nshape=circle -Npenwidth=0"


rule plot_colored_dag:
    input:
        pkg='package-data/all.tsv',
        dag='dag/dag.dot'
    output:
        'plots/{pkg}.dag.colored.svg'
    conda:
        'envs/analysis.yaml'
    script:
        'scripts/color-dag.py'


# rule plot_colored_dag:
#     input:
#         'dag/{pkg}.colored.dot'
#     output:
#         'plots/{pkg}.dag.colored.svg'
#     conda:
#         'envs/analysis.yaml'
#     shell:
#         "set +o pipefail; ccomps -zX#0 {input} | neato -Tsvg -o {output} "
#         '-Nlabel="" -Nstyle=filled -Nfillcolor="#7777778f" '
#         '-Ecolor="#3333335f" -Nwidth=0.2 -LC10 -Gsize="12,12" '
#         '-Earrowhead="none" -Nshape=circle -Npenwidth=0 '
#         '-Goutputorder=edgesfirst'


rule plot_downloads:
    input:
        "package-data/all.tsv",
    output:
        "plots/downloads.svg",
        "plots/downloads_violin.svg",
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-downloads.py"


rule plot_ecosystems:
    input:
        "bioconda-recipes/.git/index",
        pkg_data="package-data/all.tsv",
        bio="summary/hand-edited-summaries.tsv"
    output:
        "plots/ecosystems.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-ecosystems.py"


rule plot_comparison:
    input:
        "summary/hand-edited-summaries.tsv"
    output:
        counts="plots/pkg-count-comparison.svg",
        age="plots/age-comparison.svg",
        csv="plots/comparison.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-comparison.py"


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


rule plot_turnaround:
    input:
        "pr/all.tsv"
    output:
        "plots/turnaround.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/plot-turnaround.py"

########### Tables ##############


rule stats:
    input:
        "bioconda-recipes/.git/index",
        pkg="package-data/all.tsv"
    output:
        "package-data/stats.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/stats.py"


rule author_list:
    input:
        "resources/authors.tsv"
    output:
        tex="resources/authors.tex",
        table="resources/authors-commits.tsv"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/author-list.py"


########### Figures #############

rule fig1:
    input:
        comp="plots/pkg-count-comparison.svg",
        age="plots/age-comparison.svg",
        downloads="plots/downloads_violin.svg",
        ecosystems="plots/ecosystems.svg"
    output:
        "figs/fig1.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/fig1.py"


rule fig2:
    input:
        add_del="plots/add+del.svg",
        contributions="plots/contributions.svg",
        workflow="plots/workflow.svg",
        dag="plots/cnvkit.dag.colored.svg",
        turnaround="plots/turnaround.svg"
    output:
        "figs/fig2.svg"
    conda:
        "envs/analysis.yaml"
    script:
        "scripts/fig2.py"


rule convert_svg:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.{fmt,(pdf|png)}"
    conda:
        "envs/cairosvg.yaml"
    shell:
        "cairosvg -f {wildcards.fmt} {input} -o {output}"
