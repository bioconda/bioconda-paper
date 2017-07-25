import os

if os.path.exists('bioconda-recipes'):
    shell('cd bioconda-recipes && git pull origin master')

rule targets:
    input: 'plots/contributions.pdf'

rule clone:
    output: 'bioconda-recipes/.git/index'
    shell:
        'rm -rf bioconda-recipes && '
        'git clone https://github.com/bioconda/bioconda-recipes.git bioconda-recipes'

rule git_log:
    input: 'bioconda-recipes/.git/index'
    output: 'git.log'
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

rule contributions_plot:
    input: rules.git_log.output
    output: 'plots/contributions.pdf'
    script: 'scripts/plot-contributions.py'
# vim: ft=python
