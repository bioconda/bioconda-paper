import os
from operator import itemgetter
from itertools import filterfalse

import pandas as pd
from github import Github

authors = pd.read_table(snakemake.input[0], index_col=0)
authors["commits"] = 0
def add(commit):
    if commit.author and commit.author.login in authors.index:
        authors.loc[commit.author.login, "commits"] += 1

# add commits
github = Github(os.environ["GITHUB_TOKEN"])
repo = github.get_repo("bioconda/bioconda-recipes")
utils_repo = github.get_repo("bioconda/bioconda-utils")

for commit in repo.get_commits():
    add(commit)

for commit in utils_repo.get_commits():
    add(commit)

# order by commits
authors.sort_values("commits", inplace=True, ascending=False)

first_authors = ["bgruening", "daler"]
core_authors = ["chapmanb", "jerowe", "tomkinsc", "rvalieris", "druvus"]
last_author = ["johanneskoester"]

# put core into the right order
authors = pd.concat([
    authors.loc[first_authors],
    authors.loc[core_authors].sort_values("commits", ascending=False),
    authors.loc[~authors.index.isin(first_authors + core_authors + last_author)],
    authors.loc[last_author]])

authors.to_csv(snakemake.output.table, sep="\t")
