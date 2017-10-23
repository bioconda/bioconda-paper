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

first_authors = ["daler", "bgruening"]
core_authors = ["chapmanb", "jerowe", "tomkinsc", "rvalieris", "druvus"]
last_author = ["johanneskoester"]

# put core into the right order
authors = pd.concat([
    authors.loc[first_authors],
    authors.loc[core_authors].sort_values("commits", ascending=False),
    authors.loc[~authors.index.isin(first_authors + core_authors + last_author)],
    authors.loc[last_author]])

authors.to_csv(snakemake.output.table, sep="\t")

def parse_affiliations(a):
    return (b.strip() for b in a.split(";"))

def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

affiliations = {a: i + 1
                for i, a in enumerate(unique_everseen(
                    b for a in authors["affiliation"]
                      for b in parse_affiliations(a)))}

print(*sorted(affiliations), sep="\n")

with open(snakemake.output.tex, "w") as tex:
    for i, author in enumerate(authors.itertuples()):
        a = ",".join(str(affiliations[a]) for a in parse_affiliations(author.affiliation))
        footnote = ""
        if i == len(authors) -1:
            footnote = r"\thanks{To whom correspondence should be addressed.}\textsuperscript{,}"
        elif i == 1:
            footnote = r"\thanks{Co-first author}\textsuperscript{,}"
        print(r"\author[{}]{{{}{}}}".format(a, author.name, footnote), file=tex)

    for affiliation, i in sorted(affiliations.items(), key=itemgetter(1)):
        print(r"\affil[{}]{{{}}}".format(i, affiliation), file=tex)
