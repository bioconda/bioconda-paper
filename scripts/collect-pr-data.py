import os
import pandas as pd
from github import Github

github = Github(os.environ["GITHUB_TOKEN"])

repo = github.get_repo("bioconda/bioconda-recipes")

prs = []
titles = []
files = []
spans = []
for pr in repo.get_pulls(state="closed"):
    print(pr)
    if pr.merged:
        prs.append(pr.id)
        titles.append(pr.title)
        files.append(pr.changed_files)
        spans.append(pr.merged_at - pr.created_at)

prs = pd.DataFrame({
    "id": prs,
    "title": titles,
    "changed_files": files,
    "span": spans
})

prs.to_csv(snakemake.output[0], sep="\t", index=False)
