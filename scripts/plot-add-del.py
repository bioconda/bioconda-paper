import os
import matplotlib.pyplot as plt
import seaborn as sns
from github import Github
import matplotlib.dates as mdates

import common

github = Github(os.environ["GITHUB_TOKEN"])

repo = github.get_repo("bioconda/bioconda-recipes")

weeks = []
additions = []
deletions = []
print(repo.get_stats_participation().all)
for freq in repo.get_stats_code_frequency():
    weeks.append(freq.week)
    additions.append(freq.additions)
    deletions.append(abs(freq.deletions))


plt.figure(figsize=(4,1.2))

plt.semilogy(weeks, additions, "-", label="additions")
plt.semilogy(weeks, deletions, "-", label="deletions")
plt.ylabel("count per week")
plt.legend(bbox_to_anchor=(0.68, 0.65))

plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.xticks(rotation=45, ha="right")

sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
