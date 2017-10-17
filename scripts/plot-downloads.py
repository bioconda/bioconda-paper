import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import common

plt.figure(figsize=(4,2))
packages = pd.read_table(snakemake.input[0])
total_downloads = packages["downloads"].sum()
packages.loc[packages.ecosystem == 'Bioconductor', 'ecosystem'] = 'Bioconductor/R'
packages.loc[packages.ecosystem == 'R', 'ecosystem'] = 'Bioconductor/R'

# In case we want to filter downloads by whether or not a current recipe exists
recipes = set(map(os.path.basename, glob.glob('bioconda-recipes/recipes/*')))

sns.boxplot(x="ecosystem",
            y="downloads",
            data=packages,
            color="white",
            whis=False,
            showfliers=False,
            order=['Bioconductor/R', 'Other', 'Python', 'Perl'],
           )
sns.stripplot(x="ecosystem",
              y="downloads",
              data=packages,
              jitter=True,
              alpha=0.5,
              order=['Bioconductor/R', 'Other', 'Python', 'Perl'],
             )
plt.gca().set_yscale("log")
plt.ylabel("downloads (total: {:,})".format(total_downloads))
sns.despine()


plt.savefig(snakemake.output[0], bbox_inches="tight")

# Violin plots to see a little more structure (e.g., 3 tiers of downloads in
# Perl, BioC, R) and lower-limits (e.g., all BioC downloaded at least once, but
# some Perl, Python, R never downloaded).
#
# Take the log10 ahead of time so the KDE works well.
packages['log10 downloads'] = np.log10(packages.downloads + 1)

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)
sns.violinplot(
    x="ecosystem",
    y="log10 downloads",
    alpha=0.5,
    cut=0,
    data=packages,
    ax=ax,
    order=['Bioconductor/R', 'Other', 'Python', 'Perl'],
)
ax.text(x=0.5, y=1.0, s="Total downloads: {:,}".format(total_downloads),
         horizontalalignment="center", verticalalignment="top",
        transform=ax.transAxes)
ax.set_xlabel('')

# make a little room for the "total" text
ax.axis(ymax=6)
fig.tight_layout()
sns.despine()
plt.savefig(snakemake.output[1], bbox_inches="tight")
