import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import common

plt.figure(figsize=(4,2))
packages = pd.read_table(snakemake.input.tsv)
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
            showfliers=False)
sns.stripplot(x="ecosystem",
              y="downloads",
              data=packages,
              jitter=True,
              alpha=0.5)
plt.gca().set_yscale("log")
plt.ylabel("downloads (total: {})".format(total_downloads))
sns.despine()


plt.savefig(snakemake.output[0], bbox_inches="tight")

# Violin plots to see a little more structure (e.g., 3 tiers of downloads in
# Perl, BioC, R) and lower-limits (e.g., all BioC downloaded at least once, but
# some Perl, Python, R never downloaded).
#
# Take the log10 ahead of time so the KDE works well.
packages['log10 downloads'] = np.log10(packages.downloads + 1)

plt.figure(figsize=(4,2))
sns.violinplot(
    x="ecosystem",
    y="log10 downloads",
    alpha=0.5,
    cut=0,
    data=packages)
plt.text(plt.xlim()[1], plt.ylim()[1], "total: {}".format(total_downloads),
         horizontalalignment="right", verticalalignment="top")
sns.despine()
plt.savefig(snakemake.output[1], bbox_inches="tight")
