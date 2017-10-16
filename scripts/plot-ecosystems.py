import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import os

import common

plt.figure(figsize=(4,2))
summary = pd.read_table(snakemake.input.bio)
not_bio_related = summary['name'][summary.not_bio_related == 'x']

packages = pd.read_table(snakemake.input.pkg_data)
packages.loc[packages.ecosystem == 'Bioconductor', 'ecosystem'] = 'Bioconductor/R'
packages.loc[packages.ecosystem == 'R', 'ecosystem'] = 'Bioconductor/R'
recipes = set(map(os.path.basename, glob.glob('bioconda-recipes/recipes/*')))

packages['has_current_recipe'] = packages['package'].isin(recipes)
packages['not_bio_related'] = packages['package'].isin(not_bio_related)

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(1, 1, 1)
all_cnts = packages.ecosystem.value_counts()
bio_cnts = packages[~packages.not_bio_related].ecosystem.value_counts()
non_cnts = packages[packages.not_bio_related].ecosystem.value_counts()

x = range(len(all_cnts))
ax.bar(x=x, height=bio_cnts, color=sns.color_palette('colorblind'))
with sns.color_palette('Set1', desat=0.7) as c:
    ax.bar(x=x, height=non_cnts, bottom=bio_cnts, color=sns.color_palette(sns.color_palette('colorblind'), desat=0.5))
ax.set_ylabel('count')
ax.set_xticks(x)
ax.set_xticklabels(list(all_cnts.index))

#sns.countplot(x="ecosystem", data=packages)
plt.ylabel("count")
plt.text(plt.xlim()[1], plt.ylim()[1], "total: {}".format(packages.shape[0]),
         horizontalalignment="right", verticalalignment="top")
sns.despine()
plt.gcf().tight_layout()

plt.savefig(snakemake.output[0], bbox_inches="tight")
