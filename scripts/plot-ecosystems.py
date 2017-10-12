import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import os

import common

plt.figure(figsize=(4,2))

packages = pd.read_table(snakemake.input.pkg_data)
packages.loc[packages.ecosystem == 'Bioconductor', 'ecosystem'] = 'Bioconductor/R'
packages.loc[packages.ecosystem == 'R', 'ecosystem'] = 'Bioconductor/R'
recipes = set(map(os.path.basename, glob.glob('bioconda-recipes/recipes/*')))

packages['has_current_recipe'] = packages['package'].isin(recipes)
packages = packages[packages.has_current_recipe]

sns.countplot(x="ecosystem", data=packages)
plt.ylabel("count")
plt.text(plt.xlim()[1], plt.ylim()[1], "total: {}".format(packages.shape[0]),
         horizontalalignment="right", verticalalignment="top")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
