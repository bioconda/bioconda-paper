import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import common

plt.figure(figsize=(4,2))

packages = pd.read_table(snakemake.input[0])

sns.countplot(x="ecosystem", data=packages)
plt.ylabel("count")
plt.text(plt.xlim()[1], plt.ylim()[1], "total: {}".format(packages.shape[0]),
         horizontalalignment="right", verticalalignment="top")
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
