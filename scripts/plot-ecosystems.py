import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

packages = pd.read_table(snakemake.input[0])
sns.set_style("white")
sns.set_palette="Set1")

sns.countplot(x="ecosystem", data=packages)
plt.ylabel("count (total: {})".format(packages.shape[0]))
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
