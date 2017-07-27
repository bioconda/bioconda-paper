import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

packages = pd.read_table(snakemake.input[0])

sns.countplot(x="ecosystem", data=packages)
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
