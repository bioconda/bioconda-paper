import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

packages = pd.read_table(snakemake.input[0])

sns.set_style("white")
sns.set_palette("Set1")

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
plt.ylabel("downloads (total: {})".format(packages["downloads"].sum()))
sns.despine()

plt.savefig(snakemake.output[0], bbox_inches="tight")
