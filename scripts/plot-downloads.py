import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

packages = pd.read_table(snakemake.input[0])

sns.set_style("ticks")
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

# Violin plots to see a little more structure (e.g., 3 tiers of downloads in
# Perl, BioC, R) and lower-limits (e.g., all BioC downloaded at least once, but
# some Perl, Python, R never downloaded).
#
# Take the log10 ahead of time so the KDE works well.
packages['log10(downloads)'] = np.log10(packages.downloads + 1)
fig = plt.figure()
sns.violinplot(
    x="ecosystem",
    y="log10(downloads)",
    alpha=0.5,
    cut=0,
    data=packages)
sns.despine()
fig.savefig(snakemake.output[1], bbox_inches="tight")
