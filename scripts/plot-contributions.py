import seaborn as sns
import os
from matplotlib import pyplot as plt
import pandas as pd


df = pd.read_table(snakemake.input[0])

df[["time", "cumulative_authors", "cumulative_recipes"]].plot(x="time",
                                                              subplots=True,
                                                              sharex=True,
                                                              layout=(2, 1))
sns.despine()
plt.savefig(snakemake.output[0], bbox_inches="tight")
