import seaborn as sns
import os
from matplotlib import pyplot as plt
import datetime
import pandas as pd
import matplotlib.dates as mdates

import common

infile = snakemake.input[0]
outfile = snakemake.output[0]

df = pd.read_table(infile)
df["time"] = pd.to_datetime(df["time"])
fig = plt.figure(figsize=(4,1))
plt.semilogy('time', 'cumulative_authors', data=df, label="contributors")
plt.semilogy('time', 'cumulative_recipes', data=df, label="recipes")
plt.legend()
plt.ylabel("count")
plt.xlabel("")

# deactivate xticks because we have them in the plot below in the figure
plt.xticks([])
sns.despine()

fig.savefig(outfile, bbox_inches="tight")
