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
fig = plt.figure(figsize=(4,3))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)
ax1.set_ylabel('authors')
ax2.set_ylabel('recipes')
df.plot('time', 'cumulative_authors', ax=ax1, legend=False, color='0.3')
df.plot('time', 'cumulative_recipes', ax=ax2, legend=False, color='0.3')
ax1.set_xlabel('')
ax1.set_xticklabels([])
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
for tick in ax2.get_xticklabels():
    tick.set_rotation(45)
sns.despine()
fig.tight_layout()
fig.savefig(outfile)
