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

#ax1 = fig.add_subplot(1, 1, 1)
#ax2 = fig.add_subplot(2, 1, 2)

#df.semilogy('time', 'cumulative_authors', ax=ax1, legend=True)
#df.semilogy('time', 'cumulative_recipes', ax=ax1, legend=True)
#ax1.set_ylabel('Authors')
#ax2.set_ylabel('Recipes')
#ax2.set_xlabel('Time')
#ax1.set_xlabel('')
#ax1.set_xticklabels([])
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
#plt.gca().set_xticklabels([])
plt.xticks([])
#plt.xticks(rotation=45, ha="right")
sns.despine()

fig.savefig(outfile, bbox_inches="tight")
