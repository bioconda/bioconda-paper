import seaborn as sns
import pandas
import os
from matplotlib import pyplot as plt
import datetime

sns.set_style('white')

infile = snakemake.input[0]
outfile = snakemake.output[0]

df = pandas.read_table(infile)
fig = plt.figure(figsize=(6, 8))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)
df.plot('time', 'cumulative_authors', ax=ax1)
df.plot('time', 'cumulative_recipes', ax=ax2)
ax1.set_ylabel('Number of unique authors')
ax2.set_ylabel('Number of unique recipes')
fig.tight_layout()
fig.savefig(outfile)
