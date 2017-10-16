from datetime import timedelta
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import common

sns.set_palette("Blues")

prs = pd.read_table(snakemake.input[0])
prs.span = pd.to_timedelta(prs.span)

categories = pd.Series([timedelta(minutes=0), timedelta(minutes=30),
                        timedelta(hours=1), timedelta(hours=5),
                        timedelta(days=1),
                        timedelta(days=365)])
labels = [r"$\leq 30$ min", r"$\leq 1$ hour",
          r"$\leq 5$ hours", r"$\leq 1$ day", r"$ > 1$ day"]
binning = pd.cut(prs.span.dt.total_seconds(),
                 categories.dt.total_seconds(),
                 labels=labels)
counts = binning.value_counts()
# fix order
counts = counts[labels]

perc = counts / counts.sum()
fig = plt.figure(figsize=(5, 1.3))
ax = fig.add_subplot(1, 1, 1)
left = 0
for label, x in perc.items():
    ax.barh(y=0, width=x, left=left, label=label)
    ax.text(left + x/2, 0, label, horizontalalignment='center', verticalalignment='center')
    left += x

sns.despine(top=True, left=True, right=True, trim=True)
ax.set_xlabel('Fraction of pull requests')
ax.yaxis.set_visible(False)
fig.tight_layout()
fig.subplots_adjust(top=0.9)

#plt.pie(counts, shadow=False, labels=counts.index, autopct="%.0f%%")
plt.savefig(snakemake.output[0], bbox_inches="tight")
