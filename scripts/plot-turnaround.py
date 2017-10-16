from datetime import timedelta
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import common

plt.figure(figsize=(3,3))
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
plt.pie(counts, shadow=False, labels=counts.index, autopct="%.0f%%")

plt.savefig(snakemake.output[0], bbox_inches="tight")
