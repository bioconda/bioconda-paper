import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from datetime import datetime
import matplotlib.dates as mdates
import numpy as np

import common


summary = pd.read_table(snakemake.input[0])
bio_related = summary.shape[0] - (summary["not_bio_related"] == "x").sum()
# Counts from October 2017
data = pd.DataFrame.from_dict({
    "Bioconda": [bio_related, "2015-09"],
    "Debian Med": [882, "2002-05"],
    "Gentoo Science": [480, "2005-10"],   # category sci-biology
    "EasyBuild": [371, "2012-03"], # moduleclass bio
    "Biolinux": [308, "2006"],
    "Homebrew Science": [297, "2009-10"], # tag bioinformatics
    "GNU Guix": [254, "2014-12"],         # category bioinformatics
    "BioBuilds": [118, "2015-11"]}, orient="index").reset_index()
data.columns = ["source", "count", "date"]
data["date"] = pd.to_datetime(data["date"])
# age in years
data["age"] = pd.to_timedelta(datetime.now() - data["date"]).astype('timedelta64[M]') / 12

plt.figure(figsize=(4,1))

sns.barplot(x="source", y="count", data=data)
plt.gca().set_xticklabels([])
plt.xlabel("")
plt.ylabel("Number of explicitly\nbio-related packages")

# set maximum tick to be that of bioconda
yticks = plt.gca().get_yticks()
yticks[-1] = bio_related
plt.gca().set_yticks(yticks)

sns.despine()
plt.savefig(snakemake.output.counts, bbox_inches="tight")

plt.figure(figsize=(4,1))

sns.barplot(x="source", y="age", data=data)
plt.xlabel("")
plt.ylabel("\nage in years")
plt.xticks(rotation=45, ha="right")
#plt.gca().yaxis.set_major_formatter(mdates.AutoDateFormatter(mdates.AutoDateLocator()))


sns.despine()
plt.savefig(snakemake.output.age, bbox_inches="tight")

# store results as csv
data[["source", "count", "age"]].to_csv(snakemake.output.csv, sep="\t", index=False)
