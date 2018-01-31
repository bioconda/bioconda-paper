"""
Get number of direct dependencies of each package and plot histogram

"""
import glob
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import json
import pandas as pd

import common

plt.figure(figsize=(4,2))
packages = pd.read_table(snakemake.input[0])

deps = packages["deps"]


plt.hist(deps, range(0,30), lw=1)
plt.xlim([0,30])
plt.grid()
plt.xlabel("Package degree", fontsize=16)


plt.savefig(snakemake.output[0], bbox_inches="tight")
