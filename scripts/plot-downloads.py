import matplotlib.pyplot as plt
import seaborn as sns
import json
import pandas as pd

packages = pd.read_table(snakemake.input[0])

sns.stripplot(x="package",
              y="downloads",
              hue="ecosystem",
              data=packages,
              rasterized=True)
plt.savefig(snakemake.output[0], bbox_inches="tight")
