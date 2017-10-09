import os
import pandas as pd
import glob
import csv

packages = pd.read_table(snakemake.input[0])

# restrict to existing recipes
recipes = set(map(os.path.basename, glob.glob('bioconda-recipes/recipes/*')))
packages['has_current_recipe'] = packages['package'].isin(recipes)
packages = packages[packages.has_current_recipe]

with open(snakemake.output[0], "w") as out:
    out = csv.writer(out, delimiter="\t")
    out.writerow(["downloads", packages["downloads"].sum()])
    out.writerow(["versions", packages["versions"].sum()])
    out.writerow(["packages", packages.shape[0]])
