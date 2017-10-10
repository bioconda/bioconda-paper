import os
import pandas
from bioconda_utils import utils

repo_dir = os.path.dirname(os.path.dirname(snakemake.input[0]))
recipes = list(utils.get_recipes(os.path.join(repo_dir, 'recipes')))
config = os.path.join(repo_dir, 'config.yml')
df = []
for r in recipes:
    meta = next(utils.load_all_meta(r,config))
    d = dict(
        not_bio_related=" ",
        summary=meta.get('about', {}).get('summary', "").replace('\n', ''),
        name=meta['package']['name'],
        url=meta.get('about', {}).get('home', ""),
    )
    df.append(d)
df = pandas.DataFrame(df).drop_duplicates('name')
df = df.sort_values('name')
df = df[['not_bio_related', 'name', 'summary', 'url']]
df.to_csv(snakemake.output[0], sep='\t', index=False)
