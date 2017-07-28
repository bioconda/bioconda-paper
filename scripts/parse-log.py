"""
Parses and plots the output of::

    git log --pretty=format:"%h\t%aN\t%aI" --name-only recipes/*

which looks something like::

    <hash>    <user>   <timestamp>
    recipes/<package>/build.sh
    recipes/<package>/meta.yaml


    <hash>    <user>   <timestamp>
    recipes/<package1>/build.sh
    recipes/<package1>/meta.yaml
    recipes/<package2>/build.sh
    recipes/<package2>/meta.yaml

"""
import pandas
import os
import datetime

infile = snakemake.input[0]
outfile = snakemake.output[0]


class chunk(object):
    def __init__(self, block):
        commit, author, time = block[0].split('\t')
        self.author = author
        self.time = datetime.datetime.strptime(time.split('T')[0], "%Y-%m-%d")
        self._block = block
        self.recipes = self._parse_recipes(block[1:])

    def _parse_recipes(self, block):
        recipes = []
        for i in block:
            if not i.startswith('recipes/'):
                continue
            if os.path.basename(i) != 'meta.yaml':
                continue
            recipes.append(os.path.dirname(i.replace('recipes/', '')))
        return set(recipes)


def gen():
    lines = []
    for line in open(infile):
        line = line.strip()
        if len(line) == 0:
            yield chunk(lines)
            lines = []
            continue
        lines.append(line.strip())
    yield chunk(lines)


dfs = []
cumulative_recipes = set()
cumulative_authors = set()
for i in sorted(gen(), key=lambda x: x.time):
    if len(i.recipes) == 0:
        continue

    unique_recipes = i.recipes.difference(cumulative_recipes)
    if len(unique_recipes) > 0:
        dfs.append(
            {
                'time': i.time,
                'author': i.author,
                'recipes': unique_recipes,
                'nadded': len(unique_recipes),
                'new_author': i.author not in cumulative_authors
            },
        )
    cumulative_recipes.update(i.recipes)
    cumulative_authors.update([i.author])

df = pandas.DataFrame(dfs)
df['cumulative_authors'] = df.new_author.astype(int).cumsum()
df['cumulative_recipes'] = df.nadded.cumsum()
df["time"] = pandas.to_datetime(df["time"])
df.to_csv(outfile, sep='\t')
