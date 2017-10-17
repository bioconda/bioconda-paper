"""
Inspection script to be run interactively.

Reads the DAG of non-singleton recipes to find those with large sets of
dependencies.

Dataframes `have_all_depedencies` and `counts` are useful for inspecting.
Identified packages can be added to the list of those for which colored dags
will be created.
"""

from collections import Counter
import networkx as nx
from networkx.drawing.nx_agraph import read_dot
import pandas as pd
packages = pd.read_table('package-data/all.tsv', index_col=0)
packages.loc[packages.ecosystem == 'Bioconductor', 'ecosystem'] = 'Bioconductor/R'
packages.loc[packages.ecosystem == 'R', 'ecosystem'] = 'Bioconductor/R'
lookup = packages['ecosystem'].to_dict()
g = read_dot('dag/dag.dot')

for node in g:
    try:
        nx.set_node_attributes(g, node, lookup[node])
    except KeyError:
        nx.set_node_attributes(g, node, 'Other')

d = {}
for i in g:
    deps = set()
    for k, v in nx.dfs_predecessors(g.reverse(), i).items():
        deps.update([k, v])
    d[i] = deps

counts = []
for k, v in d.items():
    ecosystems = []
    for i in v:
        try:
            ecosystems.append(lookup[i])
        except KeyError:
            ecosystems.append('Other')
    c = dict(Counter(ecosystems))
    c['package'] = k
    counts.append(c)
counts = pd.DataFrame(counts).set_index('package').fillna(0)

have_all_ecosystems = counts[(counts>0).sum(axis=1) == 4]
