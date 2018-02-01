from snakemake.shell import shell
import matplotlib
matplotlib.use("agg")
import pandas as pd
import seaborn as sns
import networkx as nx
import glob
import os
from networkx.drawing.nx_agraph import read_dot, write_dot, graphviz_layout
from matplotlib.colors import rgb2hex
from pydotplus import parser
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

packages = pd.read_table(snakemake.input.pkg, index_col=0)
packages.loc[packages.ecosystem == 'Bioconductor', 'ecosystem'] = 'Bioconductor/R'
packages.loc[packages.ecosystem == 'R', 'ecosystem'] = 'Bioconductor/R'
lookup = packages['ecosystem'].to_dict()
colors = dict(zip(['Bioconductor/R', 'Other', 'Python', 'Perl'], sns.color_palette('colorblind')))
g = read_dot(snakemake.input.dag)
# reduce to largest connected component
g = max(nx.weakly_connected_component_subgraphs(g), key=len)

pkg = snakemake.wildcards.pkg
# obtain dependencies
deps = set(nx.ancestors(g, pkg))
sub = deps | {pkg}

pos = graphviz_layout(g, prog='neato')

plt.figure(figsize=(6,6))
# draw DAG
nx.draw_networkx_edges(g, pos, edge_color='#777777', alpha=0.5, arrows=False)
nx.draw_networkx_nodes(g, pos, node_color='#333333', alpha=0.5, node_size=6)

# draw induced subdag
nx.draw_networkx_edges(g, pos,
                       edgelist=[(u, v) for u, v in g.edges(sub) if u in sub and v in sub],
                       edge_color='k', width=3.0, arrows=False)
nx.draw_networkx_nodes(g, pos, nodelist=deps,
                       node_color=[rgb2hex(colors[lookup[v]]) for v in deps],
                       linewidths=0, node_size=120)
nx.draw_networkx_nodes(g, pos, nodelist=[pkg],
                       node_color=rgb2hex(colors[lookup[pkg]]),
                       linewidths=0, node_size=120, node_shape='s')
xs = [x for x, y in pos.values()]
ys = [y for x, y in pos.values()]
plt.xlim((min(xs) - 10, max(xs) + 10))
plt.ylim((min(ys) - 10, max(ys) + 10))
# remove whitespace
plt.axis('off')
plt.gca().xaxis.set_major_locator(NullLocator())
plt.gca().yaxis.set_major_locator(NullLocator())

plt.savefig(snakemake.output[0], bbox_inches='tight')
