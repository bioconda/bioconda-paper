from snakemake.shell import shell
import pandas as pd
import seaborn as sns
import networkx as nx
import glob
import os
from networkx.drawing.nx_agraph import read_dot, write_dot
from matplotlib.colors import rgb2hex
from pydotplus import parser

shell.prefix('')

packages = pd.read_table(snakemake.input.pkg, index_col=0)
packages.loc[packages.ecosystem == 'Bioconductor', 'ecosystem'] = 'Bioconductor/R'
packages.loc[packages.ecosystem == 'R', 'ecosystem'] = 'Bioconductor/R'
lookup = packages['ecosystem'].to_dict()
colors = dict(zip(['Bioconductor/R', 'Other', 'Python', 'Perl'], sns.color_palette('colorblind')))
g = read_dot(snakemake.input.dag)
dot = parser.parse_dot_data(open(snakemake.input.dag).read())
pkg = snakemake.wildcards.pkg
sub = set(nx.ancestors(g, pkg))
sub.add(pkg)
for v in sub:

    color = rgb2hex(colors[lookup[v]]) + "ff" # be opaque
    v_ = g.nodes[v]
    v_['fillcolor'] = color
    v_['style'] = 'filled'
    v_['width'] = 0.4
    if v == pkg:
        v_['shape'] = 'square'

for u, v in g.edges(sub):
    if u in sub and v in sub:
        e = g.edges[u, v]
        e['color'] = '#000000'
        e['penwidth'] = 5.0

write_dot(g, snakemake.output[0])
