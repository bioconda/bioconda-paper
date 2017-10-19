from snakemake.shell import shell
import pandas as pd
import seaborn as sns
import networkx as nx
import glob
import os
from networkx.drawing.nx_agraph import read_dot
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
sub = list(nx.ancestors(g, pkg))
for i in sub:
    node = dot.get_node(i)
    if not node:
        node = dot.get_node('"' + i + '"')
        if not node:
            continue
    node = node[0]
    node.set_fillcolor(rgb2hex(colors[lookup[i]]))
    node.set_width(0.4)
dot.write('dag/{}_colored.dot'.format(pkg))
