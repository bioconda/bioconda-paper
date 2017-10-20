from svgutils.compose import *
from common import label

Figure(
    "23.5cm", "6cm",
    Panel(SVG(snakemake.input.contributions), label("a")),
    Panel(SVG(snakemake.input.add_del), label("b").move(0, -10)).move(0, 90),
    Panel(SVG(snakemake.input.ecosystems), label("c")).move(285, 0),
    Panel(SVG(snakemake.input.downloads), label("d")).move(560, 0),
    # Grid(40, 40)
).save(snakemake.output[0])
