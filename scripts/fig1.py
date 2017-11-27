from svgutils.compose import *
from common import label

Figure(
    "22cm", "6cm",
    Panel(SVG(snakemake.input.ecosystems), label("a")),
    Panel(SVG(snakemake.input.downloads), label("b")).move(285, 0),
    Panel(SVG(snakemake.input.comp).scale(0.9).move(10, 0), label("c")).move(560, 0),
    Panel(SVG(snakemake.input.age).scale(0.9).move(19, 0)).move(560, 90),
    # Grid(40, 40)
).save(snakemake.output[0])
