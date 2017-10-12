from svgutils.compose import *
from common import label

Figure(
    "17cm", "6cm",
    Panel(SVG(snakemake.input.dag).scale(0.24), label("a")),
    Panel(SVG(snakemake.input.workflow).scale(0.45), label("b")).move(220, 0),
    Panel(SVG(snakemake.input.comp), label("c")).move(310, 0),
    #Grid(40, 40)
).save(snakemake.output[0])
