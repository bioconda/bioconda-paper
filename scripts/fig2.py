from svgutils.compose import *
from common import label

Figure(
    "9cm", "6cm",
    Panel(SVG(snakemake.input.dag).scale(0.24), label("a")),
    Panel(SVG(snakemake.input.workflow).scale(0.45), label("b")).move(220, 0),
    # Grid(40, 40)
).save(snakemake.output[0])
