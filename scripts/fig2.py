from svgutils.compose import *
from common import label

Figure(
    "23.1cm", "6cm",
    Panel(SVG(snakemake.input.dag).scale(0.24), label("a")),
    Panel(SVG(snakemake.input.workflow).scale(0.5), label("b")).move(220, 0),
    Panel(SVG(snakemake.input.turnaround).scale(0.9), label("c")).move(220, 50),
    Panel(SVG(snakemake.input.comp).scale(0.8), label("d")).move(570, 0),
    #Grid(40, 40)
).save(snakemake.output[0])
