from svgutils.compose import *
from common import label

Figure(
    "24cm", "6.1cm",
    Panel(SVG(snakemake.input.contributions), label("a")),
    Panel(SVG(snakemake.input.add_del), label("b").move(0, -10)).move(0, 90),
    Panel(SVG(snakemake.input.dag).scale(0.6), label("c")).move(285, 0),
    Panel(SVG(snakemake.input.workflow).scale(0.5), label("d")).move(505, 0),
    Panel(SVG(snakemake.input.turnaround).scale(0.9).move(5, 0), label("e")).move(505, 50),
    Panel(SVG(snakemake.input.usage).scale(0.5), label("f")).move(505, 130)
    #Grid(40, 40)
).save(snakemake.output[0])
