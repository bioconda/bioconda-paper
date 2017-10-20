from svgutils.compose import *
from common import label

Figure(
    "23.7cm", "6.1cm",
    Panel(SVG(snakemake.input.dag).scale(0.6), label("a")),
    Panel(SVG(snakemake.input.workflow).scale(0.5), label("b")).move(220, 0),
    Panel(SVG(snakemake.input.turnaround).scale(0.9).move(5, 0), label("c")).move(220, 50),
    Panel(SVG(snakemake.input.comp).scale(0.9).move(10, 0), label("d")).move(570, 0),
    Panel(SVG(snakemake.input.age).scale(0.9).move(19, 0), label("e")).move(570, 90),
    #Grid(40, 40)
).save(snakemake.output[0])
