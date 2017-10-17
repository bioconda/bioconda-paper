import seaborn as sns
from svgutils.compose import *


def label(text):
    return Text(text, 5, 15, size=12, weight="bold", font="sans-serif")


sns.set(style='ticks', context='paper', palette="tab10")

sns.set_palette('colorblind')
