import numpy as np
import datetime
from jug import value, set_jugdir
set_jugdir('recipe_times.jugdata')
from recipe_times import table
deltas = []
outdated = []
tracking = []
for r in table:
    try:
        package,_,timestamp, upstream = value(r)
    except:
        continue
    if upstream is None:
        continue
    if type(upstream[0]) in [float, int]:
        upstream_date = datetime.datetime.fromtimestamp(upstream[0])
    elif type(upstream[0]) == str:
        continue
    else:
        upstream_date = upstream[0]
    if upstream_date.year < 2016:
        continue
    delta = datetime.datetime.fromtimestamp(timestamp).toordinal() - upstream_date.toordinal()
    if delta < 0:
        continue
    deltas.append(delta)
    if upstream[1] is not None:
        (tracking if upstream[1] else outdated).append(package)
deltas = np.array(deltas)
with open('results.txt', 'wt') as output:
    output.write('''\
Mean (+/- std.dev.) number of days between upstream release and bioconda package: {mean} +/- {std_dev}
Median number of days between upstream release and bioconda package: {median}
Based on {N} packages.

On those packages where it could be heuristically determined, {tracking} of {can_track} are current with their upstream release.
'''.format(mean=deltas.mean(),
            std_dev=deltas.std(),
            median=np.median(deltas),
            N=len(deltas),
            tracking=len(tracking),
            can_track=len(tracking)+len(outdated)))
