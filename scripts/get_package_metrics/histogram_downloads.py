#!/usr/bin/env python

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from collections import OrderedDict

download_counts = OrderedDict()

with open("download_counts.csv","r") as inf:
    for line in inf:
        name,total_downloads,days_available,downloads_per_day,days_since_newest,has_linux_build,has_osx_build = line.split(",")
        download_counts[name] = downloads_per_day
        # remove the header line
    download_counts.pop("name")

counts = [float(v) for v in download_counts.values()]

mu, sigma = np.mean(counts), np.std(counts)
x = mu + sigma*np.array(counts)

with plt.style.context("ggplot"):
    plt.yscale('symlog')
    #plt.xscale('log')

    # the histogram of the data
    n, bins, patches = plt.hist(np.array(counts), 100, facecolor='green', alpha=0.75)

    print(list(n))
    print(list(bins))

    plt.xlabel('Downloads/day', fontsize=20)
    plt.ylabel('Number of packages', fontsize=20)

    plt.title(r'$\mathrm{Bioconda\ downloads/day:}\ \mu='+'{0:.1f}'.format(mu)+',\ \sigma='+'{0:.1f}'.format(sigma)+'$', fontsize=24)
    plt.axis([0, max(counts), 0, max(n)*1.05])
    plt.grid(True)


    #plt.show()
    plt.savefig("bioconda_downloads.pdf")