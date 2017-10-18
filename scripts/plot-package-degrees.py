"""
Fetches all dependencies (first order) of each package and 
plots the histogram of the number of dependencies (degree)


"""
import glob
import pylab
import json



degrees = []
for filename in glob.glob("package-data/*json"):
    with open(filename, "r") as fh:
        data = json.loads(fh.read())
        degrees.append(len(data['files'][0]['attrs']['depends']))

pylab.hist(degrees, range(0,30), lw=1)
pylab.xlim([0,30])
pylab.grid()
pylab.xlabel("Package degree", fontsize=16)

pylab.savefig("plots/package_degrees.svg")
#pylab.savefig("plots/package_degrees.pdf")

