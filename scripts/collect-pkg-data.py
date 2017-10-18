import json
import pandas as pd


packages = []
ndownloads = []
ecosystem = []
versions = []
deps = []
for path in snakemake.input:
    with open(path) as f:
        meta = json.load(f)
        ndownloads.append(sum(f["ndownloads"] for f in meta["files"]))
        name = meta["full_name"].split("/")[1]
        assert name not in packages, "duplicate package: {}".format(name)
        packages.append(name)
        versions.append(len(meta["versions"]))

        if name.startswith("bioconductor-"):
            ecosystem.append("Bioconductor")
        elif name.startswith("r-"):
            ecosystem.append("R")
        else:
            def check_for_dep(dep):
                for f in meta["files"]:
                    for d in f["dependencies"]["depends"]:
                        if d["name"] == dep:
                            return True
                return False
            if check_for_dep("python"):
                ecosystem.append("Python")
            elif check_for_dep("perl") or check_for_dep("perl-threaded"):
                ecosystem.append("Perl")
            else:
                ecosystem.append("Other")

        # stores number of dependencies based on the first (0) recipe
        deps.append(len(meta['files'][0]['attrs']['depends']))


packages = pd.DataFrame({
    "package": packages,
    "downloads": ndownloads,
    "ecosystem": ecosystem,
    "versions": versions,
    "deps": deps
}, columns=["package", "ecosystem", "downloads", "versions", "deps"])

packages.sort_values("downloads", ascending=False, inplace=True)

packages.to_csv(snakemake.output[0], sep="\t", index=False)
