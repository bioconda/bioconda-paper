import json
import pandas as pd


packages = []
ndownloads = []
ecosystem = []
for f in snakemake.input:
    with open(f) as f:
        meta = json.load(f)
        ndownloads.append(sum(f["ndownloads"] for f in meta["files"]))
        name = meta["full_name"].split("/")[1]
        packages.append(name)

        if name.startswith("bioconductor-"):
            ecosystem.append("bioconductor")
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
                ecosystem.append("python")
            elif check_for_dep("perl") or check_for_dep("perl-threaded"):
                ecosystem.append("perl")
            else:
                ecosystem.append("other")


packages = pd.DataFrame({
    "package": packages,
    "downloads": ndownloads,
    "ecosystem": ecosystem
})

packages.to_csv(snakemake.output[0], sep="\t", index=False)
