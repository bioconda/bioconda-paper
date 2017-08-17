from pkg_resources import parse_version # requires setuptools
from distutils.version import LooseVersion 
import pickle
import json
import difflib
import glob
import os

from metadata import parse_metadata


def load_bioconda_packages(force=False):
    if os.path.exists("bioconda_packages.pickle") and not force:
        bioconda_packages = pickle.load(open("bioconda_packages.pickle", "rb"))
        print("Found and loaded pickled bioconda packages")
        return bioconda_packages
    else:
        bioconda_packages = parse_metadata()
        pickle.dump(bioconda_packages, open("bioconda_packages.pickle", "wb"))
        return bioconda_packages


def load_ubuntu_packages(force=True):
    if os.path.exists("ubuntu_packages.pickle") and not force:
        ubuntu_packages = pickle.load(open("ubuntu_packages.pickle", "rb"))
        print("Found and loaded pickled ubuntu packages")
        return ubuntu_packages 
    else:
        ubuntu_packages = {}
        with open("../package-data/grepped_ubuntu_packages1604_description.txt", "r") as ubuntu_file:
            ubuntu_data = ubuntu_file.read()
        
        package_strings = ubuntu_data.split("\n")
        
        for index in range(len(package_strings)-1):
            pstring = package_strings[index]
            # print(pstring)
            if pstring.startswith("Package"):
                name = pstring.split(" ")[1]
            elif pstring.startswith("Version"):
                version = pstring.split(" ")[1]
            elif pstring.startswith("Description"):
                description = pstring.split(" ", 1)[1]
                ubuntu_packages[name] = (version, description)
        print("Done parsing")
        print(ubuntu_packages)
        pickle.dump(ubuntu_packages, open("ubuntu_packages.pickle", "wb"))
        return ubuntu_packages 


def normalize(vnr):
    if "-" in vnr:
        vnr = vnr.split("-")[0]
    if "~" in vnr:
        vnr = vnr.split("~")[0]
    if "+" in vnr:
        vnr = vnr.split("+")[0]
    if ":" in vnr:
        _, vnr  = vnr.split(":")
    return vnr


def compare_package_versions (bioconda_packages, ubuntu_packages):
    bigger = 0
    same = 0
    smaller = 0
    errors = 0
    not_in_ubuntu = 0
    unmatched = []
    perl_packages = []
    r_packages = []
    bioconductor_packages = []
    for bc_name, (bc_version, bc_summary) in bioconda_packages.items():

        if bc_name.startswith("perl-"):
            perl_packages.append(bc_name)
            continue
        elif bc_name.startswith("r-"):
            r_packages.append(bc_name)
            continue
        elif bc_name.startswith("bioconductor-"):
            bioconductor_packages.append(bc_name)
            continue
        elif bc_name in ubuntu_packages:
            # print(bc_name, "BC", bc_version, "Ubuntu", ubuntu_packages[bc_name], )
            try:
                # lv_bc = parse_version(bc_version).base_version
                # ub_version = ubuntu_packages[bc_name]
                # lv_ub = parse_version(ub_version).base_version
                lv_bc = LooseVersion(bc_version) 
                ub_version = ubuntu_packages[bc_name]
                lv_ub = LooseVersion(normalize(ub_version))
                if lv_bc > lv_ub:
                    bigger += 1
                elif lv_bc == lv_ub:
                    same += 1
                else:
                    smaller += 1
                    print(f"  -> Ubuntu wins: {bc_name}, {lv_bc} ({bc_version}) < {lv_ub} ({ub_version})")
            except (AttributeError, TypeError):
                print("Error:", bc_version, ubuntu_packages[bc_name])
                errors += 1
        else:
            not_in_ubuntu += 1
            unmatched.append(bc_name)
            # candidates = []
            # for ub_name in ubuntu_packages:
            #     r = difflib.SequenceMatcher(None, bc_name, ub_name).ratio()
            #     if r > 0.7:
            #         candidates.append((ub_name, r))
            # if candidates:
            #     print(f"  {bc_name} Not found in Ubuntu packages, but {candidates} are similar")
                
    
    print("\n\nStats:")
    print(f"newer in bioconda: {bigger}")
    print(f"newer in ubuntu: {smaller}")
    print(f"same: {same}")
    print(f"errors: {errors}")
    print(f"not in ubuntu: {not_in_ubuntu}")
    print(f"r-packages: {len(r_packages)}")
    print(f"perl-packages: {len(perl_packages)}")
    print(f"bioconductor-packages: {len(bioconductor_packages)}")

    # for p in sorted(unmatched):
    #     print(p)

def main():
    bioconda_packages = load_bioconda_packages() 
    ubuntu_packages = load_ubuntu_packages()
    compare_package_versions(bioconda_packages, ubuntu_packages)

if __name__ == '__main__':
    main()
    # TODO: find a way to identify mismatching packages
    # for example bioconda verse https://bioconda.github.io/recipes/verse/README.html
    # and ubuntu verse https://packages.ubuntu.com/trusty/verse
    # have nothing in common
