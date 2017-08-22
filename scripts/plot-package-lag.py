from pkg_resources import parse_version # requires setuptools
from distutils.version import LooseVersion 
import pickle
import json
import difflib
import distance
import glob
import os

from metadata import parse_metadata


def load_bioconda_packages(force=False):
    """Load bioconda package information from the meta.yaml files.

    Pickle the result to improve access for subsequent runs.

    Arguments:
        force (bool): Overwrite an existing pickled file instead of loading it.
    """
    if os.path.exists("bioconda_packages.pickle") and not force:
        bioconda_packages = pickle.load(open("bioconda_packages.pickle", "rb"))
        print("Found and loaded pickled bioconda packages")
        return bioconda_packages
    else:
        bioconda_packages = parse_metadata()
        pickle.dump(bioconda_packages, open("bioconda_packages.pickle", "wb"))
        return bioconda_packages


def load_ubuntu_packages(force=True):
    """Parse ubuntu package information from a file genrated
    with `apt-cache show | grep '^Package\|^Description-en\|^Version'` on ubuntu 16.04,

    Pickle the result to improve access for subsequent runs.

    Arguments:
        force (bool): Overwrite an existing pickled file instead of loading it.
    """
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
            # parse blocks of three subsequent lines
            pstring = package_strings[index]
            if pstring.startswith("Package"):
                name = pstring.split(" ")[1]
            elif pstring.startswith("Version"):
                version = pstring.split(" ")[1]
            elif pstring.startswith("Description"):
                description = pstring.split(" ", 1)[1]
                ubuntu_packages[name] = (version, description)
        # print("Done parsing")
        pickle.dump(ubuntu_packages, open("ubuntu_packages.pickle", "wb"))
        return ubuntu_packages 


def normalize(vnr):
    """Normalize version numbers for comparison.
    """
    if "-" in vnr:
        vnr = vnr.split("-")[0]
    if "~" in vnr:
        vnr = vnr.split("~")[0]
    if "+" in vnr:
        vnr = vnr.split("+")[0]
    if ":" in vnr:
        _, vnr  = vnr.split(":")
    return vnr


def split_lc(text):
    """Split and clean a description and converti it to lower case 
    for Jaccard similarity computation.
    """
    # define filling words that need to be removed
    exclude = ("the", "if", "with", "and", "that", "can", "a", "is",
               "to", "it", "for", "an", "of", "in", "or")
    # remove punctuation etc and split up composita
    words = [w.strip().lower() for w in text.split()]
    result = []
    for word in words:
        if "-" in word:
            result.extend(word.split("-"))
        else:
            result.append(word)
    # filter out filler words
    return [w for w in result if w not in exclude]
    

def compare_package_versions (bioconda_packages, ubuntu_packages):
    
    newer_in_bioconda = 0
    equal_in_both = 0
    newer_in_ubuntu = 0
    errors = 0
    not_in_ubuntu = 0
    contained = 0
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
            # lv_bc = parse_version(bc_version).base_version
            # ub_version = ubuntu_packages[bc_name]
            # lv_ub = parse_version(ub_version).base_version
            try:
                if isinstance(bc_version, str):
                    lv_bc = LooseVersion(bc_version)
                else:
                    lv_bc = LooseVersion(str(bc_version))
            except:
                print("This failed for", bc_version)
                print(f"Error loose versioning {bc_version} {type(bc_version)}")
                raise
            ub_version, ub_description = ubuntu_packages[bc_name]
            try:
                lv_ub = LooseVersion(normalize(ub_version))
            except:
                print("This failed for", ub_version)
                print(f"Error loose versioning {ub_version}")
                raise
            
            try:
                if lv_bc > lv_ub:
                    ...
            except TypeError:
                print(f"For {bc_name}: Compared {lv_bc}, ({type(lv_bc)}) and {lv_ub}, ({type(lv_ub)})")
                errors += 1
                continue

            if lv_bc > lv_ub:
                newer_in_bioconda += 1
                d = distance.jaccard(split_lc(bc_summary), split_lc(ub_description))
                # if r < 0.3:
                if d > 0.9:
                    print("     DISSIMMILAR")
                    print(f"    xxx> BC: {bc_summary}")
                    print(f"    xxx> UB: {ub_description}")
                    print(f"    xxx> I think they are not the equal_in_both package.")

            elif lv_bc == lv_ub:
                equal_in_both += 1
            else:
                print(f"  -> Ubuntu wins: {bc_name}, {lv_bc} ({bc_version}) < {lv_ub} ({ub_version})")
                # print(bc_summary, split_lc(bc_summary))
                # input()
                d = distance.jaccard(split_lc(bc_summary), split_lc(ub_description))
                # if r < 0.3:
                if d > 0.9:
                    print(f"    xxx> BC: {bc_summary}")
                    print(f"    xxx> UB: {ub_description}")
                    print(f"    xxx> I think they are not the equal_in_both package.")
                    not_in_ubuntu += 1
                else:
                    print(f"      -> BC: {bc_summary}")
                    print(f"      -> UB: {ub_description}")
                    print(f"      -> I think they look alike.")
                    newer_in_ubuntu += 1
            # except (AttributeError, TypeError):
            #     print("Error:", bc_version, ubuntu_packages[bc_name])
            #     raise
            #     errors += 1
        else:
            # find names that are contained within each other
            # example: httpretty (bioconda) is available as 
            # python-httpretty and python3-httpretty
            unmatched.append(bc_name)
            candidates = []
            for ub_name in ubuntu_packages:
                if bc_name in ub_name:
                    candidates.append(ub_name)
            if candidates:
                if len(candidates) <= 5:
                    print(f"Found containment candidates for {bc_name} {bioconda_packages[bc_name][0]}: {candidates}")
                prefixes = ["python-", "python3-", "r-cran-", "r-cran-r", "ruby-"]
                # print(set(candidates))
                # print(set([prefix + bc_name for prefix in prefixes]))
                shared = set([prefix + bc_name for prefix in prefixes]) & set(candidates)
                # print(shared)
                if shared:
                    # print(f"Best candidates are: {shared}")
                    print(f"Best candidates are: {[(s, ubuntu_packages[s][0]) for s in shared]}")
                    contained += 1
            else:
                not_in_ubuntu += 1
    
    print("\n\nStats:")
    print(f"newer in bioconda: {newer_in_bioconda}")
    print(f"newer in ubuntu: {newer_in_ubuntu}")
    print(f"same version in both: {equal_in_both}")
    print(f"errors: {errors}")
    print(f"contained: {contained}")
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
