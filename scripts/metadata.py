"""Extract name, version number and summary of all recipes.

The code in this file has been copied and modified from the following
files in the bioconda-utils repository:

bioconda_utils/utils.py  (get_meta)
bioconda_utils/linting.py  (flatten_dict, EnvMatrix, load_config)
"""
import os
import glob
import jinja2
import ruamel_yaml as yaml
from ruamel_yaml.scanner import ScannerError
from itertools import product, chain
from collections import defaultdict, Iterable


ENV_MATRIX = """CONDA_PY:
  - 27
  - 34
  - 35
CONDA_BOOST: "1.60"
CONDA_R: "3.3.1"
CONDA_PERL: "5.22.0"
CONDA_NPY: "110"
CONDA_NCURSES: "5.9"
CONDA_GSL: "1.16"
CONDA_GMP: "5.1"
"""

CONFIG = """env_matrix: example_env_matrix.yml
requirements: requirements.txt
blacklists:
    - r-blacklist
docker_image: "condaforge/linux-anvil"
docker_url: 'unix://var/run/docker.sock'
upload_channel: bioconda
channels:
    - anaconda
    - conda-forge
    - bioconda
"""

def flatten_dict(dict):
    for key, values in dict.items():
        if isinstance(values, str) or not isinstance(values, Iterable):
            values = [values]
        yield [(key, value) for value in values]

class EnvMatrix:
    """
    Intended to be initialized with a YAML file and iterated over to yield all
    combinations of environments.

    YAML file has the following format::
        CONDA_PY:
          - "2.7"
          - "3.5"
        CONDA_BOOST: "1.60"
        CONDA_PERL: "5.22.0"
        CONDA_NPY: "110"
        CONDA_NCURSES: "5.9"
        CONDA_GSL: "1.16"
    """

    def __init__(self, env):
        """
        Parameters
        ----------

        env : str or dict
            If str, assume it's a path to a YAML-format filename and load it
            into a dict. If a dict is provided, use it directly.
        """
        if isinstance(env, str):
            with open(env) as f:
                self.env = yaml.load(f)
        else:
            self.env = env
        for key, val in self.env.items():
            if key != "CONDA_PY" and not isinstance(val, str):
                raise ValueError(
                    "All versions except CONDA_PY must be strings.")

    def __iter__(self):
        """
        Given the YAML::

            CONDA_PY:
              - "2.7"
              - "3.5"
            CONDA_BOOST: "1.60"
            CONDA_NPY: "110"

        We get the following sets of env vars::

          [('CONDA_BOOST', '1.60'), ('CONDA_PY', '2.7'), ('CONDA_NPY', '110')]
          [('CONDA_BOOST', '1.60'), ('CONDA_PY', '3.5'), ('CONDA_NPY', '110')]

        A copy of the entire os.environ dict is updated and yielded for each of
        these sets.
        """
        for env in product(*flatten_dict(self.env)):
            yield env


def get_meta(recipe):
    """
    Given a package name, find the current meta.yaml file, parse it, and return
    the dict.

    Parameters
    ----------
    recipe : str
        Path to recipe (directory containing the meta.yaml file)
    """
    cfg = load_config(yaml.load(CONFIG))

    env = dict(next(iter(EnvMatrix(yaml.load(ENV_MATRIX)))))

    pth = os.path.join(recipe, 'meta.yaml')
    jinja_env = jinja2.Environment()
    content = jinja_env.from_string(
        open(pth, 'r', encoding='utf-8').read()).render(env)
    # meta = yaml.round_trip_load(content, preserve_quotes=True)
    meta = yaml.load(content, preserve_quotes=True)
    return meta



def load_config(path):
    """
    Parses config file, building paths to relevant blacklists and loading any
    specified env_matrix files.

    Parameters
    ----------
    path : str
        Path to YAML config file
    """
    # validate_config(path)

    if isinstance(path, dict):
        config = path
        relpath = lambda p: p
    else:
        config = yaml.load(open(path))
        relpath = lambda p: os.path.join(os.path.dirname(path), p)


def parse_metadata(path="../bioconda-recipes/recipes/*"):
    recipes = glob.glob("../bioconda-recipes/recipes/*")
    package_info = {}
    for recipe in recipes:
        try:
            recipe_meta = get_meta(recipe)
        except FileNotFoundError:
            f = glob.glob(recipe + "/*/meta.yaml")
            if len(f) == 1:
                meta_file = f[0]
                recipe_meta = get_meta(os.path.split(meta_file)[0])
            else:
                print("Error finding meta.yaml for {}".format(recipe))
        except ScannerError:
            print("Error, could not parse: {}".format(recipe))
            
            continue
        try:
            name = recipe_meta["package"]["name"] 
            version = recipe_meta["package"]["version"]
            summary = recipe_meta["about"]["summary"]
            package_info[name] = (version, summary)
        except KeyError:
            print("Cound not find one of the keys for {}".format(recipe))
    return package_info


def main():
    for name, (version, summary) in parse_metadata().items():
        print("{:<30} {:>20}   {}".format(name, version, summary[0:80]))

if __name__ == "__main__":
    main()
