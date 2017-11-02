import git
from collections import namedtuple
from os import environ
from jug import TaskGenerator

RecipeInfo = namedtuple('RecipeInfo', ['pname', 'version', 'url', 'source'])


def version_greater(va, vb):
    '''Compare two versions'''
    def norm_version(v):
        return tuple([vb for vb in str(v).split('.')])

    nva = norm_version(va)
    nvb = norm_version(vb)
    for ea, eb in zip(nva, nvb):
        if ea == eb:
            continue
        try:
            return int(ea) > int(eb)
        except:
            return ea > eb
    return False

def extract_meta(meta):
    '''Tries to parse the contents of a meta.yaml file to build a RecipeInfo
    object'''
    import yaml
    import jinja2
    try:
        data = yaml.load(jinja2.Template(meta.decode('utf-8')).render())
    except:
        print("Could not parse")
        return None
    if data is None:
        return None
    if 'source' not in data or not data['source']:
        return None
    if 'version' not in data['package']:
        return None
    pname = data['package']['name']
    version = data['package']['version']
    source = data['source']
    if 'url' in source:
        url = source['url']
    elif 'git_url' in source:
        url = source['git_url'] + '#'
        if 'git_rev' in source:
            url += str(source['git_rev'])
        elif 'git_tag' in source:
            url += str(source['git_tag'])
        else:
            url += 'HEAD'
            print(source)
    elif 'hg_url' in source:
        url = source['hg_url']
        if 'hg_tag' in source:
            url += '#' + str(source['hg_tag'])
        else:
            print(data)
    else:
        print(data)
        return None
    return RecipeInfo(pname, version, url, source)

def extract_recipes(c, cache):
    '''Extracts all RecipeInfo objects'''
    from six import BytesIO
    try:
        trees = c.tree['recipes'].trees
    except:
        try:
            trees = c.trees
        except:
            return []
    res = []
    while trees:
        st = trees.pop()
        trees.extend(st.trees)
        for b in st.blobs:
            if b.name == 'meta.yaml':
                # both the hex hash and the actual data are kept in the cache.
                if b.hexsha not in cache:
                    bout = BytesIO()
                    b.stream_data(bout)
                    data = bout.getvalue()
                    if data not in cache:
                        m = extract_meta(data)
                        cache[data] = m
                    else:
                        m = cache[data]
                    cache[b.hexsha] = m
                reptime = cache[b.hexsha]
                if reptime is not None:
                    res.append(reptime)
    return res


def walk_repo(repo, init_commit=None):
    '''Iterate over all commits upstream of `init_commit``'''
    if init_commit is None:
        queue = [repo.head.commit]
    else:
        queue = [repo.commit(init_commit)]
    seen = set()
    while queue:
        c = queue.pop()
        if c not in seen:
            seen.add(c)
            queue.extend(c.parents)
            yield c

@TaskGenerator
def retrieve_earliest_commit():
    repo = git.Repo(environ["BIOCONDA_RECIPES_DIRECTORY"])
    cache = {}
    earliest = {}
    initial = True
    for c in walk_repo(repo, '6ef812f1b808ac7dc39313ff133ed42f027bebdb'):
        timestamp = c.committed_date
        for repinfo in extract_recipes(c, cache):
            if repinfo.pname in earliest:
                cur = earliest[repinfo.pname]
                if cur[0].version == repinfo.version:
                    if cur[1] > timestamp:
                        cur[1] = timestamp
                elif version_greater(repinfo.version, cur[0].version):
                    if initial:
                        earliest[repinfo.pname] = [repinfo, timestamp]
            else:
                earliest[repinfo.pname] = [repinfo, timestamp]
        initial = False
    data = [(p.pname, p.version, p.url, p.source, t) for p, t in earliest.values()]
    return data

if "BIOCONDA_RECIPES_DIRECTORY" not in environ:
    from sys import exit
    print("Please set the environmental variable BIOCONDA_RECIPES_DIRECTORY to point to the directory where bioconda-recipes is checked out")
    exit(1)

earliest = retrieve_earliest_commit()
