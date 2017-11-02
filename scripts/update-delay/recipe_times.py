import re
import git
from collections import namedtuple
from os import environ
from jug import TaskGenerator, bvalue

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

@TaskGenerator
def get_pypi_metadata(package, version):
    import requests
    import json
    import time
    time.sleep(.1) # Do not bombard the server
    return json.loads(requests.get('https://pypi.python.org/pypi/{}/{}/json'.format(package, version)).text)

@TaskGenerator
def get_github_releases_metadata(owner, repo):
    import requests
    import json
    import time
    time.sleep(.1) # Do not bombard the server
    if 'GITHUB_USERNAME' not in environ or \
        'GITHUB_TOKEN' not in environ:
        from sys import exit
        print("Please set the environmental variables GITHUB_USERNAME and GITHUB_TOKEN to query Github's API")
        exit(1)

    r = requests.get('https://api.github.com/repos/{}/{}/releases'.format(owner, repo),
                auth=(environ['GITHUB_USERNAME'], environ['GITHUB_TOKEN']))
    return json.loads(r.text)

@TaskGenerator
def find_github_date(gdata, url):
    import iso8601
    dates = [(e['tag_name'],e['published_at']) for e in gdata]
    tag = re.match(r'^https?://github.com/[^/]+/[^/]+/archive/(.*)\.(tar.gz|tar.bz2|tgz|zip)$', url).groups()[0]
    for t,d in dates:
        if t == tag:
            upload_date = iso8601.parse_date(d)
            break
    else:
        raise ValueError("Could not find tag: " + tag)
    is_uptodate = True
    for t,d in dates:
        d = iso8601.parse_date(d)
        if d > upload_date:
            is_uptodate = False
    return upload_date, is_uptodate

@TaskGenerator
def find_pypi_release_date(pypimeta, ver):
    import iso8601
    if ver not in pypimeta['releases']:
        raise NotImplementedError("")
    upload_time = iso8601.parse_date(pypimeta['releases'][ver][0]['upload_time'])
    is_uptodate = True
    for rs in pypimeta['releases'].values():
        for r in rs:
            if iso8601.parse_date(r['upload_time']) > upload_time:
                is_uptodate = False
    return upload_time, is_uptodate

@TaskGenerator
def find_last_mtime(pname, ver, url):
    from os import stat, path
    from os import walk
    import subprocess
    from os import makedirs
    from glob import glob

    if not (url.endswith('.tar.gz') or
            url.endswith('.tgz') or
            url.endswith('.tar.bz2') or
            url.endswith('.tar.xz') or
            url.endswith('.zip')):
        raise NotImplementedError("Cannot handle: " + url)
    wdir = 'work-space/{}-{}'.format(pname, ver)
    try:
        makedirs(wdir)
        makedirs(wdir + '/expanded')
    except:
        pass
    val = subprocess.call(['wget', url, '-t', '5', '-P', wdir])
    if val:
        return 'wget failed'
    tgzs = glob(wdir + '/*.tar.gz') + glob(wdir + '/*.tgz')
    bgzs = glob(wdir + '/*.tar.bz2')
    xgzs = glob(wdir + '/*.tar.xz')
    zips = glob(wdir + '/*.zip')
    if len(tgzs) == 1:
        [t] = tgzs
        val = subprocess.check_call(['tar', 'xzvf', t, '-C', wdir + '/expanded'])
    elif len(bgzs) == 1:
        [t] = bgzs
        val = subprocess.check_call(['tar', 'xjvf', t, '-C', wdir + '/expanded'])
    elif len(xgzs) == 1:
        [t] = xgzs
        val = subprocess.check_call(['tar', 'xJvf', t, '-C', wdir + '/expanded'])
    elif len(zips) == 1:
        [z] = zips
        val = subprocess.check_call(['unzip', '-d', wdir + '/expanded', z])
    else:
        raise IOError(wdir)
    if val:
        raise IOError()
    last_time = 0
    for basedir, _, fs in walk(wdir + '/expanded'):
        for f in fs:
            try:
                mt = stat(path.join(basedir,f)).st_mtime
                if mt > last_time:
                    last_time = mt
            except:
                if not path.islink(path.join(basedir, f)):
                    raise
    return last_time

if "BIOCONDA_RECIPES_DIRECTORY" not in environ:
    from sys import exit
    print("Please set the environmental variable BIOCONDA_RECIPES_DIRECTORY to point to the directory where bioconda-recipes is checked out")
    exit(1)

earliest = bvalue(retrieve_earliest_commit())

gdates = []
pdates = []

table = []
for pname, ver, url, source, timestamp in earliest:
    if type(url) == list:
        url = url[0]
    if url.startswith('http://hgdownload.cse.ucsc.edu/admin/exe/macOSX'):
        continue
    pypi_match1  = re.match(r'https?://pypi.python.org/packages/source/.?/([^/]+)/([^/]+)', url)
    pypi_match2  = re.match(r'https?://pypi.io/packages/source/.?/([^/]+)/([^/]+)', url)
    github_match = re.match(r'^https?://github.com/([-._a-zA-Z0-9]+)/([-._a-zA-Z0-9]+)/archive(/release)?/(.*)$', url)
    if pypi_match1 or pypi_match2:
        match = pypi_match1 if pypi_match1 else pypi_match2
        package, fname = match.groups()
        pypimeta = get_pypi_metadata(package, ver)
        pdate = find_pypi_release_date(pypimeta, ver)
        table.append((pname, ver, timestamp, pdate))
    elif github_match:
        owner, repo, _, uver = github_match.groups()
        if uver.endswith('.tar.gz'):
            uver = uver[:-len('.tar.gz')]
        elif uver.endswith('.tar.bz2'):
            uver = uver[:-len('.tar.bz2')]
        meta = get_github_releases_metadata(owner, repo)
        gdate = find_github_date(meta, url)
        table.append((pname, ver, timestamp, gdate))
    elif type(url) == str:
        mtime = find_last_mtime(pname, ver, url)
        table.append((pname, ver, timestamp, (mtime, None)))
    else:
        raise ValueError("Cannot handle this type of URL: "+ url)

