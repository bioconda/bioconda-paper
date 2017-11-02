# Update Delay

Parse bioconda packages to discover the last update time of the package and the
upstream release.

## Environmental variables

The following variables must be set for the script to run:

- BIOCONDA_RECIPES_DIRECTORY: a git checkout of the bioconda repository
- GITHUB_USERNAME/GITHUB_TOKEN: a github username/authentication token to allow for API queries without hitting the rate limit

## Notes

- several packages get a master branch. Thus, it is not clear what "last release" would mean.

