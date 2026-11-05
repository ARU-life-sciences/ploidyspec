# Bioconda recipe

`ploidyspec/meta.yaml` is the recipe to submit to
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) as
`recipes/ploidyspec/meta.yaml`. For each release:

1. tag and push: `git tag -a vX.Y.Z -m ... && git push origin vX.Y.Z`;
2. set `version` and the tarball's sha256:
   `curl -sL https://github.com/ARU-life-sciences/ploidyspec/archive/refs/tags/vX.Y.Z.tar.gz | sha256sum`;
3. open the pull request on bioconda-recipes. After the first release, the
   bioconda autobump bot opens version updates by itself.

`.gitattributes` keeps `results/`, `docs/` and the other panel data out of the
tarball (a few hundred kB instead of 900 MB).
