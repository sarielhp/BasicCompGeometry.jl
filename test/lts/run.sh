#!/usr/bin/env bash
# Run the test suite on the oldest supported Julia (the 1.10 LTS, juliaup channel
# "1.10") from a mirror of the working tree under output/tmp/, so the 1.10
# resolution never rewrites the main Manifest.toml. Part of tools/gate --full:
# newer compilers union-split and elide allocations that 1.10 does not.
set -euo pipefail
cd "$(dirname "$0")/../.."
if ! julia +1.10 --version >/dev/null 2>&1; then
  echo "test/lts/run.sh: Julia 1.10 is required (install with: juliaup add 1.10)" >&2
  exit 1
fi
mirror=output/tmp/lts/repo
mkdir -p "$mirror"
rsync -a --delete --exclude /Manifest.toml --exclude /output --exclude /.git \
  --exclude /docs/build --exclude /test/jet/Manifest.toml ./ "$mirror/"
cd "$mirror"
exec julia +1.10 --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
