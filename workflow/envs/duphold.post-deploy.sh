#!/usr/bin/env bash
set -euo pipefail

readonly VERSION="0.2.3"
readonly SHA256="7c2f207271182953e5e5a139a5e4ad8190835fcea48fada99243ca23be92fabb"
readonly URL="https://github.com/brentp/duphold/releases/download/v${VERSION}/duphold"
readonly DEST="${CONDA_PREFIX:?CONDA_PREFIX is not set}/bin/duphold"

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

curl --fail --location --retry 3 --output "$tmp" "$URL"
echo "${SHA256}  ${tmp}" | sha256sum --check --status
install -m 0755 "$tmp" "$DEST"
