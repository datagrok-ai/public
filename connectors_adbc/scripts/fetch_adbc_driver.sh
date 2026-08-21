#!/usr/bin/env bash
set -euo pipefail

# Fetch pre-built Apache Arrow ADBC driver libraries into lib/.
#
#   ./scripts/fetch_adbc_driver.sh snowflake bigquery
#
# Both are Go builds: neither is a crates.io crate we could link, and only Snowflake has a
# Debian package, so PyPI's wheels are the one distribution channel common to both. The wheel
# names its payload lib<driver>.so on every platform; the copy step gives it the extension the
# local loader expects. Override the destination with ADBC_LIB_DIR.

VERSION=1.11.0

if [ $# -lt 1 ]; then
  echo "usage: $0 <driver>...   (e.g. snowflake bigquery)" >&2
  exit 1
fi

case "$(uname -s)-$(uname -m)" in
  Linux-x86_64)
    PLATFORM="manylinux1_x86_64.manylinux_2_28_x86_64.manylinux_2_5_x86_64"
    LIBEXT=.so
    ;;
  Linux-aarch64)
    PLATFORM="manylinux2014_aarch64.manylinux_2_17_aarch64.manylinux_2_28_aarch64"
    LIBEXT=.so
    ;;
  Darwin-x86_64)
    PLATFORM="macosx_10_15_x86_64"
    LIBEXT=.dylib
    ;;
  Darwin-arm64)
    PLATFORM="macosx_11_0_arm64"
    LIBEXT=.dylib
    ;;
  *)
    echo "unsupported host $(uname -sm)" >&2
    exit 1
    ;;
esac

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LIB_DIR="${ADBC_LIB_DIR:-${SCRIPT_DIR}/../lib}"
mkdir -p "${LIB_DIR}"

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

for name in "$@"; do
  driver="adbc_driver_${name}"
  target="${LIB_DIR}/lib${driver}${LIBEXT}"
  if [ -f "${target}" ]; then
    echo "already fetched: ${target}"
    continue
  fi

  url="https://files.pythonhosted.org/packages/py3/a/${driver}/${driver}-${VERSION}-py3-none-${PLATFORM}.whl"
  echo "downloading ${url}"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "${url}" -o "${tmp}/wheel.zip"
  elif command -v wget >/dev/null 2>&1; then
    wget -q "${url}" -O "${tmp}/wheel.zip"
  else
    echo "need curl or wget" >&2
    exit 1
  fi

  unzip -qo "${tmp}/wheel.zip" "${driver}/lib${driver}.so" -d "${tmp}"
  cp "${tmp}/${driver}/lib${driver}.so" "${target}"
  echo "installed: ${target}"
done
