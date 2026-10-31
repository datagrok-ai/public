#!/usr/bin/env python
"""Turns `micromamba list --json` on stdin into upper-bound pip constraints on stdout.

Conda versions are not always PEP 440 (tzdata `2026e`), and uv rejects the whole
constraints file on the first invalid one, so those packages are left unconstrained.
"""
import json
import sys

from packaging.version import InvalidVersion, Version


def constraints(packages):
    lines = []
    for p in packages:
        if p['name'].startswith('_'):
            continue
        try:
            version = Version(p['version'])
        except InvalidVersion:
            continue
        lines.append(p['name'] + '<=' + version.public)
    return lines


if __name__ == '__main__':
    print('\n'.join(constraints(json.load(sys.stdin))))
