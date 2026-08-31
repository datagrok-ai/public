#!/usr/bin/env python3
"""Bind the local-install compose files to a release train, and verify the binding.

A compose file and the images it starts must come from the same release train. The compose
file carries GROK_PARAMETERS, which the datlas build inside `datagrok/datagrok` parses; the
two evolve together. Pairing master's compose with a release image (or the reverse) makes
the server exit before it starts -- e.g. master moved `connectorsSettings` from an object to
an array, and 1.27.x datlas dies with `Invalid argument(s): useGrokConnect`.

The binding is expressed in the compose files themselves:

  * `x-datagrok-release` names the train -- `bleeding-edge` on master, `X.Y.Z` on
    `release/X.Y.Z`.
  * every `datagrok/*` image that ships on that train defaults to a tag from it, through a
    per-service `*_VERSION` variable.

`check` enforces both against the branch. `stamp` rewrites them when a release branch is cut
or promoted.

    python stamp_release.py check
    python stamp_release.py check --expect 1.27.8
    python stamp_release.py stamp --version 1.27.8
    python stamp_release.py stamp --version 1.27.8 --set grok_pipe=1.22.0 --dry-run

Exits non-zero on any violation, so it can gate CI and the release pipeline.
"""

import argparse
import glob
import io
import json
import os
import re
import subprocess
import sys
import urllib.request

# Images built from this repo's release train. Each maps to the variable that overrides its
# tag; the variable's default is what binds the file to a train. Anything not listed here
# (rabbitmq, pgvector, demo_db_*, mlflow) versions independently and is pinned outright.
TRAIN_SERVICES = {
    'datagrok': 'DATAGROK_VERSION',
    'grok_pipe': 'GROK_PIPE_VERSION',
    'grok_spawner': 'GROK_SPAWNER_VERSION',
    'grok_connect': 'GROK_CONNECT_VERSION',
    'grok_connect_adbc': 'GROK_CONNECT_ADBC_VERSION',
    'jkg_python': 'JKG_PYTHON_VERSION',
    'jkg_r': 'JKG_R_VERSION',
    'jkg_octave': 'JKG_OCTAVE_VERSION',
    'jkg_julia': 'JKG_JULIA_VERSION',
    'jkg_nodejs': 'JKG_NODEJS_VERSION',
}

COMPOSE_GLOB = 'localhost*.docker-compose.yaml'
INSTALL_SCRIPTS = [
    'datagrok-install-local.sh',
    'datagrok-install-local-macos-silicon.sh',
    'datagrok-install-local.cmd',
]

BLEEDING_EDGE = 'bleeding-edge'
SEMVER_RE = re.compile(r'^\d+\.\d+\.\d+$')
# These files are CRLF, so every end-of-line class has to admit a stray \r.
MARKER_RE = re.compile(r'^x-datagrok-release:[ \t]*(\S+)[ \t\r]*$', re.M)
IMAGE_RE = re.compile(r'^([ \t]*image:[ \t]*)datagrok/([a-z_0-9]+):(\S+?)[ \t\r]*$', re.M)
# `set datagrok_default_version=X` (cmd) or `datagrok_default_version="X"` (sh)
SH_DEFAULT_RE = re.compile(r'^(datagrok_default_version=")([^"]*)(")[ \t\r]*$', re.M)
CMD_DEFAULT_RE = re.compile(r'^(set datagrok_default_version=)(\S*?)[ \t\r]*$', re.M)


def read(path):
    return io.open(path, encoding='utf-8', newline='').read()


def write(path, text):
    io.open(path, 'w', encoding='utf-8', newline='').write(text)


def compose_files(d):
    return sorted(glob.glob(os.path.join(d, COMPOSE_GLOB)))


def current_branch(d):
    """The branch being built. CI usually checks out a detached HEAD, so fall back to the
    name the CI system exports."""
    try:
        branch = subprocess.check_output(
            ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
            cwd=d, stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        branch = ''
    if branch and branch != 'HEAD':
        return branch
    for var in ('BRANCH_NAME', 'GIT_BRANCH', 'GITHUB_REF_NAME', 'CHANGE_BRANCH'):
        v = os.environ.get(var, '').strip()
        if v and v != 'HEAD':
            return re.sub(r'^(origin|refs/heads)/', '', v)
    return None


def current_train(d):
    """The train this working tree should declare, inferred from the branch."""
    branch = current_branch(d)
    if branch is None:
        return None
    m = re.match(r'^release/(\d+\.\d+\.\d+)$', branch)
    return m.group(1) if m else (BLEEDING_EDGE if branch == 'master' else None)


def declared_train(text):
    m = MARKER_RE.search(text)
    return m.group(1) if m else None


def script_default(text):
    for rx in (SH_DEFAULT_RE, CMD_DEFAULT_RE):
        m = rx.search(text)
        if m:
            return m.group(2)
    return None


def latest_semver(repo):
    """Newest X.Y.Z tag for datagrok/<repo> on Docker Hub, or None if it has none."""
    url = ('https://hub.docker.com/v2/repositories/datagrok/%s/tags'
           '?page_size=100&ordering=last_updated' % repo)
    best = None
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode())
    except Exception as e:
        raise SystemExit('error: could not query Docker Hub for %s: %s' % (repo, e))
    for t in data.get('results', []):
        name = t.get('name', '')
        if SEMVER_RE.match(name):
            key = tuple(int(p) for p in name.split('.'))
            if best is None or key > best[0]:
                best = (key, name)
    return best[1] if best else None


def tag_exists(repo, tag):
    url = 'https://hub.docker.com/v2/repositories/datagrok/%s/tags/%s' % (repo, tag)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return r.status == 200
    except Exception:
        return False


def resolve_moving_tag(repo, tag):
    """Which X.Y.Z tag a moving tag currently points at, by matching manifest digests.

    The install scripts resolve a tag by reading the image's `branch` label, which needs the
    image pulled. Here we only need to confirm a moving default is usable, so digest matching
    against the registry is enough and costs nothing.
    """
    url = ('https://hub.docker.com/v2/repositories/datagrok/%s/tags'
           '?page_size=100&ordering=last_updated' % repo)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read().decode())
    except Exception:
        return None
    digests = {}
    for t in data.get('results', []):
        if t.get('digest'):
            digests.setdefault(t['digest'], []).append(t['name'])
    for names in digests.values():
        if tag in names:
            semver = sorted(n for n in names if SEMVER_RE.match(n))
            if semver:
                return semver[-1]
    return None


def branch_exists(branch):
    try:
        out = subprocess.check_output(
            ['git', 'ls-remote', '--heads', 'https://github.com/datagrok-ai/public.git', branch],
            stderr=subprocess.DEVNULL, timeout=60).decode()
        return bool(out.strip())
    except Exception:
        return False


def check(d, expect, online=False, compose_only=False):
    train = expect or current_train(d)
    if train is None:
        print('error: cannot infer the release train from the current branch; pass --expect')
        return 1

    problems = []
    files = compose_files(d)
    if not files:
        print('error: no %s found in %s' % (COMPOSE_GLOB, d))
        return 1

    for path in files:
        name = os.path.basename(path)
        text = read(path)
        got = declared_train(text)
        if got is None:
            problems.append('%s: no x-datagrok-release marker' % name)
        elif got != train:
            problems.append('%s: declares train %r, branch is %r' % (name, got, train))

        for m in IMAGE_RE.finditer(text):
            repo, tag = m.group(2), m.group(3)
            if repo not in TRAIN_SERVICES:
                if '$' in tag:
                    problems.append('%s: %s is not a release-train image but its tag is a '
                                    'variable (%s)' % (name, repo, tag))
                continue
            var = TRAIN_SERVICES[repo]
            default = tag_default(tag, var)
            if default is None:
                problems.append('%s: %s tag is %s, expected ${%s:-<tag>}' % (name, repo, tag, var))
                continue
            if default == 'latest':
                problems.append('%s: %s defaults to the floating "latest" tag, which is not '
                                'bound to any train' % (name, repo))
            elif repo == 'datagrok':
                # The core image *is* the train, so its tag has to name it exactly.
                if default != train:
                    problems.append('%s: datagrok defaults to %r, but the train is %r'
                                    % (name, default, train))
            elif train == BLEEDING_EDGE:
                if default != BLEEDING_EDGE:
                    problems.append('%s: %s defaults to %r on the bleeding-edge train'
                                    % (name, repo, default))
            elif default == BLEEDING_EDGE:
                # Companion services version independently of the core, so their tags need
                # not equal the train -- but on a release train they must be pinned. One
                # with no published release may stay floating only if it is also kept out
                # of the "all" profile, so a stable install never pulls it by accident.
                if in_all_profile(text, repo):
                    problems.append('%s: %s is unpinned (bleeding-edge) on release train %s '
                                    'and still rides the "all" profile' % (name, repo, train))

    # On a release train the install scripts must default to that release. On master they
    # default to the newest release -- the script URL is stable, only the compose moves.
    # Release branches cut before the binding existed carry the old scripts, which have no
    # default at all; --compose-only lets the release pipeline stamp those branches without
    # the legacy scripts failing the run.
    defaults = {}
    for fn in ([] if compose_only else INSTALL_SCRIPTS):
        path = os.path.join(d, fn)
        if not os.path.exists(path):
            problems.append('%s: missing' % fn)
            continue
        got = script_default(read(path))
        if got is None:
            problems.append('%s: no datagrok_default_version' % fn)
            continue
        defaults[fn] = got
        # The default is a Docker tag, not a train: `latest` is the normal choice, since the
        # scripts resolve whatever it points at and fetch that release's compose file. An
        # exact X.Y.Z pins the scripts to one release. On a release branch the scripts ship
        # alongside that release, so they must install it and not drift to a newer one.
        if train != BLEEDING_EDGE:
            if got not in (train, 'latest'):
                problems.append('%s: installs %r by default, branch is %r' % (fn, got, train))
        elif got != 'latest' and got != BLEEDING_EDGE and not SEMVER_RE.match(got):
            problems.append('%s: default %r is not a tag the scripts can resolve '
                            '(latest, X.Y.Z, or bleeding-edge)' % (fn, got))
    if len(set(defaults.values())) > 1:
        problems.append('install scripts disagree on the default version: %s'
                        % json.dumps(defaults, sort_keys=True))

    # The default is only useful if both halves of the pair it names actually exist: the
    # image on Docker Hub and the branch the compose file is fetched from. A version that
    # has been branched for RC but not yet promoted satisfies neither.
    if online:
        for tag in sorted(set(defaults.values())):
            if tag == BLEEDING_EDGE:
                continue
            if not tag_exists('datagrok', tag):
                problems.append('install default %s: datagrok/datagrok:%s is not published'
                                % (tag, tag))
                continue
            # A moving tag is only as good as what it points at, so check that release.
            version = tag if SEMVER_RE.match(tag) else resolve_moving_tag('datagrok', tag)
            if version is None:
                problems.append('install default %s: cannot tell which release it points at'
                                % tag)
            elif not branch_exists('release/%s' % version):
                problems.append('install default %s resolves to %s, but branch release/%s does '
                                'not exist, so the compose file cannot be downloaded'
                                % (tag, version, version))
            elif tag != version:
                print('install default %s currently resolves to %s' % (tag, version))

    print('release train: %s' % train)
    print('compose files: %s' % ', '.join(os.path.basename(f) for f in files))
    if problems:
        print('\n%d problem(s):' % len(problems))
        for p in problems:
            print('  - %s' % p)
        return 1
    print('\nOK: every compose file and install script is bound to %s' % train)
    return 0


def stamp(d, version, overrides, allow_floating, dry_run, verify_tags, force):
    if not SEMVER_RE.match(version) and version != BLEEDING_EDGE:
        raise SystemExit('error: --version must be X.Y.Z or %s' % BLEEDING_EDGE)

    # Stamping only pins tags -- it cannot make a compose body compatible with an older
    # release. Pinning master's body to X.Y.Z produces exactly the pairing this whole
    # mechanism exists to prevent, so stamping is confined to the branch that owns the
    # train: cut release/X.Y.Z from master first, then stamp it there.
    branch_train = current_train(d)
    if not force and branch_train is not None and branch_train != version:
        raise SystemExit(
            'error: this checkout is on the %r train but --version is %r.\n'
            '       Stamping would pin %r images to a compose body written for %r.\n'
            '       Cut and check out release/%s first, then stamp it.\n'
            '       Use --force only if you are certain the body already matches.'
            % (branch_train, version, version, branch_train, version))

    # Resolve the tag each service gets on this train. Companion services version
    # independently of the core (grok_connect 2.6.x, grok_pipe 1.22.x), so each is resolved
    # on its own -- pinning them all to the core version names tags that don't exist.
    resolved = {}
    for repo in sorted(TRAIN_SERVICES):
        if repo in overrides:
            resolved[repo] = overrides[repo]
        elif version == BLEEDING_EDGE:
            resolved[repo] = BLEEDING_EDGE
        elif repo == 'datagrok':
            resolved[repo] = version
        else:
            tag = latest_semver(repo)
            if tag is None:
                if repo in allow_floating:
                    print('warning: datagrok/%s has no released X.Y.Z tag; leaving it on %s'
                          % (repo, BLEEDING_EDGE))
                    resolved[repo] = BLEEDING_EDGE
                else:
                    raise SystemExit(
                        'error: datagrok/%s has no released X.Y.Z tag on Docker Hub.\n'
                        '       Publish one, or re-run with --allow-floating %s to keep it '
                        'on %s.' % (repo, repo, BLEEDING_EDGE))
            else:
                resolved[repo] = tag

    if verify_tags:
        for repo, tag in sorted(resolved.items()):
            if not tag_exists(repo, tag):
                raise SystemExit('error: datagrok/%s:%s does not exist on Docker Hub' % (repo, tag))

    print('stamping train %s' % version)
    for repo in sorted(resolved):
        flag = '  (floating)' if resolved[repo] == BLEEDING_EDGE and version != BLEEDING_EDGE else ''
        print('  %-24s -> %s%s' % (repo, resolved[repo], flag))
    print('')

    changed = []
    for path in compose_files(d):
        text = read(path)
        nl = '\r\n' if '\r\n' in text else '\n'
        out = text

        def repl(m):
            indent, repo = m.group(1), m.group(2)
            if repo not in resolved:
                return m.group(0)
            return '%sdatagrok/%s:${%s:-%s}' % (indent, repo, TRAIN_SERVICES[repo], resolved[repo])

        out = IMAGE_RE.sub(repl, out)

        if MARKER_RE.search(out):
            out = MARKER_RE.sub('x-datagrok-release: %s' % version, out, count=1)
        else:
            out = ('x-datagrok-release: %s%s%s' % (version, nl, nl)) + out

        # A service kept on a floating tag must not ride the `all` profile on a release
        # train, or a stable install silently pulls an unversioned image.
        for repo, tag in resolved.items():
            if tag == BLEEDING_EDGE and version != BLEEDING_EDGE:
                out = demote_from_all_profile(out, repo)

        if out != text:
            changed.append(os.path.basename(path))
            if not dry_run:
                write(path, out)

    for fn in INSTALL_SCRIPTS:
        path = os.path.join(d, fn)
        if not os.path.exists(path):
            continue
        text = read(path)
        out = SH_DEFAULT_RE.sub(lambda m: m.group(1) + version + m.group(3), text)
        out = CMD_DEFAULT_RE.sub(lambda m: m.group(1) + version, out)
        if out != text:
            changed.append(fn)
            if not dry_run:
                write(path, out)

    if not changed:
        print('nothing to change')
    else:
        print('%s %d file(s):' % ('would update' if dry_run else 'updated', len(changed)))
        for c in changed:
            print('  %s' % c)
    return 0


def set_install_default(d, version, dry_run):
    """Point the install scripts at a release without touching the compose files.

    Run on master when a release is promoted: the scripts keep a stable download URL and
    install the newest release, while the compose file still comes from that release's
    branch. Stamping the compose files on master would be wrong -- master's body belongs
    to the bleeding-edge train.
    """
    if not SEMVER_RE.match(version):
        raise SystemExit('error: --version must be X.Y.Z')
    changed = []
    for fn in INSTALL_SCRIPTS:
        path = os.path.join(d, fn)
        if not os.path.exists(path):
            raise SystemExit('error: %s is missing' % fn)
        text = read(path)
        out = SH_DEFAULT_RE.sub(lambda m: m.group(1) + version + m.group(3), text)
        out = CMD_DEFAULT_RE.sub(lambda m: m.group(1) + version, out)
        if out == text:
            continue
        changed.append(fn)
        if not dry_run:
            write(path, out)
    print('%s install default -> %s' % ('would set' if dry_run else 'set', version))
    for c in changed:
        print('  %s' % c)
    if not changed:
        print('  (already current)')
    return 0


def tag_default(tag, var):
    """The default inside `${VAR:-default}`, or None if the tag isn't in that form."""
    m = re.match(r'^\$\{%s:-([^}]*)\}$' % re.escape(var), tag)
    return m.group(1) if m else None


def service_block_re(repo):
    # CRLF throughout, so line breaks are matched as \r?\n rather than \n.
    return re.compile(r'(^[ \t]*%s:\r?\n(?:[ \t]+.*\r?\n|\r?\n)*?[ \t]+profiles:[ \t]*\[)'
                      r'([^\]]*)(\])' % re.escape(repo), re.M)


def in_all_profile(text, repo):
    m = service_block_re(repo).search(text)
    if not m:
        return False
    return 'all' in [i.strip().strip('"\'') for i in m.group(2).split(',')]


def demote_from_all_profile(text, repo):
    """Drop "all" from the profiles list of the service block named `repo`."""
    def repl(m):
        items = [i.strip().strip('"\'') for i in m.group(2).split(',') if i.strip()]
        if 'all' not in items:
            return m.group(0)
        items = [i for i in items if i != 'all']
        return m.group(1) + ' ' + ', '.join('"%s"' % i for i in items) + ' ' + m.group(3)

    return service_block_re(repo).sub(repl, text)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dir', default=os.path.dirname(os.path.abspath(__file__)),
                   help='directory holding the compose files (default: this script\'s)')
    sub = p.add_subparsers(dest='cmd')

    c = sub.add_parser('check', help='verify the binding matches the branch')
    c.add_argument('--expect', help='train to require instead of inferring from the branch')
    c.add_argument('--online', action='store_true',
                   help='also confirm the install default has a published image and a branch')
    c.add_argument('--compose-only', action='store_true',
                   help='check only the compose files, not the install scripts')

    s = sub.add_parser('stamp', help='rewrite the binding for a release')
    s.add_argument('--version', required=True, help='X.Y.Z, or bleeding-edge')
    s.add_argument('--set', action='append', default=[], metavar='SVC=TAG',
                   help='pin one service explicitly instead of resolving from Docker Hub')
    s.add_argument('--allow-floating', action='append', default=[], metavar='SVC',
                   help='permit a service with no released tag to stay on bleeding-edge')
    s.add_argument('--verify-tags', action='store_true',
                   help='confirm every resolved tag exists on Docker Hub')
    s.add_argument('--force', action='store_true',
                   help='stamp a train the current branch does not own (see the guard)')
    s.add_argument('--dry-run', action='store_true')

    i = sub.add_parser('set-install-default',
                       help='point the install scripts at a release (run on master)')
    i.add_argument('--version', required=True, help='X.Y.Z')
    i.add_argument('--dry-run', action='store_true')

    a = p.parse_args()
    if a.cmd == 'check':
        return check(a.dir, a.expect, a.online, a.compose_only)
    if a.cmd == 'set-install-default':
        return set_install_default(a.dir, a.version, a.dry_run)
    if a.cmd == 'stamp':
        overrides = {}
        for item in a.set:
            if '=' not in item:
                raise SystemExit('error: --set expects SVC=TAG, got %r' % item)
            k, v = item.split('=', 1)
            if k not in TRAIN_SERVICES:
                raise SystemExit('error: unknown service %r; known: %s'
                                 % (k, ', '.join(sorted(TRAIN_SERVICES))))
            overrides[k] = v
        return stamp(a.dir, a.version, overrides, set(a.allow_floating), a.dry_run,
                     a.verify_tags, a.force)
    p.print_help()
    return 1


if __name__ == '__main__':
    sys.exit(main())
