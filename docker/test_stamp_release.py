#!/usr/bin/env python3
"""Tests for stamp_release.py -- the compose/release binding.

Every case builds a throwaway git checkout from the real files in this directory, mutates
one thing, and asserts the checker's verdict. No network and no Docker images, so this runs
anywhere:

    python3 test_stamp_release.py

The cases that need release images (starting a stamped compose file against the datlas build
it was stamped for) are not here -- see the release pipeline's post-stamp check.
"""

import io
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = 'stamp_release.py'
COPY = ['stamp_release.py', 'datagrok-install-local.sh',
        'datagrok-install-local-macos-silicon.sh', 'datagrok-install-local.cmd']

failures = []


def make_tree(branch):
    d = tempfile.mkdtemp(prefix='binding-test-')
    for fn in COPY:
        shutil.copy(os.path.join(HERE, fn), d)
    for fn in os.listdir(HERE):
        if fn.startswith('localhost') and fn.endswith('.docker-compose.yaml'):
            shutil.copy(os.path.join(HERE, fn), d)
    q = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=d)
    subprocess.check_call(['git', 'init', '-q', '.'], **q)
    subprocess.check_call(['git', 'checkout', '-q', '-b', branch], **q)
    subprocess.check_call(['git', 'add', '-A'], **q)
    subprocess.check_call(['git', '-c', 'user.email=t@t', '-c', 'user.name=t',
                           'commit', '-qm', 'base'], **q)
    return d


def edit(d, fn, old, new):
    p = os.path.join(d, fn)
    t = io.open(p, encoding='utf-8', newline='').read()
    if old not in t:
        raise AssertionError('%s: pattern not present, test is stale: %r' % (fn, old))
    io.open(p, 'w', encoding='utf-8', newline='').write(t.replace(old, new))


def run(d, *args):
    r = subprocess.run([sys.executable, TOOL] + list(args), cwd=d,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return r.returncode, r.stdout.decode('utf-8', 'replace')


def case(name, branch, mutate, want_rc, args=('check',)):
    d = make_tree(branch)
    try:
        if mutate:
            mutate(d)
        rc, out = run(d, *args)
        if rc != want_rc:
            failures.append('%s: exit %d, wanted %d\n%s' % (name, rc, want_rc, out))
            print('FAIL %s' % name)
        else:
            print('ok   %s' % name)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    # The tree as committed must be self-consistent on master.
    case('master tree is bound', 'master', None, 0)

    # Each way the binding can rot.
    case('datagrok back on :latest', 'master',
         lambda d: edit(d, 'localhost.docker-compose.yaml',
                        'datagrok/datagrok:${DATAGROK_VERSION:-bleeding-edge}',
                        'datagrok/datagrok:latest'), 1)
    case('marker names a release on master', 'master',
         lambda d: edit(d, 'localhost.docker-compose.yaml',
                        'x-datagrok-release: bleeding-edge',
                        'x-datagrok-release: 1.27.7'), 1)
    case('marker deleted', 'master',
         lambda d: _strip_marker(d, 'localhost.docker-compose.yaml'), 1)
    case('companion pinned to a release on master', 'master',
         lambda d: edit(d, 'localhost.docker-compose.yaml',
                        'GROK_PIPE_VERSION:-bleeding-edge',
                        'GROK_PIPE_VERSION:-1.22.2'), 1)
    case('install scripts disagree', 'master',
         lambda d: edit(d, 'datagrok-install-local.cmd',
                        'set datagrok_default_version=',
                        'set datagrok_default_version=9.9.9 REM '), 1)
    def set_default(value):
        def mutate(d):
            # Read the current default once -- editing the first file moves what _default sees.
            was = _default(d)
            if was == value:
                raise AssertionError('default is already %r, test proves nothing' % value)
            for f in ('datagrok-install-local.sh',
                      'datagrok-install-local-macos-silicon.sh'):
                edit(d, f, 'datagrok_default_version="%s"' % was,
                     'datagrok_default_version="%s"' % value)
            edit(d, 'datagrok-install-local.cmd',
                 'set datagrok_default_version=%s' % was,
                 'set datagrok_default_version=%s' % value)
        return mutate

    # The default is a Docker tag the scripts can resolve, not necessarily a version.
    case('install default X.Y.Z', 'master', set_default('1.27.7'), 0)
    case('install default bleeding-edge', 'master', set_default('bleeding-edge'), 0)
    case('install default is unresolvable', 'master', set_default('stable'), 1)

    # A release branch must not carry master's bleeding-edge binding.
    case('unstamped tree on a release branch', 'release/1.28.0', None, 1)

    # Stamping a train the branch does not own would pin new config to old images.
    case('stamp refuses the wrong train', 'master', None, 1,
         args=('stamp', '--version', '1.27.7', '--dry-run'))

    # The full cut-and-stamp flow.
    def stamp_it(d):
        rc, out = run(d, 'stamp', '--version', '1.28.0', '--set', 'datagrok=1.28.0',
                      '--set', 'grok_pipe=1.22.2', '--set', 'grok_connect=2.7.0',
                      '--set', 'grok_spawner=2.23.0',
                      '--set', 'jkg_python=1.1.0', '--set', 'jkg_r=1.1.0',
                      '--set', 'jkg_octave=1.1.0', '--set', 'jkg_julia=1.1.0',
                      '--set', 'jkg_nodejs=1.1.0',
                      '--allow-floating', 'grok_connect_adbc')
        if rc != 0:
            raise AssertionError('stamp failed:\n%s' % out)
        t = io.open(os.path.join(d, 'localhost.docker-compose.yaml'),
                    encoding='utf-8', newline='').read()
        assert 'x-datagrok-release: 1.28.0' in t, 'marker not stamped'
        assert '${DATAGROK_VERSION:-1.28.0}' in t, 'core tag not pinned'
        assert '${GROK_PIPE_VERSION:-1.22.2}' in t, 'companion tag not pinned'
        # A service left floating must not ride the default profile.
        m = re.search(r'^[ \t]*grok_connect_adbc:.*?profiles:[ \t]*\[([^\]]*)\]',
                      t, re.M | re.S)
        assert m and 'all' not in m.group(1), 'floating adbc still in the all profile'
    case('cut + stamp release/1.28.0', 'release/1.28.0', stamp_it, 0)

    print('')
    if failures:
        print('%d failure(s):\n' % len(failures))
        for f in failures:
            print(f)
        return 1
    print('all tests passed')
    return 0


def _strip_marker(d, fn):
    p = os.path.join(d, fn)
    t = io.open(p, encoding='utf-8', newline='').read()
    io.open(p, 'w', encoding='utf-8', newline='').write(
        re.sub(r'^x-datagrok-release:.*\r?\n', '', t, flags=re.M))


def _default(d):
    t = io.open(os.path.join(d, 'datagrok-install-local.sh'),
                encoding='utf-8', newline='').read()
    return re.search(r'datagrok_default_version="([^"]*)"', t).group(1)


if __name__ == '__main__':
    sys.exit(main())
