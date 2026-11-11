import {describe, it, expect, beforeAll, afterAll} from 'vitest';
import {spawnSync} from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';

const grokJs = path.join(__dirname, '..', 'grok.js');
const who = {GIT_AUTHOR_NAME: 't', GIT_AUTHOR_EMAIL: 't@t', GIT_COMMITTER_NAME: 't', GIT_COMMITTER_EMAIL: 't@t'};
// Submodules over file:// are refused by default since git 2.38.1.
const env = {...process.env, ...who, GIT_CONFIG_NOSYSTEM: '1', GIT_CONFIG_GLOBAL: os.devNull,
  GIT_CONFIG_COUNT: '1', GIT_CONFIG_KEY_0: 'protocol.file.allow', GIT_CONFIG_VALUE_0: 'always'};

function git(cwd: string, ...args: string[]): string {
  const r = spawnSync('git', args, {cwd, env, encoding: 'utf-8'});
  if (r.status !== 0)
    throw new Error(`git ${args.join(' ')}: ${r.stderr}`);
  return r.stdout.trim();
}

const grok = (cwd: string, ...args: string[]) => spawnSync(process.execPath, [grokJs, 'src', ...args], {cwd, env, encoding: 'utf-8'});

function emptyDirs(dir: string): string[] {
  const kids = fs.readdirSync(dir, {withFileTypes: true});
  const below = kids.filter((k) => k.isDirectory() && k.name !== '.git').flatMap((k) => emptyDirs(path.join(dir, k.name)));
  return kids.length ? below : [dir];
}

// A bare repo at <tmp>/<name>.git with the given files and gitlinks ({path: [repo, sha]}).
function repo(tmp: string, name: string, files: Record<string, string>, links: Record<string, [string, string]> = {}): string {
  const work = path.join(tmp, 'work', name);
  fs.mkdirSync(work, {recursive: true});
  git(work, 'init', '-q', '-b', 'main');
  const mods = Object.entries(links).map(([p, [r]]) => `[submodule "${p}"]\n\tpath = ${p}\n\turl = ../${r}.git\n\tbranch = main\n`);
  for (const [f, text] of Object.entries({...files, ...mods.length ? {'.gitmodules': mods.join('')} : {}})) {
    fs.mkdirSync(path.dirname(path.join(work, f)), {recursive: true});
    fs.writeFileSync(path.join(work, f), text);
  }
  git(work, 'add', '.');
  for (const [p, [, sha]] of Object.entries(links))
    git(work, 'update-index', '--add', '--cacheinfo', `160000,${sha},${p}`);
  git(work, 'commit', '-q', '-m', 'init');
  git(tmp, 'clone', '-q', '--bare', work, `${name}.git`);
  return git(work, 'rev-parse', 'HEAD');
}

describe('grok src help', () => {
  it.each(['--help', 'help', '', 'clone --help'])('`grok src %s` prints the help and exits 0', (args) => {
    const r = grok(os.tmpdir(), ...args.split(' ').filter(Boolean));
    expect(r.status, r.stderr).toBe(0);
    expect(r.stdout).toContain('Usage: grok src <command>');
  });

  it.each(['nope', 'at'])('`grok src %s` is a usage error', (cmd) => {
    expect(grok(os.tmpdir(), cmd).status).toBe(1);
  });

  it('clones the umbrella over HTTPS by default', () => {
    const r = grok(os.tmpdir(), 'clone', '--explain', '--dir', 'x');
    expect(r.stdout).toContain('https://github.com/datagrok-ai/datagrok.git');
  });
});

describe('grok src: no empty folders for children not cloned', () => {
  let tmp: string;
  let url: string;

  beforeAll(() => {
    tmp = fs.realpathSync(fs.mkdtempSync(path.join(os.tmpdir(), 'grok-src-')));
    const a = repo(tmp, 'a', {'a.txt': 'a'});
    const b = repo(tmp, 'b', {'b.txt': 'b'});
    const n = repo(tmp, 'n', {'n.txt': 'n'}, {x: ['a', a], 'deep/y': ['b', b]});
    repo(tmp, 'u', {'README.md': 'u', 'lib/notes.md': 'notes', '.gitprofiles': '[profile "default"]\n\tpath = a\n'},
      {a: ['a', a], 'lib/b': ['b', b], n: ['n', n], bad: ['missing', a]});
    url = path.join(tmp, 'u.git');
  });

  afterAll(() => fs.rmSync(tmp, {recursive: true, force: true}));

  it('clone takes profile "default" and leaves the other children out', () => {
    const r = grok(tmp, 'clone', '--url', url, '--dir', 'c');
    expect(r.status, r.stdout + r.stderr).toBe(0);
    const c = path.join(tmp, 'c');
    expect(fs.existsSync(path.join(c, 'a', 'a.txt'))).toBe(true);
    expect(fs.existsSync(path.join(c, 'lib', 'notes.md'))).toBe(true);
    expect(['lib/b', 'n', 'bad'].filter((p) => fs.existsSync(path.join(c, p)))).toEqual([]);
    expect(emptyDirs(c)).toEqual([]);
  });

  it('a plain git pull takes a new pin of a child that is not cloned', () => {
    const b = path.join(tmp, 'work', 'b');
    git(b, 'commit', '-q', '--allow-empty', '-m', 'next');
    git(b, 'push', '-q', path.join(tmp, 'b.git'), 'HEAD:main');
    const u = path.join(tmp, 'work', 'u');
    git(u, 'update-index', '--cacheinfo', `160000,${git(b, 'rev-parse', 'HEAD')},lib/b`);
    git(u, 'commit', '-q', '-m', 'bump');
    git(u, 'push', '-q', url, 'HEAD:main');
    const c = path.join(tmp, 'c');
    git(c, 'pull', '-q');
    expect(git(c, 'rev-parse', 'HEAD')).toBe(git(u, 'rev-parse', 'HEAD'));
    expect(fs.existsSync(path.join(c, 'lib', 'b'))).toBe(false);
  });

  it('a child that fails to clone leaves no folder', () => {
    const c = path.join(tmp, 'c');
    const r = grok(c, 'add', 'bad');
    expect(r.status).toBe(1);
    expect(r.stdout).toContain('Failed (see the git errors above): bad');
    expect(fs.existsSync(path.join(c, 'bad'))).toBe(false);
  });

  it('add of a nested path clones only that grandchild', () => {
    const c = path.join(tmp, 'c');
    expect(grok(c, 'add', 'n/deep/y').status).toBe(0);
    expect(fs.existsSync(path.join(c, 'n', 'deep', 'y', 'b.txt'))).toBe(true);
    expect(fs.existsSync(path.join(c, 'n', 'x'))).toBe(false);
    expect(emptyDirs(c)).toEqual([]);
    expect(grok(c, 'doctor').status).toBe(0);
  });

  it('wt add makes no folders for children left out of the worktree', () => {
    const c = path.join(tmp, 'c');
    const r = grok(c, 'wt', 'add', 'w', 'a');
    expect(r.status, r.stdout + r.stderr).toBe(0);
    const w = path.join(tmp, 'c-worktrees', 'w');
    expect(fs.existsSync(path.join(w, 'a', 'a.txt'))).toBe(true);
    expect(fs.existsSync(path.join(w, 'n'))).toBe(false);
    expect(emptyDirs(w)).toEqual([]);
  });

  it('doctor --fix moves a clone made without sparse patterns to them', () => {
    git(tmp, 'clone', '-q', url, 'old');
    const old = path.join(tmp, 'old');
    git(old, 'submodule', 'update', '-q', '--init', '--', 'a');
    git(path.join(old, 'a'), 'switch', '-q', 'main');
    expect(emptyDirs(old).length).toBeGreaterThan(0);
    const r = grok(old, 'doctor');
    expect(r.status).toBe(1);
    expect(r.stdout).toContain('(umbrella): empty folders for the children not cloned');
    const fix = grok(old, 'doctor', '--fix');
    expect(fix.status, fix.stdout).toBe(0);
    expect(emptyDirs(old)).toEqual([]);
    expect(fs.existsSync(path.join(old, 'a', 'a.txt'))).toBe(true);
  });
});
