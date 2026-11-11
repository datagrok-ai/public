import fs from 'fs';
import path from 'path';
import {spawnSync} from 'child_process';
import * as color from '../utils/color-utils';

interface SrcArgs {_: string[]; explain?: boolean; fix?: boolean; url?: string; dir?: string}
interface Child {name: string; path: string; url: string; branch?: string; dir: string}

const arity: Record<string, [number, number]> = {clone: [0, Infinity], add: [1, Infinity], sparse: [2, Infinity],
  switch: [0, 0], status: [0, 0], pull: [0, 0], at: [1, 1], bump: [0, Infinity], wt: [3, Infinity], doctor: [0, 0]};

const env: NodeJS.ProcessEnv = {...process.env, ...process.stdin.isTTY ? {} : {GIT_TERMINAL_PROMPT: '0'}};
for (const key of ['GIT_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE', 'GIT_PREFIX'])
  delete env[key];

let explain = false;

const quote = (s: string) => /^[\w@%+=:,./-]+$/.test(s) ? s : `'${s.replace(/'/g, `'\\''`)}'`;

function git(dir: string, args: string[]) {
  const r = spawnSync('git', ['-C', dir, ...args], {encoding: 'utf-8', env, maxBuffer: 64 * 1024 * 1024});
  return {ok: r.status === 0, fatal: r.status === 128, out: (r.stdout ?? '').trim(), err: (r.stderr ?? '').trim()};
}

function step(dir: string, args: string[]): boolean {
  const rel = path.relative(process.cwd(), dir);
  console.log(`$ git ${rel ? `-C ${quote(rel)} ` : ''}${args.map(quote).join(' ')}`);
  return explain || spawnSync('git', ['-C', dir, ...args], {stdio: 'inherit', env}).status === 0;
}

function findRoot(): string | undefined {
  let dir = git(process.cwd(), ['rev-parse', '--show-toplevel']);
  if (!dir.ok)
    throw new Error(dir.err);
  for (let up; (up = git(dir.out, ['rev-parse', '--show-superproject-working-tree'])).out;)
    dir = up;
  return git(dir.out, ['ls-files', '--', '.gitmodules']).out ? dir.out : undefined;
}

function children(dir: string): Child[] {
  if (!git(dir, ['ls-files', '--', '.gitmodules']).out)
    return [];
  const file = fs.existsSync(path.join(dir, '.gitmodules')) ? ['-f', '.gitmodules'] : ['--blob', ':./.gitmodules'];
  const r = git(dir, ['config', ...file, '--get-regexp', '^submodule\\..*\\.(path|url|branch)$']);
  if (r.fatal)
    throw new Error(`${path.join(dir, '.gitmodules')}: ${r.err}`);
  const mods = new Map<string, Record<string, string>>();
  for (const line of r.out.split('\n')) {
    const m = /^submodule\.(.+)\.(path|url|branch) (.*)$/.exec(line);
    if (m && !m[3].startsWith('-'))
      mods.set(m[1], {...mods.get(m[1]), [m[2]]: m[3]});
  }
  const links = new Set(git(dir, ['ls-files', '-z', '--stage']).out.split('\0')
    .filter((l) => l.startsWith('160000 ')).map((l) => l.slice(l.indexOf('\t') + 1)));
  return [...mods].filter(([, c]) => links.has(c.path) && c.url)
    .map(([name, c]) => ({name, path: c.path, url: c.url, branch: c.branch, dir: path.join(dir, c.path)}));
}

const cloned = (dir: string) => fs.existsSync(path.join(dir, '.git')) &&
  (fs.statSync(path.join(dir, '.git')).isFile() || !git(dir, ['rev-parse', '--show-prefix']).out);
const empty = (dir: string) => !fs.existsSync(dir) || !fs.readdirSync(dir).length;

function* repos(dir: string, c?: Child, owner = dir): Generator<{dir: string, branch?: string, c?: Child, owner: string}> {
  yield {dir, branch: c?.branch, c, owner};
  for (const k of children(dir).filter((x) => cloned(x.dir)))
    yield* repos(k.dir, k, dir);
}

function state(dir: string) {
  const r = git(dir, ['status', '--porcelain=v2', '--branch', '--untracked-files=no', '--ignore-submodules=all']);
  const head = /^# branch\.head (.*)$/m.exec(r.out)?.[1] ?? '';
  return {error: r.ok ? '' : r.err.split('\n')[0], branch: head === '(detached)' ? '' : head,
    upstream: /^# branch\.ab /m.test(r.out), ahead: Number(/^# branch\.ab \+(\d+)/m.exec(r.out)?.[1] ?? 0),
    changed: r.out.split('\n').filter((l) => l && !l.startsWith('#')).length};
}

const orphaned = (dir: string) => !!git(dir, ['rev-list', '-n1', 'HEAD', '--not', '--branches', '--remotes', '--tags']).out;

const label = (root: string, dir: string) => path.relative(root, dir) || '(umbrella)';
const gitC = (dir: string) => `git -C ${quote(path.relative(process.cwd(), dir) || '.')}`;

function locate(root: string, rel: string): Child[] | undefined {
  const chain: Child[] = [];
  for (let dir = root, rest = rel.replace(/\/+$/, ''); ;) {
    const c = children(dir).find((x) => rest === x.path || rest.startsWith(x.path + '/'));
    if (!c)
      return undefined;
    chain.push(c);
    if (rest === c.path)
      return chain;
    [dir, rest] = [c.dir, rest.slice(c.path.length + 1)];
  }
}

// Git makes an empty folder for every gitlink it checks out.
function fencing(dir: string, open: string[] = []) {
  const links = !cloned(dir) ? [] : git(dir, ['ls-tree', '-r', '-z', '--full-tree', 'HEAD']).out.split('\0')
    .filter((l) => l.startsWith('160000 ')).map((l) => l.slice(l.indexOf('\t') + 1));
  const shut = links.filter((p) => !open.includes(p) && !cloned(path.join(dir, p)));
  const exclude = (p: string) => `!/${p.replace(/[*?[\\]/g, '\\$&')}`;
  const sparse = git(dir, ['config', '--type=bool', 'core.sparseCheckout']).out === 'true';
  const now = sparse ? git(dir, ['sparse-checkout', 'list']).out : '/*';
  // Kept so a plain checkout of another commit makes no empty folder.
  const gone = now.split('\n').filter((l) => l.startsWith('!/') && !links.some((p) => exclude(p) === l));
  const want = ['/*', ...[...shut.map(exclude), ...gone].sort()];
  // Else a fetch of a new pin of a child with no folder fails: "Could not access submodule".
  const quiet = git(dir, ['config', 'fetch.recurseSubmodules']).out === 'false';
  return {sparse, want, quiet, ok: !links.length || (now === want.join('\n') && (quiet || !shut.length))};
}

function fence(dir: string, open: string[] = [], start = false): boolean {
  const f = fencing(dir, open);
  return f.ok || (!f.sparse && !start) || (step(dir, ['sparse-checkout', 'set', '--no-cone', ...f.want]) &&
    (f.quiet || step(dir, ['config', 'fetch.recurseSubmodules', 'false'])));
}

const modules = (owner: string, c: Child) => path.resolve(owner, git(owner, ['rev-parse', '--git-path', `modules/${c.name}`]).out);

function rollback(owner: string, c: Child, made?: string): void {
  if (cloned(c.dir) || empty(c.dir))
    step(owner, ['submodule', 'deinit', '--force', '--', c.path]);
  step(owner, ['config', `submodule.${c.name}.active`, 'false']);
  if (made)
    fs.rmSync(made, {recursive: true, force: true});
  fence(owner);
}

function cloneChild(root: string, owner: string, c: Child, failed: string[]): boolean {
  if (cloned(c.dir))
    return true;
  const gitdir = modules(owner, c);
  const made = fs.existsSync(gitdir) ? undefined : gitdir;
  if (!fence(owner, [c.path]) ||
    !step(owner, ['submodule', 'update', '--init', '--require-init', '--checkout', '--filter=blob:none', '--', c.path])) {
    rollback(owner, c, made);
    failed.push(label(root, c.dir));
    return false;
  }
  if (c.branch && !step(c.dir, ['switch', c.branch]))
    failed.push(label(root, c.dir));
  if (fs.existsSync(path.join(c.dir, '.gitmodules')))
    step(c.dir, ['config', 'push.recurseSubmodules', 'check']);
  fence(c.dir, [], true);
  return true;
}

function cloneTree(root: string, owner: string, c: Child, failed: string[]): void {
  const fresh = !cloned(c.dir);
  if (!cloneChild(root, owner, c, failed) || !git(c.dir, ['ls-files', '--', '.gitmodules']).out)
    return;
  const cur = fresh ? c.branch : git(c.dir, ['branch', '--show-current']).out;
  if (cur === c.branch || (!cur && git(c.dir, ['rev-parse', 'HEAD']).out === git(owner, ['rev-parse', `:${c.path}`]).out))
    children(c.dir).forEach((k) => cloneTree(root, c.dir, k, failed));
}

function addSpec(root: string, dir: string, spec: string, failed: string[]): boolean {
  const kids = children(dir);
  const parent = kids.find((c) => spec.startsWith(c.path + '/'));
  if (parent) {
    return !cloneChild(root, dir, parent, failed) || (explain && !cloned(parent.dir)) ||
      addSpec(root, parent.dir, spec.slice(parent.path.length + 1), failed);
  }
  const listed = git(dir, ['ls-files', '-z', '--', spec]).out.split('\0');
  const targets = kids.filter((c) => listed.includes(c.path));
  targets.forEach((c) => cloneTree(root, dir, c, failed));
  return targets.length > 0;
}

function add(root: string, sel: string[]): boolean {
  const specs = (sel.length ? sel : ['default']).flatMap((s) => {
    const profile = git(root, ['config', '-f', '.gitprofiles', '--get-all', `profile.${s}.path`]);
    if (profile.fatal)
      throw new Error(`${path.join(root, '.gitprofiles')}: ${profile.err}`);
    return profile.ok ? profile.out.split('\n') : sel.length ? [s] : ['.'];
  });
  const failed: string[] = [];
  const unknown = specs.map((s) => s.replace(/^\.\//, '').replace(/\/+$/, '')).filter((s) => !addSpec(root, root, s, failed));
  unknown.forEach((s) => color.error(`${s}: no profile or submodule by that name`));
  if (failed.length)
    color.error(`Failed (see the git errors above): ${failed.join(', ')}`);
  return !unknown.length && !failed.length;
}

function clone(sel: string[], url: string, dir: string): boolean {
  if (!step(process.cwd(), ['clone', '-c', 'push.recurseSubmodules=check', '--filter=blob:none', '--no-checkout', '--', url, dir]))
    return false;
  if (!explain)
    return fence(dir, [], true) && step(dir, ['checkout']) && add(dir, sel);
  console.log(`$ ${gitC(dir)} sparse-checkout set --no-cone '/*' '!/<child>'...   # every child`);
  console.log(`$ ${gitC(dir)} checkout`);
  console.log(`# then, in ${dir}: grok src add ${sel.join(' ') || '(profile "default", else every child)'}`);
  return true;
}

function sparse(root: string, rel: string, globs: string[]): boolean {
  const chain = locate(root, rel);
  if (!chain)
    return (color.error(`${rel}: not a submodule, or its parent is not cloned`), false);
  const c = chain[chain.length - 1];
  const owner = chain.length > 1 ? chain[chain.length - 2].dir : root;
  const absorb = ['submodule', 'absorbgitdirs', '--', c.path];
  const set = ['sparse-checkout', 'set', '--no-cone', ...globs];
  if (cloned(c.dir))
    return step(owner, absorb) && step(c.dir, set);
  const gitdir = modules(owner, c);
  if (fs.existsSync(gitdir))
    return (color.error(`${gitdir} exists from an earlier clone: grok src add ${rel}, then sparse again`), false);
  const ok = fence(owner, [c.path]) && step(owner, ['submodule', 'init', '--', c.path]) &&
    step(owner, ['clone', '--filter=blob:none', '--no-checkout', '--',
      git(owner, ['config', `submodule.${c.name}.url`]).out || '<url>', c.path]) &&
    step(owner, absorb) && step(c.dir, set) &&
    step(c.dir, ['checkout', c.branch ?? git(owner, ['rev-parse', `:${c.path}`]).out, '--']);
  if (!ok)
    rollback(owner, c, gitdir);
  return ok;
}

function switchAll(root: string): boolean {
  let ok = true;
  for (const {dir, branch} of repos(root)) {
    const s = state(dir);
    if (!branch || s.branch === branch)
      continue;
    if (s.error || s.changed || (!s.branch && orphaned(dir))) {
      color.warn(`${label(root, dir)}: ${s.error || 'uncommitted changes or commits on no branch'}, not switched to ${branch}`);
      ok = false;
    } else
      ok = step(dir, ['switch', branch]) && ok;
  }
  return ok;
}

function status(root: string): boolean {
  const head = ['REPO', 'BRANCH', 'UNPUSHED', 'CHANGED'];
  const rows = [...repos(root)].map(({dir}) => {
    const s = state(dir);
    const unpushed = !s.branch ? '' : !s.upstream ? 'no upstream' : s.ahead ? `${s.ahead}` : '';
    return [label(root, dir), s.error ? '(error: grok src doctor)' : s.branch || '(detached)', unpushed, s.changed ? `${s.changed}` : ''];
  });
  const widths = head.map((_, i) => Math.max(...[head, ...rows].map((r) => r[i].length)));
  for (const r of [head, ...rows])
    console.log(r.map((x, i) => x.padEnd(widths[i])).join('  ').trimEnd());
  return !rows.some((r) => r[1].startsWith('(error: '));
}

function pull(root: string): boolean {
  const skipped: string[] = [];
  const failed: string[] = [];
  for (const {dir, branch} of repos(root)) {
    const s = state(dir);
    if (s.error)
      failed.push(label(root, dir));
    else if (!s.branch || !s.upstream) {
      if (branch || dir === root)
        skipped.push(label(root, dir));
    } else if (!step(dir, ['pull', '--ff-only', '--no-recurse-submodules']) || !fence(dir))
      failed.push(label(root, dir));
  }
  if (skipped.length)
    color.warn(`Skipped (detached or no upstream): ${skipped.join(', ')}`);
  if (failed.length)
    color.error(`Not pulled: ${failed.join(', ')}`);
  return !failed.length;
}

function at(root: string, ref: string): boolean {
  const all = [...repos(root)];
  const loose = all.flatMap(({c, owner}) => c && !git(owner, ['config', `submodule.${c.name}.url`]).ok ?
    [`${gitC(owner)} submodule init -- ${quote(c.path)}`] : []);
  if (loose.length)
    return (color.error(`Cloned but not initialized, so git would leave them off the pins: ${loose.join('; ')}`), false);
  const blocked = all.filter(({dir}) => {
    const s = state(dir);
    return s.error || s.changed || (!s.branch && orphaned(dir)) ||
      !git(dir, ['diff', '--cached', '--quiet', '--ignore-submodules=none']).ok;
  }).map(({dir}) => label(root, dir));
  if (blocked.length)
    return (color.error(`Uncommitted changes or commits on no branch in ${blocked.join(', ')}: commit, stash or branch first`), false);
  return step(root, ['checkout', ref, '--']) && step(root, ['submodule', 'update', '--recursive', '--checkout']) &&
    [...repos(root)].map(({dir}) => fence(dir)).every(Boolean);
}

function bump(root: string, rels: string[]): boolean {
  const failed: string[] = [];
  const owners = new Map<string, {ok: boolean, links: {path: string, tip: string}[]}>();
  let checked = 0;
  const skipped: string[] = [];
  const visit = (c: Child, owner: string) => {
    let o = owners.get(owner);
    if (!o) {
      const staged = !git(owner, ['diff', '--cached', '--quiet', '--ignore-submodules=none']).ok;
      const detached = !state(owner).branch;
      if (staged || (detached && rels.length))
        failed.push(`${label(root, owner)} (${staged ? 'staged changes' : 'not on a branch'})`);
      else if (detached)
        skipped.push(label(root, owner));
      owners.set(owner, o = {ok: !staged && !detached, links: []});
    }
    const name = label(root, c.dir);
    const branch = c.branch === '.' ? state(owner).branch : c.branch;
    if (!o.ok || !branch) {
      if (o.ok && rels.length)
        failed.push(`${name} (no branch in .gitmodules)`);
      return;
    }
    if (git(c.dir, ['config', 'remote.origin.url']).out !== git(owner, ['config', `submodule.${c.name}.url`]).out)
      return failed.push(`${name} (origin differs from .gitmodules: grok src doctor --fix)`);
    checked++;
    const ref = `refs/remotes/origin/${branch}`;
    if (!step(c.dir, ['fetch', '--no-recurse-submodules', 'origin', `+refs/heads/${branch}:${ref}`]))
      return failed.push(name);
    const tip = explain ? git(c.dir, ['ls-remote', 'origin', `refs/heads/${branch}`]).out.split('\t')[0] :
      git(c.dir, ['rev-parse', '--verify', '-q', `${ref}^{commit}`]).out;
    const pin = git(owner, ['rev-parse', `HEAD:${c.path}`]).out;
    if (!tip)
      return failed.push(`${name} (origin has no ${branch})`);
    if (tip !== pin && !explain && !git(c.dir, ['merge-base', '--is-ancestor', pin, tip]).ok)
      return failed.push(`${name} (the pin is not on origin/${branch}; move it by hand)`);
    const s = state(c.dir);
    const head = git(c.dir, ['rev-parse', 'HEAD']).out;
    const moved = head === tip || (!s.changed && (s.branch === branch ? step(c.dir, ['merge', '--ff-only', ref]) :
      !s.branch && git(c.dir, ['merge-base', '--is-ancestor', head, tip]).ok && step(c.dir, ['checkout', '--detach', tip])));
    if (tip === pin)
      return;
    if (!moved)
      color.warn(`${name}: checkout left as it is (local changes, or not on ${branch} or the pin); pinned to origin/${branch}`);
    o.links.push({path: c.path, tip});
  };
  const wanted = [...new Set(rels.map((r) => r.replace(/\\/g, '/').replace(/^(\.\/)+/, '').replace(/\/+$/, '')))];
  for (const rel of wanted.sort((a, b) => a.split('/').length - b.split('/').length)) {
    const chain = locate(root, rel);
    const c = chain?.[chain.length - 1];
    if (!chain || !c || !cloned(c.dir))
      failed.push(chain ? `${rel} (not cloned: grok src add ${rel})` : `${rel} (not a submodule)`);
    else
      visit(c, chain.length > 1 ? chain[chain.length - 2].dir : root);
  }
  if (!rels.length) {
    for (const {c, owner} of repos(root)) {
      if (c)
        visit(c, owner);
    }
  }
  let bumped = false;
  for (const [owner, {links}] of owners) {
    if (!links.length)
      continue;
    bumped = true;
    const paths = links.map((l) => l.path);
    if (step(owner, ['update-index', ...links.flatMap((l) => ['--cacheinfo', `160000,${l.tip},${l.path}`])]) &&
      step(owner, ['commit', '-m', `Bump ${paths.join(', ')}`])) {
      if (owner !== root)
        console.log(`${label(root, owner)}: push it, then grok src bump ${path.relative(root, owner).split(path.sep).join('/')}`);
    } else {
      if (!explain)
        git(owner, ['reset', '-q', '--', ...paths]);
      failed.push(label(root, owner));
    }
  }
  if (skipped.length)
    color.warn(`Not on a branch, so their children were not bumped: ${skipped.join(', ')}`);
  else if (!bumped && !failed.length)
    console.log(checked ? 'Every pin is at the tip of its branch.' : 'No child with a branch to bump.');
  if (failed.length)
    color.error(`Not bumped: ${failed.join(', ')}`);
  return !failed.length;
}

function wtAdd(root: string, name: string, rels: string[], destArg?: string): boolean {
  const dest = path.resolve(destArg ?? `${root}-worktrees/${name}`);
  const dirs = new Set([root]);
  for (const rel of rels) {
    const chain = locate(root, rel);
    if (!chain || !cloned(chain[chain.length - 1].dir))
      return (color.error(`${rel}: not a cloned submodule`), false);
    chain.forEach((c) => dirs.add(c.dir));
  }
  const todo = [...dirs].map((dir) => ({dir, to: path.join(dest, path.relative(root, dir))}));
  for (const {dir, to} of todo) {
    const elsewhere = git(dir, ['for-each-ref', '--format=%(worktreepath)', `refs/heads/${name}`]).out;
    if (cloned(to) ? git(to, ['branch', '--show-current']).out !== name : elsewhere) {
      const where = elsewhere && (git(dir, [`--git-dir=${elsewhere}`, 'rev-parse', '--show-toplevel']).out || elsewhere);
      return (color.error(cloned(to) ? `${to} exists and is not on ${name}` :
        `${label(root, dir)}: ${name} is checked out in ${where} (if that is gone: ${gitC(dir)} worktree prune)`), false);
    }
  }
  for (const {dir, to} of todo.filter((t) => !cloned(t.to))) {
    const local = git(dir, ['show-ref', '--verify', '--quiet', `refs/heads/${name}`]).ok;
    const remote = !local && git(dir, ['show-ref', '--verify', '--quiet', `refs/remotes/origin/${name}`]).ok;
    const start = local ? ['--', to, name] : remote ? ['--track', '-b', name, '--', to, `origin/${name}`] : ['-b', name, '--', to];
    if (!step(dir, ['worktree', 'add', ...start]))
      return false;
    const own = explain ? '<worktree gitdir>/config.worktree' :
      git(to, ['--work-tree=.', 'rev-parse', '--path-format=absolute', '--git-path', 'config.worktree']).out;
    if (dir !== root && !step(dir, ['config', '-f', own, 'core.worktree', to]))
      return false;
    if (!fence(to, [...dirs].map((d) => path.relative(dir, d).split(path.sep).join('/')), true))
      return false;
  }
  console.log(`Worktrees: ${dest}`);
  return true;
}

function doctor(root: string, fix: boolean): boolean {
  const problems: string[] = [];
  const synced: boolean[] = [];
  if (process.platform === 'win32' && git(root, ['config', '--type=bool', 'core.longpaths']).out !== 'true')
    problems.push('core.longpaths is not set: git config --global core.longpaths true');
  for (const {dir, branch, c, owner} of repos(root)) {
    const s = state(dir);
    const name = label(root, dir);
    const at = gitC(dir);
    if (s.error) {
      problems.push(`${name}: ${s.error}`);
      continue;
    }
    const orphan = !s.branch && orphaned(dir);
    if (!s.branch && (branch || orphan))
      problems.push(`${name}: detached HEAD: ${orphan ? `${at} switch -c <name>` : 'grok src switch'}`);
    if (s.ahead)
      problems.push(`${name}: ${s.ahead} unpushed commit(s): ${at} push`);
    const local = s.branch && !s.upstream ? Number(git(dir, ['rev-list', '--count', 'HEAD', '--not', '--remotes']).out) : 0;
    if (local)
      problems.push(`${name}: ${local} commit(s) on no remote: ${at} push -u origin ${s.branch}`);
    const f = fencing(dir);
    if (!f.ok) {
      problems.push(`${name}: ${f.sparse ? 'sparse checkout does not match the cloned children' :
        'empty folders for the children not cloned'}: grok src doctor --fix`);
      synced.push(!fix || fence(dir, [], true));
    }
    if (!c)
      continue;
    const base = git(owner, ['config', 'remote.origin.url']).out.replace(/\/+$/, '');
    const url = c.url.startsWith('../') ? base && base.replace(/[^/:]+$/, '') + c.url.slice(3) : c.url;
    if (!url) {
      problems.push(`${name}: ${label(root, owner)} has no origin remote`);
      continue;
    }
    const known = git(owner, ['config', `submodule.${c.name}.url`]);
    if (known.out !== url || git(dir, ['config', 'remote.origin.url']).out !== url) {
      problems.push(`${name}: URL differs from .gitmodules (${url}): grok src doctor --fix`);
      synced.push(!fix || ((known.ok || step(owner, ['submodule', 'init', '--', c.path])) &&
        step(owner, ['config', `submodule.${c.name}.url`, url]) && step(dir, ['remote', 'set-url', 'origin', url])));
    }
  }
  if (fix && !explain && synced.length && synced.every(Boolean))
    return doctor(root, false);
  problems.forEach((p) => console.log(p));
  if (!problems.length)
    color.success('No problems found.');
  return !problems.length || (fix && explain && synced.length > 0);
}

function run(cmd: string, rest: string[], argv: SrcArgs): boolean {
  if (cmd === 'clone') {
    const url = argv.url ?? 'https://github.com/datagrok-ai/datagrok.git';
    return clone(rest, url, path.resolve(argv.dir ?? path.basename(url).replace(/\.git$/, '')));
  }
  const root = findRoot();
  if (!root)
    return (color.error('Not inside an umbrella checkout (a repo with .gitmodules)'), false);
  const [gitDir, commonDir] = git(root, ['rev-parse', '--path-format=absolute', '--git-dir', '--git-common-dir']).out.split('\n');
  if (gitDir !== commonDir && (['add', 'sparse', 'switch', 'at'].includes(cmd) || (cmd === 'doctor' && argv.fix)))
    return (color.error(`grok src ${cmd} changes what the main checkout shares (config, branches): run it in the main checkout`), false);
  switch (cmd) {
    case 'add': return add(root, rest);
    case 'sparse': return sparse(root, rest[0], rest.slice(1));
    case 'switch': return switchAll(root);
    case 'status': return status(root);
    case 'pull': return pull(root);
    case 'at': return at(root, rest[0]);
    case 'bump': return bump(root, rest);
    case 'wt': return wtAdd(root, rest[1], rest.slice(2), argv.dir);
    default: return doctor(root, !!argv.fix);
  }
}

export async function src(argv: SrcArgs): Promise<boolean> {
  const [cmd, ...rest] = argv._.slice(1);
  const options = Object.keys(argv).filter((k) => !['_', 'explain', 'fix', 'url', 'dir'].includes(k));
  if (options.length || [argv.url, argv.dir].some((v) => v !== undefined && (typeof v !== 'string' || !v))) {
    color.error(options.length ? `Unknown option --${options[0]}` : '--url and --dir take one value');
    return false;
  }
  if (!cmd || cmd === 'help') {
    console.log((await import('./help')).help.src);
    return true;
  }
  if (!Object.hasOwn(arity, cmd) || rest.length < arity[cmd][0] || rest.length > arity[cmd][1] || (cmd === 'wt' && rest[0] !== 'add')) {
    color.error(`Unknown command or wrong arguments: grok src ${cmd}`);
    return false;
  }
  const [major, minor] = (/(\d+)\.(\d+)/.exec(git('.', ['--version']).out) ?? []).slice(1).map(Number);
  if (!(major * 100 + minor >= 236)) {
    color.error('grok src needs git 2.36 or newer');
    process.exitCode = 1;
    return true;
  }
  explain = !!argv.explain;
  try {
    if (!run(cmd, rest, argv))
      process.exitCode = 1;
  } catch (e) {
    color.error(e instanceof Error ? e.message : String(e));
    process.exitCode = 1;
  }
  return true;
}
