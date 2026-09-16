import fs from 'fs';
import path from 'path';
import {spawnSync} from 'child_process';
import * as color from '../utils/color-utils';
import {findWorkspaceRoot} from './build';

interface SetupArgs {
  _: string[];
  check?: boolean;
  global?: boolean;
  verbose?: boolean;
}

/**
 * `grok setup`: one command to get (or keep) a public/ checkout ready to build.
 *   1. Node 20+ (22 recommended); pnpm through corepack, at the version package.json pins.
 *   2. Removes per-package node_modules left by the npm era and stray package-lock.json files.
 *   3. `pnpm install` at the workspace root.
 *   4. Reports a global `grok` older than the workspace one (`--global` updates it).
 * `--check` only reports.
 */
export async function setup(args: SetupArgs): Promise<boolean> {
  const root = findWorkspaceRoot(process.cwd());
  if (!root) {
    color.error('Not inside the public/ workspace (no pnpm-workspace.yaml above the current directory).');
    return false;
  }
  const check = !!args.check;
  const rootPkg = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
  const pinnedPnpm = (rootPkg.packageManager || 'pnpm@10').split('@')[1];
  let ok = true;

  // 1. Node and pnpm
  const nodeMajor = parseInt(process.versions.node.split('.')[0], 10);
  if (nodeMajor < 20) {
    color.error(`Node ${process.versions.node}: 20 or later is required (22 recommended).`);
    return false;
  }
  color.info(`Node ${process.versions.node}`);
  let pnpmVersion = run('pnpm', ['--version']);
  if (pnpmVersion !== pinnedPnpm) {
    color.warn(`pnpm ${pnpmVersion || 'not found'}; the workspace pins ${pinnedPnpm}.`);
    if (!check) {
      if (run('corepack', ['--version']) === null)
        color.warn('corepack not found: `npm install -g corepack`, then rerun `grok setup`.');
      else {
        spawnSync('corepack', ['enable'], {stdio: 'inherit', shell: true});
        spawnSync('corepack', ['prepare', `pnpm@${pinnedPnpm}`, '--activate'], {stdio: 'inherit', shell: true});
        pnpmVersion = run('pnpm', ['--version']);
      }
    }
  }
  if (pnpmVersion === pinnedPnpm)
    color.info(`pnpm ${pnpmVersion}`);
  else
    ok = false;

  // 2. leftovers from the npm era
  const projectDirs = ['js-api', 'tools', 'build-config'].map((d) => path.join(root, d));
  for (const group of ['packages', 'libraries'])
    for (const d of fs.readdirSync(path.join(root, group)))
      projectDirs.push(path.join(root, group, d));
  const staleModules = projectDirs.filter((d) => isLegacyNodeModules(path.join(d, 'node_modules')));
  const staleLocks = projectDirs.map((d) => path.join(d, 'package-lock.json')).filter((f) => fs.existsSync(f));
  if (staleModules.length || staleLocks.length) {
    color.warn(`${staleModules.length} per-package node_modules from npm and ${staleLocks.length} package-lock.json file(s)` +
      (check ? ' would be removed' : ': removing'));
    if (!check) {
      for (const d of staleModules) fs.rmSync(path.join(d, 'node_modules'), {recursive: true, force: true});
      for (const f of staleLocks) fs.rmSync(f, {force: true});
    }
  }
  else
    color.info('no npm-era leftovers');

  // 3. install
  if (!check && pnpmVersion === pinnedPnpm) {
    color.log('pnpm install ...');
    const r = spawnSync('pnpm', ['install', '--frozen-lockfile'], {stdio: 'inherit', cwd: root, shell: true});
    if (r.status !== 0) {
      color.warn('frozen install failed (lockfile out of date?): retrying without --frozen-lockfile');
      if (spawnSync('pnpm', ['install'], {stdio: 'inherit', cwd: root, shell: true}).status !== 0) {
        color.error('pnpm install failed');
        return false;
      }
    }
  }

  // 4. the global grok
  const workspaceGrok = JSON.parse(fs.readFileSync(path.join(root, 'tools', 'package.json'), 'utf8')).version;
  // datagrok-tools before 6.6 has no --version and prints its usage instead
  const raw = run('grok', ['--version']);
  const globalGrok = raw === null ? null : /^\d+\.\d+\.\d+/.test(raw) ? raw : 'older than 6.6';
  if (globalGrok && globalGrok !== workspaceGrok) {
    color.warn(`global grok ${globalGrok}; the workspace has ${workspaceGrok}` +
      (args.global && !check ? ': updating' : ' (run `npm install -g datagrok-tools@' + workspaceGrok + '`, or `grok setup --global`)'));
    if (args.global && !check)
      spawnSync('npm', ['install', '-g', `datagrok-tools@${workspaceGrok}`], {stdio: 'inherit', shell: true});
  }
  else if (globalGrok)
    color.info(`grok ${globalGrok}`);

  if (ok && !check)
    color.success('Ready: `grok build` in any package, `grok build --all` for everything.');
  return ok;
}

/** npm's node_modules holds real directories; pnpm's per-package one holds only symlinks (and .bin). */
function isLegacyNodeModules(dir: string): boolean {
  if (!fs.existsSync(dir))
    return false;
  for (const e of fs.readdirSync(dir, {withFileTypes: true})) {
    if (e.name === '.bin' || e.name === '.cache' || e.name.startsWith('.'))
      continue;
    const p = path.join(dir, e.name);
    if (e.isSymbolicLink())
      continue;
    if (e.isDirectory() && e.name.startsWith('@')) {
      if (fs.readdirSync(p, {withFileTypes: true}).some((s) => s.isDirectory() && !s.isSymbolicLink()))
        return true;
      continue;
    }
    if (e.isDirectory())
      return true;
  }
  return false;
}

function run(cmd: string, args: string[]): string | null {
  const r = spawnSync(cmd, args, {encoding: 'utf8', shell: true});
  return r.status === 0 ? r.stdout.trim().split(/\r?\n/).pop() || '' : null;
}
