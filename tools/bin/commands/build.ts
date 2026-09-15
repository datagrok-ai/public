import fs from 'fs';
import path from 'path';
import * as readline from 'readline';
import {spawnSync} from 'child_process';
import * as utils from '../utils/utils';
import * as color from '../utils/color-utils';
import {toolchain, copyLibraryAssets, ensureLibraryExports} from '../utils/toolchain';

interface BuildArgs {
  _: string[];
  all?: boolean;
  recursive?: boolean;
  affected?: boolean;
  typecheck?: boolean;
  filter?: string;
  parallel?: number;
  verbose?: boolean;
  force?: boolean;
  local?: boolean;
  'skip-check'?: boolean;
  [key: string]: any;
}

/** Walks up from `dir` to the pnpm workspace root (the directory holding pnpm-workspace.yaml). */
export function findWorkspaceRoot(dir: string): string | null {
  let d = path.resolve(dir);
  for (;;) {
    if (fs.existsSync(path.join(d, 'pnpm-workspace.yaml')))
      return d;
    const parent = path.dirname(d);
    if (parent === d)
      return null;
    d = parent;
  }
}

/**
 * `grok build` in a workspace drives Turborepo: the package in cwd and everything it depends on,
 * in dependency order, cached (`--all`, `--affected`, `--typecheck`, `--filter`). Inside a Turborepo
 * task (TURBO_HASH is set), in a standalone package, or with `--local`, it builds just the package in
 * cwd: an rspack bundle plus `grok check` for a plugin (src/package.ts or rspack.config.js), a
 * TypeScript emit into dist/ for a library.
 */
export async function build(args: BuildArgs): Promise<boolean> {
  const cwd = process.cwd();
  const root = findWorkspaceRoot(cwd);
  if (args.local || process.env.TURBO_HASH || !root)
    return buildLocal(cwd, args);
  return buildWithTurbo(cwd, root, args);
}

async function buildWithTurbo(cwd: string, root: string, args: BuildArgs): Promise<boolean> {
  const filters: string[] = [];
  if (args.affected)
    filters.push('--filter=...[origin/master]');
  else if (!(args.all || args.recursive)) {
    const pkgPath = path.join(cwd, 'package.json');
    if (!fs.existsSync(pkgPath)) {
      color.error('Not a package directory (no package.json). Use --all or --affected from anywhere in the workspace.');
      return false;
    }
    filters.push(`--filter=${JSON.parse(fs.readFileSync(pkgPath, 'utf8')).name}...`);
  }
  if (args.filter)
    filters.push(`--filter=${args.filter}`);

  const tasks = ['build'];
  if (args.typecheck)
    tasks.push('typecheck');
  const cmd = ['exec', 'turbo', 'run', ...tasks, ...filters, `--concurrency=${args.parallel ?? 3}`, '--continue'];
  if (!args.verbose)
    cmd.push('--output-logs=errors-only');
  if (args.force)
    cmd.push('--force');

  const r = spawnSync('pnpm', cmd, {stdio: 'inherit', cwd: root, shell: true, env: {...process.env, TURBO_TELEMETRY_DISABLED: '1'}});
  return r.status === 0;
}

async function buildLocal(dir: string, args: BuildArgs): Promise<boolean> {
  if (!fs.existsSync(path.join(dir, 'package.json'))) {
    color.error('Not a package directory (no package.json).');
    return false;
  }
  const isPlugin = fs.existsSync(path.join(dir, 'rspack.config.js')) ||
    ['ts', 'js'].some((e) => fs.existsSync(path.join(dir, 'src', `package.${e}`)));
  return isPlugin ? bundle(dir, args) : emitLibrary(dir);
}

async function bundle(dir: string, args: BuildArgs): Promise<boolean> {
  const {bundler, rspack} = toolchain(dir).buildConfig;
  // --key=value arguments reach a config exported as a function, e.g. `grok build --only=browser`
  const env: Record<string, any> = {};
  for (const [k, v] of Object.entries(args))
    if (!['_', 'local', 'verbose', 'skip-check', 'typecheck', 'all', 'affected', 'filter', 'parallel', 'force'].includes(k))
      env[k] = v;
  const cfgPath = path.join(dir, 'rspack.config.js');
  let config = fs.existsSync(cfgPath) ? require(cfgPath) : bundler({dir});
  if (typeof config === 'function')
    config = config(env);
  if (config && typeof config.then === 'function')
    config = await config;

  if (!args['skip-check'] && fs.existsSync(path.join(dir, 'src', 'package.ts'))) {
    const {check} = require('./check');
    if (check({_: ['check'], soft: true, 'no-exit': true}) === false) {
      color.error('grok check failed');
      return false;
    }
  }

  const t0 = Date.now();
  const ok = await new Promise<boolean>((resolve, reject) => {
    rspack(config, (err: any, stats: any) => {
      if (err) return reject(err);
      const text = stats.toString({colors: process.stdout.isTTY, preset: 'errors-warnings'});
      if (text.trim()) console.log(text);
      resolve(!stats.hasErrors());
    });
  });
  const dist = path.join(dir, 'dist', 'package.js');
  const size = fs.existsSync(dist) ? `${(fs.statSync(dist).size / 1024).toFixed(0)} KB` : '';
  console.log(`${ok ? 'bundled' : 'FAILED'} ${path.basename(dir)} in ${((Date.now() - t0) / 1000).toFixed(1)}s ${size}`);
  return ok;
}

/** A library: `tsc -p tsconfig.json` into dist/, then the css/wasm/json assets its sources import. */
function emitLibrary(dir: string): boolean {
  ensureLibraryExports(dir);
  const r = spawnSync(process.execPath, [toolchain(dir).tsc, '-p', 'tsconfig.json'], {stdio: 'inherit', cwd: dir});
  if (r.status !== 0)
    return false;
  copyLibraryAssets(dir);
  return true;
}

export interface PackageInfo {
  dir: string;
  name: string;
  friendlyName: string;
  version: string;
  packageJson: any;
}

/** Packages directly under `baseDir` (used by `grok testall`). */
export function discoverPackages(baseDir: string): PackageInfo[] {
  const packages: PackageInfo[] = [];
  for (const entry of fs.readdirSync(baseDir)) {
    if (entry.startsWith('.') || entry === 'node_modules')
      continue;
    const dir = path.join(baseDir, entry);
    const packageJsonPath = path.join(dir, 'package.json');
    if (!fs.existsSync(packageJsonPath))
      continue;
    try {
      const packageJson = JSON.parse(fs.readFileSync(packageJsonPath, 'utf-8'));
      packages.push({dir, name: packageJson.name || entry, friendlyName: packageJson.friendlyName || packageJson.name || entry,
        version: packageJson.version || '0.0.0', packageJson});
    }
    catch (_) { /* not a package */ }
  }
  return packages.sort((a, b) => a.friendlyName.localeCompare(b.friendlyName));
}

/** `field:regex && field:regex` filter over package.json fields. */
export function applyFilter(packages: PackageInfo[], filterStr: string): PackageInfo[] {
  const conditions = filterStr.split('&&').map((s) => s.trim()).map((cond) => {
    const colonIdx = cond.indexOf(':');
    return colonIdx === -1 ? {field: cond, pattern: new RegExp('.')} :
      {field: cond.substring(0, colonIdx).trim(), pattern: new RegExp(cond.substring(colonIdx + 1).trim())};
  });
  return packages.filter((pkg) => conditions.every((cond) => {
    const value = getNestedValue(pkg.packageJson, cond.field);
    return value !== undefined && cond.pattern.test(String(value));
  }));
}

export function getNestedValue(obj: any, path: string): any {
  let current = obj;
  for (const part of path.split('.')) {
    if (current == null || typeof current !== 'object')
      return undefined;
    current = current[part];
  }
  return current;
}

export function confirm(message: string): Promise<boolean> {
  const rl = readline.createInterface({input: process.stdin, output: process.stdout});
  return new Promise((resolve) => {
    rl.question(`${message} [Y/n] `, (answer) => {
      rl.close();
      resolve(answer.trim().toLowerCase() !== 'n');
    });
  });
}
