/* A bdd project: where its features, bindings and generated specs live, and which binding tiers it
   loads. A package keeps everything under `<package>/bdd/`; the library itself is a project too
   (its own smoke features under `features/`). Every project sees the base tiers (`common`,
   `platform`), the tiers its `bdd.config.json` names, and its own `bindings/`. */
import {existsSync, readdirSync, readFileSync, statSync} from 'node:fs';
import {basename, dirname, join, relative, resolve, sep} from 'node:path';
import {fileURLToPath} from 'node:url';

export const PACKAGE_NAME = '@datagrok-libraries/bdd';
/** The directory holding this build's `bindings/`: `dist/` when running compiled, the library root
 * when running from source. */
export const LIB_DIR = resolve(dirname(fileURLToPath(import.meta.url)), '..');
export const LIB_ROOT = basename(LIB_DIR) === 'dist' ? dirname(LIB_DIR) : LIB_DIR;
export const BASE_TIERS = ['common', 'platform'];
export const CONFIG_FILE = 'bdd.config.json';

export interface BindingSource {
  dir: string;
  /** What a generated spec imports for a module of this source: a package subpath for the
   * library's bindings, the file itself (made relative by the compiler) for a project's own. */
  specifierFor(file: string): string;
}

export interface ProjectConfig {
  /** Library tiers beyond the base ones, by directory name under `bindings/tiers/`. */
  tiers?: string[];
}

export interface Project {
  /** The directory holding `features/`, `bindings/` and `generated/`. */
  root: string;
  /** The package (or library) the project belongs to — where `node_modules` is. */
  packageDir: string;
  name: string;
  tiers: string[];
  sources: BindingSource[];
  featuresDir: string;
  generatedDir: string;
}

function posix(p: string): string {
  return p.split(sep).join('/');
}

export function findRoot(cwd: string): string {
  for (const candidate of [join(cwd, 'bdd'), cwd]) {
    if (existsSync(join(candidate, 'features')))
      return candidate;
  }
  throw new Error(`no bdd project at ${cwd}: expected ${join(cwd, 'bdd', 'features')} or ${join(cwd, 'features')}`);
}

export function readConfig(root: string): ProjectConfig {
  const file = join(root, CONFIG_FILE);
  return existsSync(file) ? JSON.parse(readFileSync(file, 'utf8')) as ProjectConfig : {};
}

export function availableTiers(): string[] {
  const dir = join(LIB_DIR, 'bindings', 'tiers');
  return existsSync(dir) ? readdirSync(dir).filter((n) => statSync(join(dir, n)).isDirectory()).sort() : [];
}

export function librarySource(tier: string): BindingSource {
  const bindings = join(LIB_DIR, 'bindings');
  return {
    dir: join(bindings, tier),
    specifierFor: (file) => `${PACKAGE_NAME}/bindings/${posix(relative(bindings, file)).replace(/\.[cm]?[jt]s$/, '')}`,
  };
}

function packageName(dir: string): string {
  const file = join(dir, 'package.json');
  if (!existsSync(file))
    return basename(dir);
  return (JSON.parse(readFileSync(file, 'utf8')) as {name?: string}).name ?? basename(dir);
}

export function loadProject(cwd: string): Project {
  const root = findRoot(cwd);
  const isLibrary = resolve(root) === resolve(LIB_ROOT);
  const packageDir = !isLibrary && basename(root) === 'bdd' ? dirname(root) : root;
  const tiers = readConfig(root).tiers ?? [];
  for (const tier of tiers) {
    if (BASE_TIERS.includes(tier))
      throw new Error(`${CONFIG_FILE}: tier "${tier}" is always loaded — remove it`);
    if (!existsSync(join(LIB_DIR, 'bindings', 'tiers', tier)))
      throw new Error(`${CONFIG_FILE}: unknown tier "${tier}" (available: ${availableTiers().join(', ') || 'none'})`);
  }
  const sources = [...BASE_TIERS.map(librarySource), ...tiers.map((t) => librarySource(join('tiers', t)))];
  const local = join(root, 'bindings');
  if (!isLibrary && existsSync(local))
    sources.push({dir: local, specifierFor: (file) => file});
  return {root, packageDir, name: packageName(packageDir), tiers, sources,
    featuresDir: join(root, 'features'), generatedDir: join(root, 'generated')};
}
