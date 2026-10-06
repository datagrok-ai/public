import fs from 'fs';
import path from 'path';

/**
 * The shared toolchain (`@datagrok/build-config`: rspack factory, TypeScript) as seen from `dir`:
 * the workspace root's copy inside public/, the package's own devDependency in a standalone package.
 */
export function toolchain(dir: string): {buildConfig: any; buildConfigDir: string; tsc: string} {
  let pkgJson: string;
  try {
    pkgJson = require.resolve('@datagrok/build-config/package.json', {paths: [dir]});
  }
  catch {
    throw new Error('@datagrok/build-config is not installed: run `pnpm install` at the root of public/, ' +
      'or add it as a devDependency of a standalone package');
  }
  const buildConfigDir = path.dirname(pkgJson);
  const tsPkg = require.resolve('typescript/package.json', {paths: [buildConfigDir]});
  return {
    buildConfig: require(buildConfigDir),
    buildConfigDir,
    tsc: path.join(path.dirname(tsPkg), 'bin', 'tsc'),
  };
}

const ASSET = /\.(css|scss|wasm|js|mjs|cjs|json|png|jpe?g|gif|svg|txt|csv|md|html|sdf|mol)$/i;

/**
 * A library's dist/ must be self-contained: the css, wasm, json and plain-js files its sources
 * import relative to themselves are mirrored next to the compiled output. A .js next to a .ts of
 * the same name is a stale in-place compile, not an asset.
 */
export function copyLibraryAssets(dir: string): number {
  if (!fs.existsSync(path.join(dir, 'dist')))
    return 0;
  let n = 0;
  const walk = (rel: string) => {
    for (const e of fs.readdirSync(path.join(dir, rel), {withFileTypes: true})) {
      if (e.name === 'node_modules' || e.name === 'dist' || e.name.startsWith('.'))
        continue;
      const r = rel ? `${rel}/${e.name}` : e.name;
      if (e.isDirectory())
        walk(r);
      else if (rel && ASSET.test(e.name) && e.name !== 'package.json' && !isCompiledSibling(dir, r)) {
        const to = path.join(dir, 'dist', r);
        fs.mkdirSync(path.dirname(to), {recursive: true});
        const from = path.join(dir, r);
        if (!fs.existsSync(to) || fs.statSync(to).mtimeMs < fs.statSync(from).mtimeMs) {
          fs.copyFileSync(from, to);
          n++;
        }
      }
    }
  };
  walk('');
  if (n)
    console.log(`copied ${n} asset file(s) into dist/`);
  return n;
}

/**
 * Keeps a library's package.json `exports` in step with its source tree. A library that ships dist/
 * has an `exports` map ending in the `./*` wildcard; the wildcard cannot resolve a directory import
 * (`@datagrok-libraries/bio/src/trees`), so every directory with an index.ts needs an explicit entry.
 * Adds the missing ones and leaves everything else (hand-written entries, the wildcards, entries for
 * directories that no longer exist) untouched. Packages without the wildcard are not managed.
 */
export function ensureLibraryExports(dir: string): boolean {
  const pj = path.join(dir, 'package.json');
  const p = JSON.parse(fs.readFileSync(pj, 'utf8'));
  const ex = p.exports;
  if (!ex || typeof ex !== 'object' || !ex['./*'])
    return false;
  const indexDirs = new Set<string>();
  const walk = (d: string, rel: string) => {
    for (const e of fs.readdirSync(d, {withFileTypes: true})) {
      if (e.name === 'node_modules' || e.name === 'dist' || e.name.startsWith('.'))
        continue;
      const r = rel ? `${rel}/${e.name}` : e.name;
      if (e.isDirectory())
        walk(path.join(d, e.name), r);
      else if (e.name === 'index.ts' && rel)
        indexDirs.add(rel);
    }
  };
  walk(dir, '');
  // Add-only: hand-written entries (aliases such as js-api's "./u2core") are never touched.
  const next: Record<string, any> = {...ex};
  let changed = false;
  for (const rel of [...indexDirs].sort()) {
    const key = `./${rel}`;
    if (!(key in next)) {
      next[key] = {types: `./dist/${rel}/index.d.ts`, default: `./dist/${rel}/index.js`};
      changed = true;
    }
  }
  if (!changed)
    return false;
  // wildcards last: exports are matched by specificity, but keeping them last reads better
  const wild = ['./*.js', './*'].filter((k) => k in next);
  const ordered: Record<string, any> = {};
  for (const [k, v] of Object.entries(next)) if (!wild.includes(k)) ordered[k] = v;
  for (const k of wild) ordered[k] = next[k];
  p.exports = ordered;
  fs.writeFileSync(pj, JSON.stringify(p, null, 2) + '\n');
  console.log(`package.json exports updated for ${path.basename(dir)}`);
  return true;
}

function isCompiledSibling(dir: string, rel: string): boolean {
  const m = rel.match(/^(.*)\.(js|mjs|cjs)$/);
  return !!m && ['.ts', '.tsx'].some((e) => fs.existsSync(path.join(dir, m[1] + e)));
}
