import {cp, mkdir, readFile, readdir, rm, writeFile} from 'node:fs/promises';
import {execFileSync} from 'node:child_process';
import {resolve} from 'node:path';
import {fileURLToPath} from 'node:url';

const source = process.argv[2];
if (!source)
  throw new Error('Usage: node scripts/update-crux.mjs /path/to/crux-js (run npm run build there first)');
const repo = resolve(source);
const target = fileURLToPath(new URL('../wasm/crux/', import.meta.url));
const revision = (path) => execFileSync('git', ['-C', path, 'rev-parse', 'HEAD'], {encoding: 'utf8'}).trim();
const manifest = JSON.parse(await readFile(resolve(repo, 'package.json'), 'utf8'));
const provenance = {
  name: manifest.name,
  version: manifest.version,
  license: manifest.license,
  cruxJs: revision(repo),
  cruxCore: revision(resolve(repo, '../crux-core')),
};
// Preserve the directory structure used by the ESM worker and wasm-bindgen loader.
await mkdir(target, {recursive: true});
for (const file of await readdir(target)) {
  if (file !== 'README.md')
    await rm(resolve(target, file), {recursive: true, force: true});
}
await cp(resolve(repo, 'dist'), target, {recursive: true});
await writeFile(resolve(target, 'provenance.json'), JSON.stringify(provenance, null, 2) + '\n');
console.log(`Copied ${manifest.name} ${manifest.version} (${provenance.cruxJs}) to ${target}`);
