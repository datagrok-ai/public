/// Packages and libraries from their package.json (build-plan.md WO-3a): `package` and `library` nodes,
/// `depends-on` between them, and the semantic types a package.json declares.
import * as fs from 'fs';
import * as path from 'path';
import {Emitter} from '../../emitter';
import {Row} from '../../normalize';
import {BuildContext, Extractor} from '../../registry';
import {pkgId, libId, semtypeId, posix} from '../../ids';

export interface PackageFolder {
  /** The folder name, which is the platform package name (`Chem`, `TensorFlow.js`). */
  folder: string;
  /** Posix path from the monorepo root. */
  dir: string;
  json: Record<string, any>;
}

const DEPENDENCY_KINDS: [string, string][] = [['dependencies', 'runtime'], ['devDependencies', 'dev'], ['peerDependencies', 'peer'], ['optionalDependencies', 'optional']];

export const packagesExtractor: Extractor = {
  name: 'ts-packages',
  layer: 'public',
  modes: ['full', 'public'],
  run(ctx: BuildContext, emitter: Emitter): void {
    const packages = listPackages(ctx.repoRoot);
    const libraries = listLibraries(ctx.repoRoot);
    const npm = new Map<string, string>();
    for (const p of packages) if (typeof p.json.name === 'string') npm.set(p.json.name, pkgId(p.folder));
    for (const l of libraries) if (typeof l.json.name === 'string') npm.set(l.json.name, libId(l.folder));
    npm.set('datagrok-api', libId('js-api'));
    for (const p of packages) {
      const id = pkgId(p.folder);
      const {json} = p;
      emitter.node({type: 'package', id, name: p.folder, description: json.description, friendly_name: json.friendlyName, version: json.version,
        category: json.category, author: json.author?.name, service: json.servicePackage === true ? true : undefined,
        settings: Array.isArray(json.properties) ? json.properties.map((s: any) => s?.name).filter((n: unknown) => typeof n === 'string') : undefined,
        sources: json.sources, npm: json.name, language: ['src/package.ts', 'src/package.g.ts', 'tsconfig.json'].some((f) => fs.existsSync(path.join(ctx.repoRoot, p.dir, f))) ? 'ts' : 'js',
        path: p.dir, provenance: 'registry', source_layer: 'public'});
      for (const st of json.meta?.semanticTypes ?? json.semanticTypes ?? []) {
        if (typeof st?.semType !== 'string') continue;
        emitter.node({type: 'semantic-type', id: semtypeId(st.semType), name: st.semType, description: st.description, language: 'other', declared_in: id, provenance: 'registry', source_layer: 'public'});
        emitter.edge({type: 'declares', from: id, to: semtypeId(st.semType), derived_by: 'registry', confidence: 1, evidence: [`${p.dir}/package.json`]});
      }
      emitDependencies(emitter, id, p, npm);
    }
    for (const l of libraries) {
      const id = libId(l.folder);
      emitter.node({type: 'library', id, name: l.folder, description: l.json.description, npm: l.json.name, version: l.json.version, language: 'ts',
        path: l.dir, provenance: 'registry', source_layer: 'public'});
      emitDependencies(emitter, id, l, npm);
    }
  },
};

/** Every `package.json` one level under `public/packages`, sorted by folder. */
export function listPackages(repoRoot: string): PackageFolder[] {
  return readFolders(repoRoot, 'public/packages');
}

/** Every `package.json` one level under `public/libraries`, and the JS API itself as `js-api`. */
export function listLibraries(repoRoot: string): PackageFolder[] {
  const libraries = readFolders(repoRoot, 'public/libraries');
  const jsApi = readJson(path.join(repoRoot, 'public', 'js-api', 'package.json'));
  return jsApi ? [...libraries, {folder: 'js-api', dir: 'public/js-api', json: jsApi}] : libraries;
}

function readFolders(repoRoot: string, dir: string): PackageFolder[] {
  const root = path.join(repoRoot, dir);
  if (!fs.existsSync(root)) return [];
  const out: PackageFolder[] = [];
  for (const folder of fs.readdirSync(root).sort()) {
    const json = readJson(path.join(root, folder, 'package.json'));
    if (json) out.push({folder, dir: posix(`${dir}/${folder}`), json});
  }
  return out;
}

function readJson(file: string): Record<string, any> | null {
  if (!fs.existsSync(file)) return null;
  try {
    return JSON.parse(fs.readFileSync(file, 'utf8'));
  } catch {
    return null;
  }
}

/** `@datagrok/*`, `@datagrok-libraries/*` and `datagrok-api` dependencies as depends-on; other npm packages are not in the graph. */
function emitDependencies(emitter: Emitter, from: string, p: PackageFolder, npm: Map<string, string>): void {
  for (const [key, kind] of DEPENDENCY_KINDS) {
    const deps = p.json[key];
    if (!deps || typeof deps !== 'object') continue;
    for (const [name, range] of Object.entries(deps as Record<string, unknown>)) {
      if (!/^(@datagrok\/|@datagrok-libraries\/|datagrok-api$)/.test(name)) continue;
      const to = npm.get(name);
      if (!to) {
        emitter.problem('unresolved_ids', `${p.dir}/package.json: ${key} ${name} is not under public/packages or public/libraries`);
        continue;
      }
      const row: Row = {type: 'depends-on', from, to, kind, derived_by: 'registry', confidence: 1, evidence: [`${p.dir}/package.json`]};
      if (typeof range === 'string') row.range = range;
      emitter.edge(row);
    }
  }
}
