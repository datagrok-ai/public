import fs from 'fs';
import path from 'path';
import {findDomainManifests, loadDomainManifest} from './api';
import {createClient} from '../utils/server-client';
import * as color from '../utils/color-utils';
import {buildSnapshot, describeChange, diff, emptySnapshot, hashOf, scaffold, serialize,
  Change, Snapshot} from '../utils/domain-snapshot';

const VERBS = ['seal', 'check', 'diff', 'migrate'];

interface SchemaArgs {
  _: string[];
  help?: boolean;
  dir?: string;
  name?: string;
  /** `--server` alone = the default host; `--server <alias|url>` = that one. */
  server?: string | boolean;
}

function short(hash: string | undefined): string {
  return hash == null ? '(none)' : hash.substring(0, 12);
}

function printChanges(changes: Change[]): void {
  for (const c of changes)
    console.log(`  * ${describeChange(c)}`);
}

function readSealed(sealedPath: string): Snapshot | null {
  return fs.existsSync(sealedPath) ? JSON.parse(fs.readFileSync(sealedPath, 'utf8')) : null;
}

/** The snapshot recorded for [schema] on the server, or null when none is (404). */
async function fetchRecorded(host: string | boolean, schema: string): Promise<Snapshot | null> {
  const client = await createClient(typeof host === 'string' ? host : undefined);
  let body: any;
  try {
    body = await client.get(`/domains/schemas/${encodeURIComponent(schema)}/snapshot`);
  } catch (x: any) {
    if (x.apiError?.errorCode === 404)
      return null;
    throw x;
  }
  if (body?.['#type'] === 'ApiError')
    throw new Error(body.message ?? body.error);
  if (typeof body !== 'object' || body?.tables == null)
    throw new Error(`Unexpected response for the recorded snapshot of '${schema}'`);
  return body;
}

function seal(sealedPath: string, relPath: string, snapshot: Snapshot): void {
  fs.writeFileSync(sealedPath, serialize(snapshot), 'utf8');
  console.log(`Sealed ${relPath} (${short(snapshot.hash)})`);
}

/** The next migration number: one past the highest `NNNN_` prefix in [dir], from 0001. */
function nextMigrationNumber(dir: string): number {
  let max = 0;
  if (fs.existsSync(dir)) {
    for (const f of fs.readdirSync(dir)) {
      const m = /^(\d+)_.*\.sql$/.exec(f);
      if (m != null)
        max = Math.max(max, Number(m[1]));
    }
  }
  return max + 1;
}

async function run(argv: SchemaArgs, verb: string, packageDir: string): Promise<boolean> {
  const manifestPaths = argv.dir
    ? [path.join(packageDir, String(argv.dir), 'schema.json')] : findDomainManifests(packageDir);
  if (manifestPaths.length === 0 || !fs.existsSync(manifestPaths[0])) {
    color.error(argv.dir ? `No schema.json in ${argv.dir}` :
      'No databases/<schema>/schema.json manifests found. Run the command from the package directory');
    return false;
  }
  let ok = true;
  for (const manifestPath of manifestPaths) {
    const relPath = path.relative(packageDir, manifestPath);
    const manifest = loadDomainManifest(manifestPath, relPath);
    if (manifest == null)
      return false;
    let current: Snapshot;
    try {
      current = buildSnapshot(manifest);
    } catch (x: any) {
      color.error(`${relPath}: ${x.message}`);
      return false;
    }
    const dir = path.dirname(manifestPath);
    const sealedPath = path.join(dir, 'snapshot.json');
    const sealedRel = path.relative(packageDir, sealedPath);
    if (verb === 'seal') {
      seal(sealedPath, sealedRel, current);
      continue;
    }

    const baseline = argv.server ? await fetchRecorded(argv.server, manifest.name) : readSealed(sealedPath);
    const source = argv.server ? 'recorded on the server' : `sealed in ${sealedRel}`;
    const baseHash = baseline == null ? undefined : hashOf(baseline);
    const changes = diff(baseline ?? emptySnapshot(manifest.name), current);
    const span = `${short(baseHash)} -> ${short(current.hash)}`;

    if (verb === 'check') {
      if (baseline == null) {
        color.error(`${manifest.name}: no snapshot ${source}; run 'grok schema seal'`);
        ok = false;
      } else if (baseHash === current.hash)
        console.log(`${manifest.name}: up to date (${short(current.hash)})`);
      else {
        color.error(`${manifest.name}: the manifest differs from the snapshot ${source} (${span}):`);
        printChanges(changes);
        ok = false;
      }
      continue;
    }

    if (verb === 'diff') {
      if (baseline == null)
        color.warn(`${manifest.name}: no snapshot ${source}; every table is new`);
      if (changes.length === 0)
        console.log(`${manifest.name}: no changes (${short(current.hash)})`);
      else {
        console.log(`${manifest.name}: ${changes.length} change(s), ${span}:`);
        printChanges(changes);
      }
      continue;
    }

    // migrate
    if (changes.length === 0) {
      console.log(`${manifest.name}: no changes (${short(current.hash)}), nothing to migrate`);
      continue;
    }
    const physical = changes.filter((c) => c.physical);
    if (physical.length === 0) {
      console.log(`${manifest.name}: nothing physical to migrate — registry-only changes, ` +
        'converged by the next deploy:');
      printChanges(changes);
      seal(sealedPath, sealedRel, current);
      continue;
    }
    const migrationsDir = path.join(dir, 'migrations');
    const downDir = path.join(migrationsDir, 'down');
    const file = `${String(nextMigrationNumber(migrationsDir)).padStart(4, '0')}_${argv.name}.sql`;
    const scripts = scaffold(changes, manifest.name, {from: baseHash, to: current.hash,
      title: `${manifest.name}: ${argv.name}`, generatedBy: 'grok schema migrate'});
    fs.mkdirSync(downDir, {recursive: true});
    fs.writeFileSync(path.join(migrationsDir, file), scripts.up, 'utf8');
    fs.writeFileSync(path.join(downDir, file), scripts.down, 'utf8');
    console.log(`${manifest.name}: ${physical.length} physical change(s), ${span}:`);
    printChanges(changes);
    console.log(`Wrote ${path.relative(packageDir, path.join(migrationsDir, file))} and ` +
      `${path.relative(packageDir, path.join(downDir, file))} — review before committing`);
    seal(sealedPath, sealedRel, current);
  }
  if (!ok)
    process.exitCode = 1;
  return true;
}

export async function schema(argv: SchemaArgs, packageDir: string = process.cwd()): Promise<boolean> {
  const verb: string | undefined = argv._[1];
  if (argv.help || verb == null || !VERBS.includes(verb))
    return false;
  if (verb === 'migrate' && !/^[A-Za-z0-9_-]+$/.test(String(argv.name ?? ''))) {
    color.error('migrate needs --name <x> (letters, digits, _ and -)');
    return false;
  }
  try {
    return await run(argv, verb, packageDir);
  } catch (x: any) {
    // A runtime failure is not a usage error: report it and exit non-zero without the help block.
    color.error(x.message ?? String(x));
    process.exitCode = 1;
    return true;
  }
}
