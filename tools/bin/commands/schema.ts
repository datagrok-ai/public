import fs from 'fs';
import path from 'path';
import {findDomainManifests, loadDomainManifest} from './api';
import {createClient} from '../utils/server-client';
import * as color from '../utils/color-utils';
import {buildSnapshot, chainGap, describeChange, diff, emptySnapshot, hashOf, parseMigrationHeader, scaffold,
  serialize, squashScripts, Change, MigrationScript, Snapshot} from '../utils/domain-snapshot';

const VERBS = ['seal', 'check', 'diff', 'migrate', 'squash'];
const scriptName = /^(\d+)_.*\.sql$/;

interface SchemaArgs {
  _: string[];
  help?: boolean;
  dir?: string;
  name?: string;
  /** `--server` alone = the default host; `--server <alias|url>` = that one. */
  server?: string | boolean;
  all?: boolean;
  from?: number | string;
  to?: number | string;
  discard?: boolean;
  yes?: boolean;
  'force-missing-down'?: boolean;
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
  for (const f of fs.readdirSync(dir)) {
    const m = scriptName.exec(f);
    if (m != null)
      max = Math.max(max, Number(m[1]));
  }
  return max + 1;
}

/** The schema dir's `NNNN_*.sql` scripts carrying a transition header, in name order,
 * each with its `down/` twin when there is one. */
function readMigrationScripts(dir: string): (MigrationScript & {number: number})[] {
  const scripts: (MigrationScript & {number: number})[] = [];
  for (const file of fs.readdirSync(dir).sort()) {
    const m = scriptName.exec(file);
    if (m == null)
      continue;
    const up = fs.readFileSync(path.join(dir, file), 'utf8');
    const header = parseMigrationHeader(up);
    if (header == null)
      continue;
    const downPath = path.join(dir, 'down', file);
    scripts.push({file, number: Number(m[1]), from: header.from, to: header.to, up,
      down: fs.existsSync(downPath) ? fs.readFileSync(downPath, 'utf8') : undefined});
  }
  return scripts;
}

/** `squash`: folds a continuous chain of scripts into one (or, with --discard, deletes them
 * all). False = a validation failure, reported. */
function squash(dir: string, relDir: string, schema: string, argv: SchemaArgs): boolean {
  const scripts = readMigrationScripts(dir);
  const remove = (s: MigrationScript) => {
    fs.rmSync(path.join(dir, s.file));
    if (s.down != null)
      fs.rmSync(path.join(dir, 'down', s.file));
  };
  if (argv.discard) {
    if (!argv.yes) {
      color.error(`${schema}: --discard deletes every migration script in ${relDir}; pass --yes to confirm`);
      return false;
    }
    for (const s of scripts)
      remove(s);
    color.warn(`${schema}: deleted ${scripts.length} migration script(s) and their down twins — the manifest is ` +
      'the baseline; stands behind the current declaration can no longer migrate by script');
    return true;
  }

  const from = argv.from == null ? -Infinity : Number(argv.from);
  const to = argv.to == null ? Infinity : Number(argv.to);
  const selected = argv.all ? scripts : scripts.filter((s) => s.number >= from && s.number <= to);
  if (selected.length === 0 && argv.all) {
    console.log(`${schema}: no migration scripts in ${relDir}, nothing to squash`);
    return true;
  }
  if (selected.length < 2) {
    color.error(`${schema}: ${selected.length} script(s) selected in ${relDir}; squash needs at least two`);
    return false;
  }
  const gap = chainGap(selected);
  if (gap != null) {
    color.error(`${schema}: the scripts do not chain — ${gap}`);
    return false;
  }
  const missing = selected.filter((s) => s.down == null).map((s) => s.file);
  if (missing.length > 0 && !argv['force-missing-down']) {
    color.error(`${schema}: no down script for ${missing.join(', ')}; write them, or pass --force-missing-down ` +
      'to leave a TODO in the squashed down script');
    return false;
  }

  // Keeps the first script's number, so the result stays in place among unrelated scripts.
  const file = `${scriptName.exec(selected[0].file)![1]}_${argv.name}.sql`;
  const out = squashScripts(selected, schema, argv.name!);
  for (const s of selected)
    remove(s);
  fs.mkdirSync(path.join(dir, 'down'), {recursive: true});
  fs.writeFileSync(path.join(dir, file), out.up, 'utf8');
  fs.writeFileSync(path.join(dir, 'down', file), out.down, 'utf8');
  console.log(`${schema}: squashed ${selected.map((s) => s.file).join(', ')} into ${file}`);
  console.log(`Wrote ${path.join(relDir, file)} and ${path.join(relDir, 'down', file)}; deleted ` +
    `${selected.length} script(s) and their down twins — review before committing`);
  return true;
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
    if (verb === 'squash') {
      if (!squash(dir, path.relative(packageDir, dir), manifest.name, argv))
        ok = false;
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
    if (physical.every((c) => c.auto))
      console.log('No change needs a migration script — `grok schema seal` is enough; scripts written for rollback.');
    // The up script sits next to schema.json, where the deployer runs NNNN_*.sql files
    // non-recursively; the down script goes in a subdirectory it skips.
    const downDir = path.join(dir, 'down');
    const file = `${String(nextMigrationNumber(dir)).padStart(4, '0')}_${argv.name}.sql`;
    const scripts = scaffold(changes, manifest.name, {from: baseHash, to: current.hash,
      title: `${manifest.name}: ${argv.name}`, generatedBy: 'grok schema migrate'});
    fs.mkdirSync(downDir, {recursive: true});
    fs.writeFileSync(path.join(dir, file), scripts.up, 'utf8');
    fs.writeFileSync(path.join(downDir, file), scripts.down, 'utf8');
    console.log(`${manifest.name}: ${physical.length} physical change(s), ${span}:`);
    printChanges(changes);
    console.log(`Wrote ${path.relative(packageDir, path.join(dir, file))} and ` +
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
  if (verb === 'squash') {
    if (argv.discard && !argv.all) {
      color.error('--discard applies to --all only');
      return false;
    }
    if (!argv.all && argv.from == null && argv.to == null) {
      color.error('squash needs --all, or --from NNNN and/or --to NNNN');
      return false;
    }
  }
  const needsName = verb === 'migrate' || (verb === 'squash' && !argv.discard);
  if (needsName && !/^[A-Za-z0-9_-]+$/.test(String(argv.name ?? ''))) {
    color.error(`${verb} needs --name <x> (letters, digits, _ and -)`);
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
