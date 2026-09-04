/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';
import * as yaml from 'js-yaml';
import {NodeDapi} from '../utils/node-dapi';
import {createClient} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {BytesKind, DEFAULT_TYPES, nqNameOf, resolveTypes} from '../utils/migrate/registry';
import * as bundle from '../utils/migrate/bundle';
import {Selection, expand, normalizeSince, pullBytes, select} from '../utils/migrate/walker';
import {ConflictPolicy, Row, plan, push, summarize} from '../utils/migrate/pusher';

const SELECTION_USAGE = '  [<nqName|id>...] [--type t,t] [--namespace ns] [--space s] [--name n] [--author login]\n' +
  '  [--tag t] [--since 2w] [--filter expr] [--no-deps] [--no-include-data] [--include-files]';
const PULL_USAGE = `Usage: grok s pull --out <dir> [--replace] [--host <alias>]\n${SELECTION_USAGE}`;
const MIGRATE_USAGE = 'Usage: grok s migrate --from <alias> --to <alias> [--dry-run] [--keep]\n' +
  '  [--on-conflict fail|skip|duplicate|adopt] [--creds <file.yaml>]\n' + SELECTION_USAGE;

export async function handleMigrate(dapi: NodeDapi, verb: string, rest: string[], argv: any,
                                    output: OutputFormat): Promise<boolean> {
  if (verb === 'pull') return await handlePull(dapi, rest, argv, output);
  if (verb === 'push') return await handlePush(dapi, rest, argv, output);
  if (verb === 'migrate') return await handleTransfer(rest, argv, output);
  if (verb === 'diff') return await handleDiff(dapi, rest, argv, output);
  if (verb === 'bundle') return handleBundle(rest, output);
  printError(new Error(`Unknown migrate verb '${verb}'. Valid: pull, push, migrate, diff, bundle ls`));
  return false;
}

const hasSelection = (rest: string[], argv: any): boolean =>
  !!rest.length || !!argv.type || ['name', 'namespace', 'space', 'author', 'tag', 'since', 'filter', 'f'].some((f) => argv[f]);

async function handlePull(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat,
                          print: boolean = true, report?: {dropped: number}): Promise<boolean> {
  const out: string = argv.out ?? '';
  if (!out || !hasSelection(rest, argv)) {
    printError(new Error(out ? `Refusing to pull the whole server — pass entity names, --type, or a filter flag.\n${PULL_USAGE}` : PULL_USAGE));
    return false;
  }

  const sel: Selection = {
    ...resolveTypes(argv.type ? String(argv.type).split(',') : DEFAULT_TYPES),
    names: rest,
    name: argv.name,
    namespace: argv.namespace,
    space: argv.space,
    author: argv.author,
    tag: argv.tag,
    since: normalizeSince(argv.since),
    filter: argv.filter ?? argv.f,
  };

  const notes: Row[] = [];
  const note = (row: Row) => notes.push(row);
  const selected = await select(dapi, sel, note);
  const entities = argv.deps !== false ? await expand(dapi, selected, note) : selected;
  const kinds: BytesKind[] = [];
  if (argv['include-data'] !== false) kinds.push('tables');
  if (argv['include-files']) kinds.push('files');
  const bytes = await pullBytes(dapi, entities, note, kinds);

  const info = await dapi.serverInfo();
  const user = await dapi.client.get('/users/current');
  bundle.write(out, entities, {
    source: {
      url: dapi.client.baseUrl,
      version: info.version,
      commit: info.commit,
      userNamespace: user?.project?.name ? `${user.project.name}:` : '',
    },
    args: process.argv.slice(3),
    packages: notes.filter((n) => n.reason === 'package_entity').map((n) => n.detail!).filter(Boolean),
  }, {replace: !!argv.replace}, bytes);

  const rows: Row[] = [...notes];
  for (const [, {type, json}] of entities)
    rows.push({name: nqNameOf(json), entityType: type, action: 'info', reason: 'pulled'});
  // An entity the server would not hand over — or whose bytes it would not — is missing from
  // the bundle, and pushing it would quietly promote less than was asked for.
  const dropped = notes.filter((n) => n.reason === 'fetch_failed' || n.reason === 'no_data').length;
  if (dropped) {
    rows.push({name: out, entityType: 'Bundle', action: 'failed', reason: 'incomplete',
      detail: `${dropped} entities could not be read in full`});
    process.exitCode = 1;
    if (report)
      report.dropped = dropped;
  }
  // `migrate` silences the pull's own report, but a refusal has to say why.
  if (print || dropped)
    printOutput(rows, output);
  return true;
}

const ENV_RE = /\$\{(\w*)\}/g;

/**
 * Target-side secrets, keyed by connection nqName. `${VAR}` is resolved from the
 * environment, as `grok publish` does for `connections/*.json` — after the YAML is parsed,
 * so a broken file is reported without a secret in the message.
 */
export function loadCreds(file?: string): Record<string, any> | undefined {
  if (!file) return undefined;
  const loaded = (yaml.load(fs.readFileSync(file, 'utf8')) ?? {}) as Record<string, any>;
  const missing: string[] = [];
  const creds: Record<string, any> = {};
  for (const [key, params] of Object.entries(loaded)) {
    if (!params || typeof params !== 'object' || Array.isArray(params))
      throw new Error(`${file}: "${key}" must be a map of connection parameters, e.g. '${key}: {password: \${VAR}}'`);
    creds[key] = {};
    for (const [name, value] of Object.entries(params as Record<string, any>))
      creds[key][name] = typeof value !== 'string' ? value : value.replace(ENV_RE, (whole, env) => {
        const resolved = process.env[env];
        if (resolved !== undefined) return resolved;
        missing.push(env);
        return whole;
      });
  }
  if (missing.length)
    throw new Error(`${file}: cannot find environment variable "${[...new Set(missing)].join('", "')}"`);
  return creds;
}

const POLICIES: ConflictPolicy[] = ['fail', 'skip', 'duplicate', 'adopt'];

function conflictPolicy(argv: any): ConflictPolicy {
  const policy: ConflictPolicy = argv['on-conflict'] ?? 'fail';
  if (!POLICIES.includes(policy))
    throw new Error(`Unsupported conflict policy '${policy}'. Valid: ${POLICIES.join(', ')}`);
  return policy;
}

async function handlePush(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const dir = rest[0];
  if (!dir) {
    printError(new Error('Usage: grok s push <bundle-dir> [--dry-run] [--on-conflict fail|skip|duplicate|adopt] [--creds <file.yaml>] [--host <alias>]'));
    return false;
  }
  const onConflict = conflictPolicy(argv);
  const creds = loadCreds(argv.creds);
  const dryRun = !!argv['dry-run'];
  const result = await push(dapi, bundle.read(dir), {dryRun, onConflict, creds}, (rows) => {
    if (output === 'table' && !dryRun) {
      console.log('Plan:');
      printOutput(rows, output);
      console.log('');
    }
  });
  if (result.items.some((r) => r.action === 'failed'))
    process.exitCode = 1;
  if (output === 'json') {
    printOutput({...result, items: result.items.map((r) => ({...r, detail: r.detail ?? ''}))}, 'json');
    return true;
  }
  if (!dryRun) console.log('Result:');
  printOutput(result.items, output);
  return true;
}

/**
 * A bundle pulled from one instance and pushed into another, with a temporary bundle
 * directory in between — the same two verbs, so nothing behaves differently.
 */
async function handleTransfer(rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  if (!argv.from || !argv.to || !hasSelection(rest, argv)) {
    printError(new Error(MIGRATE_USAGE));
    return false;
  }
  const from = new NodeDapi(await createClient(String(argv.from), !!argv.admin));
  const to = new NodeDapi(await createClient(String(argv.to), !!argv.admin));
  if (from.client.baseUrl === to.client.baseUrl)
    throw new Error(`--from and --to are the same server (${to.client.baseUrl}) — nothing to migrate`);

  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'grok-migrate-'));
  try {
    const report = {dropped: 0};
    if (!await handlePull(from, rest, {...argv, out: dir}, output, output !== 'json', report))
      return false;
    if (report.dropped) {
      printError(new Error(`Refusing to push a partial bundle: ${report.dropped} entities could not be read. ` +
        'Re-run, or pull with --keep and push the bundle yourself.'));
      return true;
    }
    return await handlePush(to, [dir], argv, output);
  } finally {
    // stderr: `--output json` must stay one parseable document.
    if (argv.keep)
      console.error(`Bundle kept at ${dir}`);
    else
      fs.rmSync(dir, {recursive: true, force: true});
  }
}

/** What a push would do, read-only: the plan plus the top-level keys that differ. */
async function handleDiff(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const dir = rest[0];
  if (!dir) {
    printError(new Error('Usage: grok s diff <bundle-dir> [--host <alias>]'));
    return false;
  }
  const read = bundle.read(dir);
  const {rows} = await plan(dapi, read, {onConflict: argv['on-conflict'] ? conflictPolicy(argv) : 'skip', idmap: {...read.idmap}});
  const items: Row[] = rows.map((r) => ({...r, detail: r.detail ?? ''}));
  printOutput(output === 'json' ? summarize(items, dapi, 'dry-run') : items, output);
  return true;
}

function handleBundle(rest: string[], output: OutputFormat): boolean {
  if (rest[0] !== 'ls' || !rest[1]) {
    printError(new Error('Usage: grok s bundle ls <bundle-dir>'));
    return false;
  }
  printOutput(bundle.list(rest[1]), output);
  return true;
}
