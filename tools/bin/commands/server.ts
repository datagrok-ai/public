/// Docs: [Grok Dapi](/docs/plans/grok-dapi/)
import * as fs from 'fs';
import * as path from 'path';
import {NodeDapi, BatchRequest, BatchOperation, InternalDataSource, PackageOpResult} from '../utils/node-dapi';
import {createClient} from '../utils/server-client';
import {printOutput, printBatchOutput, printError, setOutputFormat, OutputFormat} from '../utils/server-output';
import {handleMigrate} from './server-migrate';
import {isUuid} from '../utils/migrate/registry';
import {resolveEntity} from '../utils/migrate/walker';

/** `queries|scripts get <nqName>` resolves through `/entities`, which needs the entity type. */
const ENTITY_TYPES: Record<string, string> = {queries: 'DataQuery', scripts: 'Script', reports: 'UserReport'};

const ENTITIES = ['users', 'groups', 'functions', 'connections', 'queries', 'scripts', 'packages', 'reports', 'files', 'tables'];
const VERBS = ['list', 'get', 'delete'];

export async function server(argv: any): Promise<boolean> {
  const args: string[] = argv['_'].slice(1);
  const entity: string | undefined = args[0];
  const verb: string | undefined = args[1];
  const rest: string[] = args.slice(2);

  const output: OutputFormat = argv.output ?? argv.o ?? 'table';
  setOutputFormat(output);
  const limit: number = Number(argv.limit ?? argv.l ?? 50);
  const offset: number = Number(argv.offset ?? 0);
  const filter: string = argv.filter ?? argv.f ?? '';
  const host: string | undefined = argv.host;
  const recursive: boolean = !!(argv.r ?? argv.recursive);

  if (!entity || argv.help) {
    console.log(HELP_SERVER);
    return true;
  }

  let client;
  try {
    client = await createClient(host);
  } catch (err: any) {
    printError(err);
    return false;
  }
  const dapi = new NodeDapi(client);

  try {
    // await each handler so a rejection lands in the catch below instead of
    // escaping to grok.js as an unhandled rejection (stack trace + exit 255)
    if (['pull', 'push', 'migrate', 'diff', 'bundle'].includes(entity))
      return await handleMigrate(dapi, entity, [verb, ...rest].filter(Boolean), argv, output);
    if (entity === 'batch') return await handleBatch(dapi, argv, verb, rest, output);
    if (entity === 'raw') return await handleRaw(dapi, verb, rest, output);
    if (entity === 'describe') return await handleDescribe(dapi, verb ?? rest[0], output);
    if (entity === 'healthcheck') return await handleHealthcheck(dapi, argv, output);
    if (entity === 'sync') return await handleSync(dapi, verb, rest, argv, output);
    if (entity === 'functions' && verb === 'run') return await handleFuncRun(dapi, rest, argv, output);
    if (entity === 'functions' && verb === 'list') return await handleFunctionsList(dapi, argv, limit, offset, filter, output);
    if (entity === 'files' && verb === 'list') {
      const path = rest[0] ?? '';
      const result = await dapi.files.list(path, recursive);
      printOutput(result, output);
      return true;
    }
    if (entity === 'files' && verb === 'get') {
      const result = await dapi.files.get(rest[0]);
      printOutput(result, output);
      return true;
    }
    if (entity === 'files' && verb === 'delete') {
      await dapi.files.delete(rest[0]);
      if (output !== 'quiet') console.log(`Deleted ${rest[0]}`);
      return true;
    }
    if (entity === 'files' && verb === 'put') return await handleFilesPut(dapi, rest, output);
    if (entity === 'shares' && verb === 'add') return await handleSharesAdd(dapi, rest, argv, output);
    if (entity === 'shares' && verb === 'list') return await handleSharesList(dapi, rest, output);
    if (entity === 'users' && verb === 'save') return await handleUserSave(dapi, argv, output);
    if (entity === 'groups' && verb === 'save') return await handleGroupSave(dapi, argv, output);
    if (entity === 'connections' && verb === 'save') return await handleConnSave(dapi, argv, output);
    if (entity === 'connections' && verb === 'test') return await handleConnTest(dapi, rest, argv, output);
    if (entity === 'groups' && verb === 'add-members') return await handleGroupAddMembers(dapi, rest, argv, output);
    if (entity === 'groups' && verb === 'remove-members') return await handleGroupRemoveMembers(dapi, rest, argv, output);
    if (entity === 'groups' && verb === 'list-members') return await handleGroupListMembers(dapi, rest, argv, output);
    if (entity === 'groups' && verb === 'list-memberships') return await handleGroupListMemberships(dapi, rest, argv, output);
    if (entity === 'users' && verb === 'block') return await handleUserBlock(dapi, rest, output);
    if (entity === 'users' && verb === 'unblock') return await handleUserUnblock(dapi, rest, output);
    if (entity === 'tables' && verb === 'download') return await handleTablesDownload(dapi, rest, argv, output);
    if (entity === 'tables' && verb === 'upload') return await handleTablesUpload(dapi, rest, output);
    if (entity === 'packages' && verb === 'install') return await handlePackagesInstall(dapi, rest, argv, output);
    if (entity === 'packages' && verb === 'uninstall') return await handlePackagesUninstall(dapi, rest, output);
    if (entity === 'packages' && verb === 'share') return await handlePackagesShare(dapi, rest, argv, output);
    if (entity === 'packages' && verb === 'versions') return await handlePackagesVersions(dapi, rest, output);
    if (entity === 'packages' && verb === 'set-version') return await handlePackagesSetVersion(dapi, rest, output);
    if (entity === 'packages' && verb === 'outdated') return await handlePackagesOutdated(dapi, output);
    if (entity === 'packages' && verb === 'update') return await handlePackagesUpdate(dapi, rest, argv, output);

    const source = (dapi as any)[entity];
    if (!source || !ENTITIES.includes(entity)) {
      printError(new Error(`Unknown entity type: '${entity}'. Valid: ${ENTITIES.join(', ')}`));
      return false;
    }
    if (!verb) {
      console.log(`Usage: grok s ${entity} <${VERBS.join('|')}> [args]`);
      return true;
    }

    if (verb === 'list') {
      const results = source instanceof InternalDataSource
        ? await source.list({text: filter || undefined, limit, page: Math.floor(offset / limit) + 1})
        : await source.filter(filter).by(limit).page(Math.floor(offset / limit)).list();
      printOutput(results, output);
      return true;
    }
    if (verb === 'get') {
      if (!rest[0]) { printError(new Error('Usage: grok s <entity> get <id>')); return false; }
      // The internal routers find by id only; a name has to be resolved through /entities first.
      const id = source instanceof InternalDataSource && !isUuid(rest[0])
        ? (await resolveEntity(dapi, rest[0], ENTITY_TYPES[entity])).id
        : rest[0];
      const result = await source.find(id);
      printOutput(result, output);
      return true;
    }
    if (verb === 'delete') {
      if (!rest[0]) { printError(new Error('Usage: grok s <entity> delete <id>')); return false; }
      await source.delete(rest[0]);
      if (output !== 'quiet') console.log(`Deleted ${entity}/${rest[0]}`);
      return true;
    }

    const extraVerbs = entity === 'packages' ? ', install, uninstall, update, outdated, versions, set-version, share' : '';
    printError(new Error(`Unknown verb: '${verb}'. Valid: ${VERBS.join(', ')}${extraVerbs}`));
    return false;
  } catch (err: any) {
    // A runtime failure is not a usage error: report it and exit non-zero without
    // making grok.js dump the help block (which it does for every `false` result).
    printError(err, {verbose: !!argv.verbose});
    process.exitCode = 1;
    return true;
  }
}

async function handleFilesPut(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  const [localPath, remotePath] = rest;
  if (!localPath || !remotePath) {
    printError(new Error('Usage: grok s files put <local-path> <remote-path>\n  e.g. grok s files put ./smiles.csv "System:DemoFiles/smiles.csv"'));
    return false;
  }
  if (!fs.existsSync(localPath)) {
    printError(new Error(`Local file not found: ${localPath}`));
    return false;
  }
  const result = await dapi.files.put(localPath, remotePath);
  if (output === 'quiet') return true;
  if (output === 'json') printOutput(result, output);
  else console.log(`Uploaded ${result.size} bytes to ${result.path}`);
  return true;
}

async function handleUserBlock(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!rest[0]) {
    printError(new Error('Usage: grok s users block <id-or-login>'));
    return false;
  }
  const user = await dapi.users.find(rest[0]);
  await dapi.users.block(user);
  if (output !== 'quiet') console.log(`Blocked ${user?.login ?? rest[0]}`);
  return true;
}

async function handleUserUnblock(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!rest[0]) {
    printError(new Error('Usage: grok s users unblock <id-or-login>'));
    return false;
  }
  const user = await dapi.users.find(rest[0]);
  await dapi.users.unblock(user);
  if (output !== 'quiet') console.log(`Unblocked ${user?.login ?? rest[0]}`);
  return true;
}

async function handleTablesDownload(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const name = rest[0];
  if (!name) {
    printError(new Error('Usage: grok s tables download <name-or-id> [-O <file>]'));
    return false;
  }
  const csv = await dapi.tables.download(name);
  const outFile: string | undefined = argv['output-file'] ?? argv.O;
  if (outFile) {
    fs.writeFileSync(outFile, csv);
    if (output !== 'quiet') console.log(`Wrote ${csv.length} bytes to ${outFile}`);
  }
  else
    process.stdout.write(csv);
  return true;
}

async function installPackages(dapi: NodeDapi, targets: {name: string; version: string}[], output: OutputFormat): Promise<PackageOpResult[]> {
  const results: PackageOpResult[] = [];
  for (const t of targets) {
    let name = t.name;
    try {
      const pkg = await dapi.packages.resolve(t.name);
      if (pkg) name = pkg.name;
      if (output === 'table') console.log(`Installing ${name} (${t.version})...`);
      const id = await dapi.packages.install(name, t.version);
      // the server returns the published id for any existing explicit version
      // (even an already-current one), so null means the version doesn't exist
      if (!id && t.version !== 'latest')
        results.push({package: name, version: t.version, status: 'error', error: `Version '${t.version}' not found`});
      else
        results.push({package: name, version: t.version, status: id ? 'installed' : 'noop', id: id ?? undefined});
    }
    catch (err: any) {
      results.push({package: name, version: t.version, status: 'error', error: (err?.message ?? String(err)).replace(/\s*\n\s*/g, '; ')});
    }
  }
  return results;
}

function printPackageOps(results: PackageOpResult[], output: OutputFormat): void {
  if (output === 'quiet') {
    for (const r of results)
      if (r.id) console.log(r.id);
  }
  else
    printOutput(results, output);
  if (results.some((r) => r.status === 'error')) process.exitCode = 1;
}

async function handlePackagesInstall(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  if (!rest.length) {
    printError(new Error('Usage: grok s packages install <name> [<name> ...] [--version <v>]'));
    return false;
  }
  if (argv.version !== undefined && rest.length > 1) {
    printError(new Error('--version applies to a single package'));
    return false;
  }
  const version = argv.version !== undefined ? String(argv.version) : 'latest';
  const results = await installPackages(dapi, rest.map((r) => ({name: String(r), version})), output);
  printPackageOps(results, output);
  return true;
}

async function handlePackagesUninstall(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!rest[0]) {
    printError(new Error('Usage: grok s packages uninstall <name-or-id>'));
    return false;
  }
  const result = await dapi.packages.uninstall(rest[0]);
  if (output === 'json' || output === 'csv') printOutput(result, output);
  else if (output === 'quiet') console.log(result.id);
  else
    console.log(result.repoBacked
      ? `Uninstalled ${result.name} (package entry kept — reinstall with 'grok s packages install ${result.name}')`
      : `Deleted local package ${result.name}`);
  return true;
}

async function handlePackagesShare(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const [name, ...groupArgs] = rest;
  if (!name || !groupArgs.length) {
    printError(new Error('Usage: grok s packages share <name-or-id> <group>[,<group>...] [--access View|Edit]'));
    return false;
  }
  const groups = groupArgs.flatMap((g) => g.split(',')).map((g) => g.trim()).filter(Boolean).join(',');
  const access = typeof argv.access === 'string' ? argv.access : 'View';
  if (access !== 'View' && access !== 'Edit') {
    printError(new Error(`Invalid --access '${access}'. Use 'View' or 'Edit'.`));
    return false;
  }
  const pkg = await dapi.packages.resolve(name);
  if (!pkg) {
    printError(new Error(`Package '${name}' not found`));
    return false;
  }
  const result = await dapi.shares.share(pkg.id, groups, access);
  printOutput(result, output);
  if (result?.status === 'failed') process.exitCode = 1;
  return true;
}

async function handlePackagesVersions(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!rest[0]) {
    printError(new Error('Usage: grok s packages versions <name-or-id>'));
    return false;
  }
  const versions = await dapi.packages.versions(rest[0]);
  const flat = versions.map((v: any) => ({
    version: v?.version ?? '',
    current: v?.isCurrent ?? false,
    latest: v?.isLatest ?? false,
    local: v?.isLocal ?? false,
    remote: v?.isRemote ?? false,
    debug: v?.debug ?? false,
    id: v?.id ?? '',
  }));
  printOutput(flat, output);
  return true;
}

async function handlePackagesSetVersion(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!rest[0] || rest[1] === undefined) {
    printError(new Error('Usage: grok s packages set-version <name> <version>'));
    return false;
  }
  const version = String(rest[1]);
  const pkg = await dapi.packages.resolve(String(rest[0]));
  const name = pkg?.name ?? String(rest[0]);
  const id = await dapi.packages.install(name, version);
  if (!id) {
    printError(new Error(`Version '${version}' not found for package '${name}'`));
    return false;
  }
  if (output === 'quiet') console.log(id);
  else if (output === 'json' || output === 'csv') printOutput({package: name, version, id}, output);
  else console.log(`${name} -> ${version} (published package ${id})`);
  return true;
}

async function handlePackagesOutdated(dapi: NodeDapi, output: OutputFormat): Promise<boolean> {
  const rows = await dapi.packages.outdated();
  if (!rows.length && output === 'table') {
    console.log('All packages up to date.');
    return true;
  }
  printOutput(rows, output);
  return true;
}

async function handlePackagesUpdate(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const all = argv.all === true;
  if (!all && !rest.length) {
    printError(new Error('Usage: grok s packages update <name> [<name> ...] | --all'));
    return false;
  }
  const outdated = await dapi.packages.outdated();
  const results: PackageOpResult[] = [];
  // a package tracking 'latest' keeps that auto-update intent; a pinned one gets
  // pinned to the concrete registry-latest version
  const targetVersion = (o: {desiredVersion: string; latest: string}) => o.desiredVersion === 'latest' ? 'latest' : o.latest;
  let targets: {name: string; version: string}[];
  if (all) {
    targets = outdated.map((o) => ({name: o.name, version: targetVersion(o)}));
    if (!targets.length && output === 'table') {
      console.log('All packages up to date.');
      return true;
    }
  }
  else {
    targets = [];
    for (const r of rest) {
      const raw = String(r);
      let pkg: any = null;
      try {
        pkg = await dapi.packages.resolve(raw);
      }
      catch (err: any) {
        results.push({package: raw, version: '', status: 'error', error: (err?.message ?? String(err)).replace(/\s*\n\s*/g, '; ')});
        continue;
      }
      if (!pkg) {
        results.push({package: raw, version: '', status: 'error', error: 'Package not found'});
        continue;
      }
      const o = outdated.find((x) => x.name === pkg.name);
      if (o) {
        targets.push({name: pkg.name, version: targetVersion(o)});
        continue;
      }
      const versions: any[] = Array.isArray(pkg.publishedVersions) ? pkg.publishedVersions : [];
      const installed = versions.some((v) => v?.isCurrent && !v?.debug);
      results.push({package: pkg.name, version: '', status: 'noop', note: installed ? 'already up to date' : 'not installed'});
    }
  }
  results.push(...await installPackages(dapi, targets, output));
  printPackageOps(results, output);
  return true;
}

async function handleTablesUpload(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  const [name, localPath] = rest;
  if (!name || !localPath) {
    printError(new Error('Usage: grok s tables upload <name> <file.csv|file.d42>'));
    return false;
  }
  if (!fs.existsSync(localPath)) {
    printError(new Error(`Local file not found: ${localPath}`));
    return false;
  }
  const ct = path.extname(localPath).toLowerCase() === '.d42'
    ? 'application/octet-stream' : 'text/csv';
  const result = await dapi.tables.upload(name, localPath, ct);
  if (output === 'quiet') console.log(result?.ID ?? result?.id ?? '');
  else printOutput(result, output);
  return true;
}

/**
 * Normalize common type aliases to the server's `source` field. The server uses
 * four discriminators: `function` (standalone), `function-package` (bundled in a
 * plugin package), `script`, `data-query`. `--type function` broadens to both
 * standalone and package funcs since that matches user intent.
 */
function funcTypeClause(t: string): string {
  const v = t.toLowerCase();
  if (v === 'query' || v === 'data-query' || v === 'dataquery') return 'source="data-query"';
  if (v === 'script') return 'source="script"';
  if (v === 'package' || v === 'package-function' || v === 'packagefunc') return 'source="function-package"';
  if (v === 'function' || v === 'func') return '(source="function" or source="function-package")';
  return `source="${v}"`;
}

async function handleFunctionsList(dapi: NodeDapi, argv: any, limit: number, offset: number,
                                   userFilter: string, output: OutputFormat): Promise<boolean> {
  const typeArg: string | undefined = argv.type;
  const language: string | undefined = argv.language;
  const pkg: string | undefined = argv.package;

  // --language only makes sense for scripts; if the user left --type off, imply 'script'.
  const effectiveType = typeArg ?? (language ? 'script' : undefined);

  const clauses: string[] = [];
  if (effectiveType) clauses.push(funcTypeClause(effectiveType));
  if (pkg) clauses.push(`package.shortName="${pkg}"`);
  if (userFilter) clauses.push(`(${userFilter})`);
  const composed = clauses.join(' and ');

  // Pull a generous batch; the server's public functions endpoint ignores limit/page
  // today, so pagination is done client-side below.
  let results: any[] = await dapi.functions.filter(composed).by(10000).list();

  // `language` lives on the Script subclass and isn't queryable via smart filter —
  // post-filter in Node. Same for any other subclass-only attribute we add later.
  if (language) {
    const want = language.toLowerCase();
    results = results.filter((f) => (f?.language ?? '').toLowerCase() === want);
  }

  if (offset > 0) results = results.slice(offset);
  if (limit > 0) results = results.slice(0, limit);

  printOutput(results, output);
  return true;
}

async function handleRaw(dapi: NodeDapi, method: string | undefined, rest: string[], output: OutputFormat): Promise<boolean> {
  if (!method || !rest[0]) {
    printError(new Error('Usage: grok s raw <METHOD> <path>'));
    return false;
  }
  const path = rest[0];
  const result = await dapi.raw(method, path);
  printOutput(result, output);
  return true;
}

async function handleSync(dapi: NodeDapi, subject: string | undefined, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  // `grok s sync` subcommands. Mirrors the API surface in
  // core/server/datlas/lib/src/routers/sync.dart — pairs / setups /
  // runs are read-only here; `run` triggers an actual push and prints
  // the per-item summary. Full reference in
  // core/docs/plans/instance-sync.md (Phase 7 / operational polish).
  //
  // `_callSync` does dual-path routing: tries `/api/sync/...` first
  // (nginx-fronted deployments) and falls back to `/sync/...` (bare
  // datlas on :8082). Same approach as the server-side handshake
  // code so the CLI works against either layout.
  const verb = rest[0];
  const callSync = async (method: string, path: string, body?: any): Promise<any> => {
    for (const prefix of ['/api', '']) {
      const r: any = body !== undefined
        ? await dapi.raw(method, `${prefix}${path}`, body)
        : await dapi.raw(method, `${prefix}${path}`);
      // raw() returns `null` for 404, but the server returns HTML for
      // 404 too — treat anything that isn't a sync-shaped object/array
      // as a miss and fall through to the alternate prefix.
      if (r && (Array.isArray(r) || typeof r === 'object') && r['#type'] !== 'ApiError') return r;
    }
    return null;
  };
  if (subject === 'pairs' && verb === 'list') {
    const status = argv.status ? `?status=${encodeURIComponent(argv.status)}` : '';
    const pairs = await callSync('GET', `/sync/pairs${status}`);
    printOutput(pairs, output);
    return true;
  }
  if (subject === 'setups' && verb === 'list') {
    const pairId = argv.pair ?? rest[1];
    if (!pairId) { printError(new Error('Usage: grok s sync setups list --pair <pair-id>')); return false; }
    const setups = await callSync('GET', `/sync/pairs/${encodeURIComponent(pairId)}/setups`);
    printOutput(setups, output);
    return true;
  }
  if (subject === 'setup' && verb === 'get') {
    const id = rest[1];
    if (!id) { printError(new Error('Usage: grok s sync setup get <setup-id>')); return false; }
    const setup = await callSync('GET', `/sync/setups/${encodeURIComponent(id)}`);
    printOutput(setup, output);
    return true;
  }
  if (subject === 'run') {
    // grok s sync run <setup-id>
    const id = rest[1] ?? verb;
    if (!id) { printError(new Error('Usage: grok s sync run <setup-id>')); return false; }
    const result = await callSync('POST', `/sync/setups/${encodeURIComponent(id)}/run`, {});
    if (output === 'json' || output === 'csv') { printOutput(result, output); return true; }
    if (output === 'quiet') { console.log(result?.runId ?? ''); return true; }
    if (!result) { printError(new Error(`Sync setup ${id} not found, or server has no /sync routes.`)); return false; }
    console.log(`Run ${result?.runId} → ${result?.remoteUrl} (${result?.status})`);
    if (result?.counts) {
      const cs = Object.keys(result.counts).filter(k => (result.counts[k] ?? 0) > 0)
        .map(k => `${result.counts[k]} ${k}`).join(', ');
      if (cs) console.log(`  ${cs}`);
    }
    const items = Array.isArray(result?.items) ? result.items : [];
    if (items.length) {
      console.log('\nItems:');
      printOutput(items, 'table');
    }
    return true;
  }
  printError(new Error(
    'Usage:\n' +
    '  grok s sync pairs list [--status active|pending|revoked]\n' +
    '  grok s sync setups list --pair <pair-id>\n' +
    '  grok s sync setup get <setup-id>\n' +
    '  grok s sync run <setup-id>'));
  return false;
}

async function handleHealthcheck(dapi: NodeDapi, argv: any, output: OutputFormat): Promise<boolean> {
  const module: string | undefined = argv.module;
  const path = module
    ? `/api/public/v1/healthcheck?module=${encodeURIComponent(module)}`
    : '/api/public/v1/healthcheck';
  const result = await dapi.raw('GET', path);

  if (output === 'json') {
    printOutput(result, output);
    return true;
  }

  if (output !== 'quiet') {
    console.log(`status:  ${result?.status ?? ''}`);
    console.log(`server:  ${result?.server ?? ''}`);
    console.log(`version: ${result?.version ?? ''}`);
    console.log(`time:    ${result?.time ?? ''}`);
    console.log('');
  }
  printOutput(result?.services ?? [], output);
  return true;
}

async function handleDescribe(dapi: NodeDapi, entityType: string | undefined, output: OutputFormat): Promise<boolean> {
  if (!entityType) {
    printError(new Error('Usage: grok s describe <entity-type>'));
    return false;
  }
  const result = await dapi.describe(entityType);
  printOutput(result, output);
  return true;
}

async function handleFuncRun(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  let funcName = rest[0];
  if (!funcName) {
    printError(new Error('Usage: grok s functions run <Name:function(args...)> [--json params.json]'));
    return false;
  }

  let params: Record<string, any> = {};

  if (argv.json) {
    const fs = require('fs');
    try {
      params = JSON.parse(fs.readFileSync(argv.json, 'utf8'));
    } catch (err: any) {
      printError(new Error(`Cannot read params file: ${err.message}`));
      return false;
    }
  } else {
    const parsed = parseFuncCall(funcName);
    funcName = parsed.name;
    params = parsed.params;
  }

  const result = await dapi.functions.run(funcName, params);
  printOutput(result, output);
  return true;
}

function readJsonFile(jsonPath: string): any {
  return JSON.parse(fs.readFileSync(jsonPath, 'utf8'));
}

function readJsonBody(argv: any, entity: string): any {
  if (!argv.json)
    throw new Error(`Usage: grok s ${entity} save --json ${entity}.json`);
  try {
    return readJsonFile(argv.json);
  } catch (err: any) {
    throw new Error(`Cannot read ${entity} file '${argv.json}': ${err.message}`);
  }
}

async function handleSharesAdd(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const [entity, ...groupArgs] = rest;
  if (!entity || !groupArgs.length) {
    printError(new Error('Usage: grok s shares add <entity-id-or-name> <group>[,<group>...] [--access View|Edit]'));
    return false;
  }
  const groups = groupArgs.flatMap((g) => g.split(',')).map((g) => g.trim()).filter(Boolean).join(',');
  const access = typeof argv.access === 'string' ? argv.access : 'View';
  if (access !== 'View' && access !== 'Edit') {
    printError(new Error(`Invalid --access '${access}'. Use 'View' or 'Edit'.`));
    return false;
  }
  const result = await dapi.shares.share(entity, groups, access);
  printOutput(result, output);
  if (result?.status === 'failed') process.exitCode = 1;
  return true;
}

async function handleSharesList(dapi: NodeDapi, rest: string[], output: OutputFormat): Promise<boolean> {
  const entityId = rest[0];
  if (!entityId) {
    printError(new Error('Usage: grok s shares list <entity-id>   (entity id must be a UUID)'));
    return false;
  }
  const perms = await dapi.shares.list(entityId);
  const flat = (Array.isArray(perms) ? perms : []).map((p: any) => ({
    group: p?.userGroup?.friendlyName ?? p?.userGroup?.name ?? p?.userGroup?.id ?? '',
    groupId: p?.userGroup?.id ?? '',
    access: p?.permission?.name ?? p?.permission?.friendlyName ?? '',
    personal: p?.userGroup?.personal ?? false,
  }));
  printOutput(flat, output);
  return true;
}

async function handleUserSave(dapi: NodeDapi, argv: any, output: OutputFormat): Promise<boolean> {
  let body: any;
  try { body = readJsonBody(argv, 'users'); }
  catch (err: any) { printError(err); return false; }
  const result = await dapi.users.save(body);
  printOutput(result, output);
  return true;
}

async function handleGroupSave(dapi: NodeDapi, argv: any, output: OutputFormat): Promise<boolean> {
  let body: any;
  try { body = readJsonBody(argv, 'groups'); }
  catch (err: any) { printError(err); return false; }
  const saveRelations = argv['save-relations'] === true || argv.saveRelations === true;
  const result = await dapi.groups.save(body, saveRelations);
  printOutput(result, output);
  return true;
}

async function handleConnSave(dapi: NodeDapi, argv: any, output: OutputFormat): Promise<boolean> {
  if (!argv.json) {
    printError(new Error('Usage: grok s connections save --json conn.json [--save-credentials]'));
    return false;
  }
  let body: any;
  try {
    body = readJsonFile(argv.json);
  } catch (err: any) {
    printError(new Error(`Cannot read connection file '${argv.json}': ${err.message}`));
    return false;
  }
  const saveCredentials = argv['save-credentials'] === true || argv.saveCredentials === true;
  const result = await dapi.connections.save(body, saveCredentials);
  printOutput(result, output);
  return true;
}

async function handleConnTest(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  let body: any;
  if (argv.json) {
    try {
      body = readJsonFile(argv.json);
    } catch (err: any) {
      printError(new Error(`Cannot read connection file '${argv.json}': ${err.message}`));
      return false;
    }
  } else if (rest[0]) {
    body = await dapi.connections.find(rest[0]);
  } else {
    printError(new Error('Usage: grok s connections test <id-or-name> | --json conn.json'));
    return false;
  }
  await dapi.connections.test(body);
  if (output !== 'quiet') console.log('ok');
  return true;
}

async function handleGroupAddMembers(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const [group, ...members] = rest;
  if (!group || !members.length) {
    printError(new Error('Usage: grok s groups add-members <group> <member> [<member> ...] [--admin] [--user]'));
    return false;
  }
  const isAdmin = argv.admin === true;
  const personalOnly = argv.user === true;
  const results = await dapi.groups.addMembers(group, members, isAdmin, personalOnly);
  printOutput(results, output);
  const anyError = results.some((r) => r.status === 'error');
  if (anyError) process.exitCode = 1;
  return true;
}

async function handleGroupRemoveMembers(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const [group, ...members] = rest;
  if (!group || !members.length) {
    printError(new Error('Usage: grok s groups remove-members <group> <member> [<member> ...] [--user]'));
    return false;
  }
  const personalOnly = argv.user === true;
  const results = await dapi.groups.removeMembers(group, members, personalOnly);
  printOutput(results, output);
  const anyError = results.some((r) => r.status === 'error');
  if (anyError) process.exitCode = 1;
  return true;
}

async function handleGroupListMembers(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  if (!rest[0]) {
    printError(new Error('Usage: grok s groups list-members <group> [--admin | --no-admin]'));
    return false;
  }
  const admin: boolean | undefined = typeof argv.admin === 'boolean' ? argv.admin : undefined;
  const result = await dapi.groups.getMembers(rest[0], admin);
  printOutput(result, output);
  return true;
}

async function handleGroupListMemberships(dapi: NodeDapi, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  if (!rest[0]) {
    printError(new Error('Usage: grok s groups list-memberships <group> [--admin | --no-admin]'));
    return false;
  }
  const admin: boolean | undefined = typeof argv.admin === 'boolean' ? argv.admin : undefined;
  const result = await dapi.groups.getMemberships(rest[0], admin);
  printOutput(result, output);
  return true;
}

async function handleBatch(dapi: NodeDapi, argv: any, verb: string | undefined, rest: string[], output: OutputFormat): Promise<boolean> {
  let request: BatchRequest;

  if (verb?.endsWith('.json') && rest.length === 0) {
    // grok s batch manifest.json
    let raw: any;
    try {
      raw = JSON.parse(fs.readFileSync(verb, 'utf8'));
    } catch (err: any) {
      printError(new Error(`Cannot read manifest file '${verb}': ${err.message}`));
      return false;
    }
    request = resolveManifestSources(raw);
  } else if (verb && ENTITIES.includes(verb)) {
    const batchVerb: string | undefined = rest[0];
    if (!batchVerb) {
      printError(new Error(`Usage: grok s batch ${verb} <verb> [args...]`));
      return false;
    }
    if (argv.json) {
      // grok s batch users create --json users.json
      let paramsArray: any[];
      try {
        const raw = JSON.parse(fs.readFileSync(argv.json, 'utf8'));
        paramsArray = Array.isArray(raw) ? raw : [raw];
      } catch (err: any) {
        printError(new Error(`Cannot read params file '${argv.json}': ${err.message}`));
        return false;
      }
      request = buildInlineManifest(verb, batchVerb, paramsArray);
    } else {
      // grok s batch files delete path1 path2 ...
      const batchArgs = rest.slice(1);
      if (!batchArgs.length) {
        printError(new Error(`Usage: grok s batch ${verb} ${batchVerb} <arg1> [arg2 ...]`));
        return false;
      }
      request = buildInlineManifest(verb, batchVerb, batchArgs);
    }
  } else {
    printError(new Error('Usage: grok s batch <entity> <verb> [args]\n       grok s batch <entity> <verb> --json params.json\n       grok s batch manifest.json'));
    return false;
  }

  const result = await dapi.batch(request);
  printBatchOutput(result, output);
  return result.summary.failed === 0 && result.summary.partial === 0;
}

export function buildInlineManifest(entity: string, verb: string, args: any[]): BatchRequest {
  const action = `${entity}.${verb}`;
  let operations: BatchOperation[];

  if (args.length > 0 && typeof args[0] === 'object') {
    // Array of param objects (from --json array)
    operations = args.map((p, i) => ({id: `op${i}`, action, params: p}));
  } else {
    // Array of string args — map to appropriate param key
    const paramKey = entity === 'files' ? 'path' : 'id';
    operations = (args as string[]).map((arg, i) => ({
      id: `op${i}`,
      action,
      params: {[paramKey]: arg},
    }));
  }

  return {operations};
}

export function resolveManifestSources(manifest: any): BatchRequest {
  if (!manifest.operations) return manifest;
  const operations = manifest.operations.map((op: any) => {
    if (op.action === 'files.put' && op.params?.source) {
      const sourcePath: string = op.params.source;
      const content = fs.readFileSync(sourcePath).toString('base64');
      const {source: _dropped, ...rest} = op.params;
      return {...op, params: {...rest, content}};
    }
    return op;
  });
  return {...manifest, operations};
}

export function parseFuncCall(expr: string): {name: string; params: Record<string, any>} {
  const parenIdx = expr.indexOf('(');
  if (parenIdx === -1) return {name: expr, params: {}};

  const name = expr.slice(0, parenIdx);
  const inner = expr.slice(parenIdx + 1, expr.lastIndexOf(')')).trim();

  if (!inner) return {name, params: {}};

  if (inner.startsWith('{')) {
    try {
      const params = JSON.parse(inner.replace(/([{,]\s*)([a-zA-Z_]\w*)(\s*:)/g, '$1"$2"$3'));
      return {name, params};
    } catch {
      return {name, params: {}};
    }
  }

  const positional = inner.split(',').map((s) => {
    s = s.trim();
    if ((s.startsWith('"') && s.endsWith('"')) || (s.startsWith("'") && s.endsWith("'")))
      return s.slice(1, -1);
    const n = Number(s);
    return isNaN(n) ? s : n;
  });
  return {name, params: Object.fromEntries(positional.map((v, i) => [String(i), v]))};
}

export const HELP_SERVER = `
Usage: grok server <entity> <verb> [args] [options]
       grok s <entity> <verb> [args] [options]

Manage a Datagrok server from the command line.

Entities:
  users, groups, functions, connections, queries, scripts, packages, reports, files, tables

Verbs:
  list      List entities
  get       Get a single entity by ID or name
  delete    Delete an entity by ID

Special commands:
  grok s functions run <Name:func(args)>              Call a function
  grok s functions list [--type <t>] [--language <l>] [--package <p>] [--filter <expr>]
                                                      Type: script|query|function|package
                                                      Language applies to scripts (python, r, julia, nodejs, octave, grok)
  grok s files list <path> [-r]                       List files (recursive with -r)
  grok s files get <path>                             Download a file (returns bytes)
  grok s files delete <path>                          Delete a file
  grok s files put <local> <remote>                   Upload a local file
  grok s raw <METHOD> <path>                          Hit any API endpoint
  grok s describe <entity-type>                       Show entity JSON schema
  grok s healthcheck [--module <name>]                Check server + per-module health
  grok s shares add <entity> <group>[,<group>...] [--access View|Edit]
                                                      Share an entity with one or more groups
  grok s shares list <entity-id>                      List who an entity (UUID) is shared with
  grok s users save --json user.json                  Create or update a user from a JSON file
  grok s groups save --json group.json [--save-relations]
                                                      Create or update a group from a JSON file
  grok s connections save --json conn.json [--save-credentials]
                                                      Create or update a connection from a JSON file
  grok s connections test <id-or-name>                Test connectivity of an existing connection
  grok s connections test --json conn.json            Test connectivity of a connection defined in JSON
  grok s groups add-members <group> <m>... [--admin]  Add one or more users/groups as members
  grok s groups remove-members <group> <m>...         Remove members (no-op if not a member)
  grok s groups list-members <group> [--admin]        List members (optionally filter by admin)
  grok s groups list-memberships <group> [--admin]    List parent groups
  grok s users block <id-or-login>                    Block a user from the platform
  grok s users unblock <id-or-login>                  Unblock a previously blocked user
  grok s tables upload <name> <file.csv|file.d42>     Upload a CSV or d42 binary as a Datagrok table
  grok s tables download <name-or-id> [-O <file>]     Download a table as CSV (stdout by default)
  grok s packages install <name>... [--version <v>]   Install latest (or pinned) versions from the registry
  grok s packages uninstall <name>                    Uninstall a package (repository entry is kept)
  grok s packages update <name>... | --all            Upgrade to the latest registry version
  grok s packages outdated                            Show installed vs latest registry versions
  grok s packages versions <name>                     List published versions (current/latest/debug flags)
  grok s packages set-version <name> <version>        Activate a specific published version
  grok s packages share <name> <group>[,...] [--access View|Edit]
                                                      Share a package with one or more groups
  grok s pull [<nqName|id>...] --out <dir>            Export entities into a bundle directory
  grok s push <bundle-dir> [--dry-run]                Import a bundle into the target server
  grok s migrate <selection> --from <a> --to <b>      Pull from one instance and push into another
  grok s diff <bundle-dir>                            What a push would change (read-only, plans with skip)
  grok s bundle ls <bundle-dir>                       List what a bundle contains
  grok s batch <entity> <verb> arg1 [arg2 ...]        Batch operation (one round-trip)
  grok s batch <entity> <verb> --json params.json     Batch from JSON array
  grok s batch manifest.json                          Run a workflow manifest
  grok s sync pairs list [--status <s>]               List cross-instance sync pairs (status: active|pending|revoked)
  grok s sync setups list --pair <pair-id>            List the named sync setups under a pair
  grok s sync setup get <setup-id>                    Inspect a setup (selections, direction, last run)
  grok s sync run <setup-id>                          Trigger a push run; prints per-item outcome

Pull / push / migrate options:
  --out <dir>           Bundle directory to write (pull; merges into an existing bundle)
  --type <t,t>          conn, query, script, project, dashboard, space, view, layout, table,
                        file, group, job, notebook, model
  --name <glob>         Entity name glob: 'Cereal*' matches by prefix, '*demo*' anywhere
  --namespace <ns>      Everything under a namespace, recursively
  --space <nqName>      The space itself plus everything under it
  --author <login>      Entities authored by a login
  --tag <tag>           Entities carrying a tag
  --since <2w|date>     Entities updated since ('2w' means '-2w'; also --since=-2w)
  --filter <expr>       Smart-filter expression ANDed with the other flags
  --no-deps             Do not follow dependencies (never-travels rules still apply)
  --no-include-data     Do not pull table data (.d42) for pulled tables
  --include-files       Also pull the bytes of pulled files
  --replace             Clear the bundle directory before writing
  --dry-run             Push: print the plan and stop
  --on-conflict <p>     Push: what to do when the name is taken on the target by another id —
                        fail (default) aborts before any write, skip leaves the target alone
                        (and fails whatever depends on it), adopt writes into the twin and
                        records idmap.json, duplicate creates the bundle entity under its own
                        id next to the twin (a UserGroup cannot be duplicated — group names
                        are unique, so that row fails)
  --creds <file.yaml>   Push: target-side connection secrets, 'Ns:Conn:' then '  password: \${VAR}'
                        (a covered connection is always written, so a push rotates the secret)
  --from <alias|url>    Migrate: source instance (pulled from, never written to)
  --to <alias|url>      Migrate: target instance
  --keep                Migrate: keep the temporary bundle and print its path on stderr

Options:
  --host <alias|url>    Server alias from config or full URL
  --output <format>     Output format: table (default), json, csv, quiet
  --filter <text>       Smart filter expression
  --limit <n>           Page size (default: 50)
  --offset <n>          Start offset (default: 0)
  -r, --recursive       Recursive (for files list)
  --verbose             Print the stack of a runtime failure, not just its message
  --json <file>         Read function parameters or batch params from JSON file
  -O, --output-file     Write table download to a file instead of stdout
  --type <t>            Function discriminator: script | query | function | package
  --language <lang>     Script language: python, r, julia, nodejs, octave, grok
  --package <name>      Restrict to functions belonging to a package (by short name)
  --version <v>         Package version for 'packages install' (default: latest)
  --all                 For 'packages update': upgrade every outdated package

Batch manifest options (in manifest.json):
  stopOnError           Stop on first failure (default: true)
  transaction           Wrap DB ops in transaction (default: false)
  concurrency           Accepted (always treated as 1)

Examples:
  grok s users list
  grok s users list --output json --limit 10
  grok s users save --json user.json
  grok s groups save --json group.json --save-relations
  grok s shares add "JohnDoe:MyConnection" Chemists,Admins --access Edit
  grok s shares list <entity-uuid>
  grok s connections list --filter "PostgreSQL"
  grok s connections get <id>
  grok s connections delete <id>
  grok s connections save --json conn.json --save-credentials
  grok s connections test "JohnDoe:MyConnection"
  grok s connections test --json conn.json
  grok s functions run 'Chem:smilesToMw("ccc")'
  grok s functions run Chem:test --json params.json
  grok s functions list --type script --language python --limit 20
  grok s functions list --package Chem --type query
  grok s functions list --package PowerPack --type package --filter 'name contains "grid"'
  grok s files list "System:AppData" -r
  grok s files put ./smiles.csv "System:DemoFiles/smiles.csv"
  grok s packages install Chem Bio PowerGrid
  grok s packages install Chem --version 1.14.0
  grok s packages update --all
  grok s packages versions Chem
  grok s packages share Chem Chemists --access View
  grok s raw GET /api/users/current
  grok s describe connections
  grok s users list --host dev
  grok s users list --host "https://mygrok.com/api"
  grok s groups add-members Admins alice bob --admin
  grok s groups remove-members Admins alice
  grok s groups list-members Admins --admin
  grok s groups list-memberships alice
  grok s pull Admin:CerealDemog --out ./bundle --host dev
  grok s pull Chemists --out ./bundle --host dev
  grok s pull --type script --author alice --no-deps --out ./bundle
  grok s bundle ls ./bundle
  grok s push ./bundle --host local --dry-run
  grok s push ./bundle --host local --on-conflict adopt
  grok s pull --type space --namespace MySpace --include-files --out ./bundle
  grok s push ./bundle --host prod --creds ./creds.yaml
  grok s migrate Chem:TargetDashboard --from dev --to prod --dry-run
  grok s diff ./bundle --host local
  grok s batch files delete "System:AppData/old.txt" "System:DemoFiles/tmp.txt"
  grok s batch users create --json users.json
  grok s batch manifest.json
  grok s batch files delete "System:AppData/a" "System:DemoFiles/b" --output json
`;
