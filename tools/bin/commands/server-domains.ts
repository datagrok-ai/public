/// `grok s domains ...` — entity-mapped domain schemas and their rows (DomainsRouter, `/domains/`).
/// Every verb addresses a schema as `<schema>` or a table as `<schema>.<table>`.
import * as fs from 'fs';
import * as path from 'path';
import {NodeDapi, DomainAddress, parseDomainAddress} from '../utils/node-dapi';
import {printOutput, printError, csvCell, OutputFormat} from '../utils/server-output';

const APPLY_KEYS = ['tables', 'extend', 'propertySchemas', 'dropTables', 'ifVersion', 'confirmDestructive'];
const PERMISSIONS = ['View', 'Edit', 'Delete', 'Share', 'Extend'];

const USAGE = `Usage: grok s domains <verb> [args]
  list [<schema>] [--filter <text>]                      Schemas, or the tables of one schema
  get <schema> | <schema.table> [<row-id>]               Manifest, a table's columns, or one row
  query <schema.table> [--filter <expr>] [--columns a,b] [--sort 'a,!b'] [--expand x] [--limit n] [--offset n]
  count <schema.table> [--filter <expr>]
  insert <schema.table> --json rows.json | col=value ... [--error-on-duplicate]
  update <schema.table> <row-id> --json values.json | col=value ... [--version n]
  delete <schema.table> <row-id> | --filter <expr> [--limit n]
  delete <schema> --force                                Purge a user-managed schema (data included)
  upload <schema.table> <file.csv|.d42|.json> [--upsert] [--no-all-or-nothing] [--error-on-duplicate]
  download <schema.table> [-O out.csv|out.d42] [--filter <expr>] [--columns a,b] [--sort 'a,!b'] [--limit n]
  aggregate <schema.table> --measures 'count,sum(x) as total' [--group-by a,b] [--filter <expr>] [--sort s] [--limit n]
  aggregate <schema.table> --json spec.json
  transaction <schema> --json ops.json
  audit <schema> | <schema.table> [<row-id>] [--limit n]
  capabilities <schema.table>
  grants <schema> | <schema.table>
  grant <schema> | <schema.table> <group>[,<group>...] [--access View|Edit|Delete|Share|Extend]
  revoke <schema> | <schema.table> <group>[,<group>...] [--access <permission>]
  create <name> [--friendly-name <text>] [--description <text>]
  apply <schema> --json manifest.json [--dry-run] [--confirm-destructive] [--if-version <v>]`;

export async function handleDomains(dapi: NodeDapi, verb: string | undefined, rest: string[], argv: any,
                                    output: OutputFormat): Promise<boolean> {
  const args = rest.map(String);
  const domains = dapi.domains;
  const filter: string = argv.filter ?? argv.f ?? '';
  const limit: number | undefined = argv.limit ?? argv.l;
  const offset: number = Number(argv.offset ?? 0);

  switch (verb) {
    case 'list': {
      if (!args[0]) {
        const schemas = await domains.schemas(filter);
        printOutput(schemas.map(schemaRow), output);
        return true;
      }
      const schema = await domains.schema(parseDomainAddress(args[0], {table: false}).schema);
      printOutput((schema.tables ?? []).map(tableRow), output);
      return true;
    }
    case 'get': {
      if (!args[0]) return usage('get <schema> | <schema.table> [<row-id>]');
      const a = parseDomainAddress(args[0]);
      if (!a.table) {
        printJson(await domains.manifest(a.schema), output, (m) => m?.name);
        return true;
      }
      if (args[1]) {
        const row = await domains.getRow(a.schema, a.table, args[1]);
        if (!row) throw new Error(`Row ${args[1]} not found in ${args[0]}`);
        printOutput(row, output);
        return true;
      }
      const manifest = await domains.manifest(a.schema);
      const table = manifest?.tables?.[a.table];
      if (!table) throw new Error(`Domain table '${args[0]}' not found`);
      if (output === 'json') printOutput(table, output);
      else printOutput(columnRows(table), output);
      return true;
    }
    case 'query': {
      if (!args[0]) return usage('query <schema.table> [--filter <expr>] [--columns a,b] [--sort s] [--limit n]');
      const a = parseDomainAddress(args[0], {table: true});
      const rows = await domains.query(a.schema, a.table!, querySpec(argv, filter, limit ?? 50, offset));
      printOutput(rows, output);
      return true;
    }
    case 'count': {
      if (!args[0]) return usage('count <schema.table> [--filter <expr>]');
      const a = parseDomainAddress(args[0], {table: true});
      const count = await domains.count(a.schema, a.table!, filter || undefined);
      if (output === 'json' || output === 'csv') printOutput({count}, output);
      else console.log(count);
      return true;
    }
    case 'insert': {
      if (!args[0]) return usage('insert <schema.table> --json rows.json | col=value ...');
      const a = parseDomainAddress(args[0], {table: true});
      const rows = argv.json ? readJson(argv.json) : parseKeyValues(args.slice(1));
      if (!argv.json && !Object.keys(rows).length) return usage('insert <schema.table> --json rows.json | col=value ...');
      printOutput(await domains.insert(a.schema, a.table!, rows, argv['error-on-duplicate'] === true), output);
      return true;
    }
    case 'update': {
      if (!args[0] || !args[1]) return usage('update <schema.table> <row-id> --json values.json | col=value ... [--version n]');
      const a = parseDomainAddress(args[0], {table: true});
      const values = argv.json ? readJson(argv.json) : parseKeyValues(args.slice(2));
      if (!argv.json && !Object.keys(values).length) return usage('update <schema.table> <row-id> --json values.json | col=value ...');
      const version = argv.version !== undefined ? Number(argv.version) : undefined;
      printOutput(await domains.update(a.schema, a.table!, args[1], values, version), output);
      return true;
    }
    case 'delete': {
      if (!args[0]) return usage('delete <schema.table> <row-id> | delete <schema.table> --filter <expr> | delete <schema> --force');
      const a = parseDomainAddress(args[0]);
      if (!a.table) {
        if (argv.force !== true)
          throw new Error(`Deleting schema '${a.schema}' purges every table and row in it; re-run with --force to confirm`);
        await domains.deleteSchema(a.schema);
        if (output !== 'quiet') console.log(`Deleted schema ${a.schema}`);
        return true;
      }
      if (args[1]) {
        await domains.deleteRow(a.schema, a.table, args[1]);
        if (output !== 'quiet') console.log(`Deleted ${args[0]}/${args[1]}`);
        return true;
      }
      if (!filter) return usage('delete <schema.table> <row-id> | delete <schema.table> --filter <expr> [--limit n]');
      const report = await domains.deleteWhere(a.schema, a.table, filter, limit);
      if (output === 'json' || output === 'csv') printOutput(report, output);
      else if (output !== 'quiet') console.log(`Deleted ${report.deleted} rows${report.hasMore ? ' (more remain — re-run to continue)' : ''}`);
      return true;
    }
    case 'upload': return await handleUpload(dapi, args, argv, output);
    case 'download': return await handleDownload(dapi, args, argv, filter, limit, offset, output);
    case 'aggregate': {
      if (!args[0]) return usage('aggregate <schema.table> --measures <m> [--group-by a,b] | --json spec.json');
      const a = parseDomainAddress(args[0], {table: true});
      const spec = argv.json ? readJson(argv.json) : aggregateSpec(argv, filter, limit);
      if (!spec.measures?.length && !argv.json) return usage('aggregate <schema.table> --measures \'count,sum(x) as total\' [--group-by a,b]');
      printOutput(await domains.aggregate(a.schema, a.table!, spec), output);
      return true;
    }
    case 'transaction': {
      if (!args[0] || !argv.json) return usage('transaction <schema> --json ops.json');
      const a = parseDomainAddress(args[0], {table: false});
      const ops = readJson(argv.json);
      printOutput(await domains.transaction(a.schema, Array.isArray(ops) ? ops : ops.ops), output);
      return true;
    }
    case 'audit': {
      if (!args[0]) return usage('audit <schema> | <schema.table> [<row-id>] [--limit n]');
      const a = parseDomainAddress(args[0]);
      const entries = !a.table ? await domains.schemaAudit(a.schema, limit)
        : args[1] ? await domains.rowAudit(a.schema, a.table, args[1])
          : await domains.tableAudit(a.schema, a.table, limit);
      printOutput(entries, output);
      return true;
    }
    case 'capabilities': {
      if (!args[0]) return usage('capabilities <schema.table>');
      const a = parseDomainAddress(args[0], {table: true});
      printOutput(flatten(await domains.capabilities(a.schema, a.table!)), output);
      return true;
    }
    case 'grants': {
      if (!args[0]) return usage('grants <schema> | <schema.table>');
      const grants = await domains.grants(await domains.entityId(parseDomainAddress(args[0])));
      printOutput((grants ?? []).map(grantRow), output);
      return true;
    }
    case 'grant':
    case 'revoke': return await handleGrant(dapi, verb, args, argv, output);
    case 'create': {
      if (!args[0]) return usage('create <name> [--friendly-name <text>] [--description <text>]');
      const res = await domains.createSchema(args[0], optString(argv['friendly-name']), optString(argv.description));
      printOutput(res, output);
      return true;
    }
    case 'apply': {
      if (!args[0] || !argv.json) return usage('apply <schema> --json manifest.json [--dry-run] [--confirm-destructive] [--if-version <v>]');
      const a = parseDomainAddress(args[0], {table: false});
      const body = applyBody(readJson(argv.json), argv);
      printJson(await domains.applySchema(a.schema, body, argv['dry-run'] === true), output, (r) => r?.version);
      return true;
    }
  }
  printError(new Error(USAGE));
  return false;
}

function usage(line: string): boolean {
  printError(new Error(`Usage: grok s domains ${line}`));
  return false;
}

function optString(v: any): string | undefined {
  return v === undefined || v === true ? undefined : String(v);
}

function readJson(file: string): any {
  try {
    return JSON.parse(fs.readFileSync(String(file), 'utf8'));
  } catch (err: any) {
    throw new Error(`Cannot read JSON file '${file}': ${err.message}`);
  }
}

/** A nested result (manifest, apply plan) is printed as JSON in every format but quiet. */
function printJson(value: any, output: OutputFormat, quietKey?: (v: any) => any): void {
  if (output === 'quiet') console.log(quietKey?.(value) ?? '');
  else printOutput(value, 'json');
}

function listFlag(v: any): string[] | undefined {
  if (v === undefined || v === true) return undefined;
  const parts = (Array.isArray(v) ? v : [v]).flatMap((x) => String(x).split(',')).map((s) => s.trim()).filter(Boolean);
  return parts.length ? parts : undefined;
}

/** Row values from `col=value` arguments; a value that parses as JSON is sent typed, anything else as a string. */
export function parseKeyValues(args: string[]): Record<string, any> {
  const values: Record<string, any> = {};
  for (const a of args) {
    const s = String(a);
    const eq = s.indexOf('=');
    if (eq <= 0) throw new Error(`Expected <column>=<value>, got '${s}'`);
    const raw = s.slice(eq + 1);
    let v: any = raw;
    try { v = JSON.parse(raw); } catch { /* a bare string */ }
    values[s.slice(0, eq)] = v;
  }
  return values;
}

/** `'count, sum(amount) as total, avg(amount)'` → aggregate measures. */
export function parseMeasures(s: string): any[] {
  return String(s).split(',').map((m) => m.trim()).filter(Boolean).map((m) => {
    const match = /^(\w+)\s*(?:\(\s*([\w.]*)\s*\))?(?:\s+as\s+(\w+))?$/i.exec(m);
    if (!match) throw new Error(`Cannot parse measure '${m}' (expected 'count', 'sum(amount)' or 'avg(amount) as mean')`);
    const [, fn, column, as] = match;
    return {fn: fn.toLowerCase(), ...(column ? {column} : {}), ...(as ? {as} : {})};
  });
}

export function querySpec(argv: any, filter: string, limit: number | undefined, offset: number): Record<string, any> {
  const spec: Record<string, any> = {};
  if (filter) spec.filter = filter;
  if (argv.sort) spec.sort = String(argv.sort);
  const columns = listFlag(argv.columns);
  if (columns) spec.columns = columns;
  const expand = listFlag(argv.expand);
  if (expand) spec.expand = expand;
  if (limit !== undefined) spec.limit = Number(limit);
  if (offset > 0) spec.offset = offset;
  return spec;
}

function aggregateSpec(argv: any, filter: string, limit: number | undefined): Record<string, any> {
  const spec: Record<string, any> = {measures: argv.measures ? parseMeasures(argv.measures) : []};
  const groupBy = listFlag(argv['group-by']);
  if (groupBy) spec.groupBy = groupBy;
  if (filter) spec.filter = filter;
  if (argv.sort) spec.sort = String(argv.sort);
  if (limit !== undefined) spec.limit = Number(limit);
  return spec;
}

/** The apply payload from a schema.json (or a partial body): only the keys the server accepts, flags folded in. */
export function applyBody(json: any, argv: any): Record<string, any> {
  const body: Record<string, any> = {};
  for (const k of APPLY_KEYS)
    if (json?.[k] !== undefined) body[k] = json[k];
  if (argv['confirm-destructive'] === true) body.confirmDestructive = true;
  if (argv['if-version'] !== undefined) body.ifVersion = String(argv['if-version']);
  return body;
}

/** Serializer omits defaults, so the table view fills them back in. */
export function schemaRow(s: any): Record<string, any> {
  return {
    name: s?.name ?? '',
    friendlyName: s?.friendlyName ?? '',
    managedBy: s?.managedBy ?? 'package',
    version: s?.version ?? '',
    tables: Array.isArray(s?.tables) ? s.tables.length : 0,
    id: s?.id ?? '',
  };
}

export function tableRow(t: any): Record<string, any> {
  return {
    name: t?.name ?? '',
    securityMode: t?.securityMode ?? 'table',
    businessKey: Array.isArray(t?.businessKey) ? t.businessKey.join(',') : '',
    nameColumn: t?.nameColumn ?? '',
    origin: t?.origin ?? 'package',
    readOnly: t?.readOnly ?? false,
    description: t?.description ?? '',
    id: t?.id ?? '',
  };
}

export function columnRows(table: any): Record<string, any>[] {
  return Object.entries(table?.columns ?? {}).map(([name, c]: [string, any]) => ({
    column: name,
    type: c?.type ?? '',
    ref: c?.ref ?? '',
    required: !!c?.required,
    unique: !!c?.unique,
    default: c?.default ?? '',
    choices: Array.isArray(c?.choices) ? c.choices.join('|') : '',
    description: c?.description ?? '',
  }));
}

function grantRow(g: any): Record<string, any> {
  return {
    group: g?.group?.friendlyName ?? g?.group?.id ?? '',
    permission: g?.permission ?? '',
    personal: g?.group?.personal ?? false,
    groupId: g?.group?.id ?? '',
  };
}

function flatten(o: any): Record<string, any> {
  const out: Record<string, any> = {};
  for (const [k, v] of Object.entries(o ?? {}))
    out[k] = Array.isArray(v) ? v.join(',') : v;
  return out;
}

/** Full-width CSV for exports: every key of every row, objects as JSON (the table printer truncates both). */
export function rowsToCsv(rows: any[]): string {
  if (!rows.length) return '';
  const keys: string[] = [];
  const seen = new Set<string>();
  for (const r of rows)
    for (const k of Object.keys(r ?? {}))
      if (!seen.has(k)) { seen.add(k); keys.push(k); }
  const cell = (v: any) => v === null || v === undefined ? '' : typeof v === 'object' ? JSON.stringify(v) : String(v);
  const lines = [keys.map(csvCell).join(',')];
  for (const r of rows)
    lines.push(keys.map((k) => csvCell(cell(r?.[k]))).join(','));
  return lines.join('\n') + '\n';
}

function contentTypeFor(file: string): string {
  const ext = path.extname(file).toLowerCase();
  if (ext === '.csv') return 'text/csv';
  if (ext === '.d42') return 'application/octet-stream';
  if (ext === '.json') return 'application/json';
  throw new Error(`Unsupported upload format '${ext || file}': use .csv, .d42 or .json`);
}

async function handleUpload(dapi: NodeDapi, args: string[], argv: any, output: OutputFormat): Promise<boolean> {
  if (!args[0] || !args[1]) return usage('upload <schema.table> <file.csv|.d42|.json> [--upsert] [--no-all-or-nothing] [--error-on-duplicate]');
  const a = parseDomainAddress(args[0], {table: true});
  if (!fs.existsSync(args[1])) throw new Error(`Local file not found: ${args[1]}`);
  const report = await dapi.domains.batch(a.schema, a.table!, fs.readFileSync(args[1]), contentTypeFor(args[1]), {
    mode: argv.upsert === true ? 'upsert' : 'insert',
    allOrNothing: argv['all-or-nothing'] !== false,
    errorOnDuplicate: argv['error-on-duplicate'] === true,
  });
  printBatchReport(report, output);
  return true;
}

export function printBatchReport(report: any, output: OutputFormat): void {
  const rows: any[] = Array.isArray(report?.rows) ? report.rows : [];
  const failed = (report?.errorCount ?? 0) > 0 || !!report?.error;
  if (output === 'json') printOutput(report, output);
  else if (output === 'csv') printOutput(rows, output);
  else if (output === 'quiet') {
    for (const r of rows)
      if (r?.id) console.log(r.id);
  }
  else {
    console.log(`inserted: ${report?.inserted ?? 0}  updated: ${report?.updated ?? 0}  skipped: ${report?.skipped ?? 0}  errors: ${report?.errorCount ?? 0}${report?.error ? `  (${report.error})` : ''}`);
    const errors = rows.filter((r) => r?.status === 'error' || r?.status === 'duplicate');
    if (errors.length)
      printOutput(errors.map((r) => ({
        index: r.index, status: r.status, id: r.id ?? r.existingId ?? '',
        detail: (r.errors ?? []).map((e: any) => `${e.column ? e.column + ': ' : ''}${e.message ?? e.code}`).join('; '),
      })), 'table');
  }
  if (failed) process.exitCode = 1;
}

async function handleDownload(dapi: NodeDapi, args: string[], argv: any, filter: string, limit: number | undefined,
                              offset: number, output: OutputFormat): Promise<boolean> {
  if (!args[0]) return usage('download <schema.table> [-O out.csv|out.d42] [--filter <expr>] [--columns a,b] [--sort s] [--limit n]');
  const a = parseDomainAddress(args[0], {table: true});
  const outFile: string | undefined = argv['output-file'] ?? argv.O;
  const spec = querySpec(argv, filter, limit, offset);
  if (outFile && path.extname(String(outFile)).toLowerCase() === '.d42') {
    const bytes = await dapi.domains.queryD42(a.schema, a.table!, spec);
    fs.writeFileSync(String(outFile), bytes);
    if (output !== 'quiet') console.log(`Wrote ${bytes.length} bytes to ${outFile}`);
    return true;
  }
  const csv = rowsToCsv(await dapi.domains.query(a.schema, a.table!, spec));
  if (outFile) {
    fs.writeFileSync(String(outFile), csv);
    if (output !== 'quiet') console.log(`Wrote ${csv.length} bytes to ${outFile}`);
  }
  else
    process.stdout.write(csv);
  return true;
}

async function handleGrant(dapi: NodeDapi, verb: string, args: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const [address, ...groupArgs] = args;
  if (!address || !groupArgs.length)
    return usage(`${verb} <schema> | <schema.table> <group>[,<group>...] [--access ${PERMISSIONS.join('|')}]`);
  const access = argv.access === undefined ? (verb === 'grant' ? 'View' : undefined) : String(argv.access);
  if (access !== undefined && !PERMISSIONS.includes(access))
    throw new Error(`Invalid --access '${access}'. Use one of ${PERMISSIONS.join(', ')}.`);
  const target: DomainAddress = parseDomainAddress(address);
  const entityId = await dapi.domains.entityId(target);
  const groups = groupArgs.flatMap((g) => g.split(',')).map((g) => g.trim()).filter(Boolean);
  const results: any[] = [];
  for (const g of groups) {
    try {
      const group = await dapi.groups.resolve(g);
      if (verb === 'grant') await dapi.domains.grant(entityId, group.id, access!);
      else await dapi.domains.revoke(entityId, group.id, access);
      results.push({group: g, permission: access ?? 'all', status: verb === 'grant' ? 'granted' : 'revoked'});
    } catch (err: any) {
      results.push({group: g, permission: access ?? 'all', status: 'error', error: err?.apiError?.error ?? err?.message ?? String(err)});
    }
  }
  printOutput(results, output);
  if (results.some((r) => r.status === 'error')) process.exitCode = 1;
  return true;
}
