/// `grok s errors ...` — platform errors as query results (ErrorsRouter, `/errors`).
import * as fs from 'fs';
import {Query} from '../utils/node-observability';
import {Connect, forEachHost, hostList, singleHost} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {cronFromSchedule, fmtDate, fmtMinutes, fmtTime, hasValue, listArg, parseTime, printBlock, shortRequestId,
  shortSig, sinceArg, slug, sparkline, truncate} from '../utils/obs-format';
import {rowsToCsv} from './server-domains';

export const ERRORS_USAGE = `Usage: grok s errors <verb> [filters] [options]
  list [filters] [--limit 50] [--host a --host b ...]
  top [filters] [--by <dim>[,<dim>[,<dim>]]] [--trend hour|day] [--limit 20] [--format csv|json|parquet] [-O file] [--host ...]
  show <signature> [--since 24h]
  diff --before <from>..<to> --after <from>..<to> [filters]
  diff --since 7d --host <a> --host <b> [filters]
  export [filters] [--by ...] --format csv|json|parquet [-O file]
  save "<name>" [filters] [--by ...] --to "<Share:path/>" [--format csv|json] [--schedule "MON 07:00" | --schedule "<cron>"]
Filters: --since 7d | --from <iso|-7d> --to <iso|-1d>, --signature s, --package p, --version v, --user login,
  --group name, --service server|client, --route "<METHOD> /api/...", --server name, --connection <Ns:Name>,
  --function <nqName>, --regressed, --min-users n, --min-count n
Dimensions: signature package version user group service route server connection function`;

export const DIMENSIONS = ['signature', 'package', 'version', 'user', 'group', 'service', 'route', 'server', 'connection', 'function'];
const FORMATS = ['csv', 'json', 'parquet'];
const STRING_FILTERS = ['signature', 'package', 'version', 'user', 'group', 'server', 'connection', 'function'];

export async function handleErrors(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                   output: OutputFormat): Promise<boolean> {
  switch (verb) {
    case 'list': {
      const q = {...errorFilters(argv), limit: argv.limit ?? 50, offset: argv.offset};
      const rows = await forEachHost(argv, connect, async (dapi) => {
        const occurrences: any[] = await dapi.errors.query(q) ?? [];
        return output === 'json' ? occurrences : occurrences.map(occurrenceRow);
      }, output);
      printOutput(rows, output);
      return true;
    }
    case 'top': {
      const by = byArg(argv.by, 'signature');
      const trend = trendArg(argv.trend);
      const q = {...errorFilters(argv), by: by.join(','), trend, limit: argv.limit ?? 20};
      if (argv.format !== undefined)
        return await writeExport(connect, argv, q, output);
      const rows = await forEachHost(argv, connect, async (dapi) => {
        const aggregates: any[] = await dapi.errors.query(q) ?? [];
        return output === 'json' ? aggregates : aggregates.map((r) => aggregateRow(r, by, trend));
      }, output);
      printOutput(rows, output);
      return true;
    }
    case 'export': {
      if (argv.format === undefined) return usage('export [filters] [--by ...] --format csv|json|parquet [-O file]');
      const by = byArg(argv.by);
      const q = {...errorFilters(argv), by: by.length ? by.join(',') : undefined,
        trend: by.length ? trendArg(argv.trend) : undefined, limit: argv.limit};
      return await writeExport(connect, argv, q, output);
    }
    case 'show': {
      if (rest[0] === undefined) return usage('show <signature> [--since 24h]');
      const since = sinceArg(argv.since ?? '24h');
      const dapi = await connect(singleHost(argv, 'errors show'));
      const doc = await dapi.errors.show(String(rest[0]), {since});
      if (output === 'table') printShow(doc);
      else printOutput(doc, output);
      return true;
    }
    case 'diff':
      return hostList(argv.host).length > 1 ? await diffHosts(connect, argv, output) : await diffWindows(connect, argv, output);
    case 'save': return await saveJob(connect, rest, argv, output);
  }
  printError(new Error(ERRORS_USAGE));
  return false;
}

function usage(line: string): boolean {
  printError(new Error(`Usage: grok s errors ${line}`));
  return false;
}

/** `"POST /api/queries/{id}/run"` → `"POST /queries/{id}/run"`: routes are recorded without the `/api` prefix. */
export function normalizeRoute(route: string): string {
  return String(route).trim().replace(/^((?:[A-Za-z]+\s+)?)\/api(?=\/|$)/, '$1').replace(/^([A-Za-z]+\s+)$/, '$1/') || '/';
}

/**
 * The shared filter flags as query parameters. [time] false leaves out the window (`diff` has its
 * own); [range] false refuses `--from`/`--to` (`save`, where `--to` is the destination).
 */
export function errorFilters(argv: any, opts: {time?: boolean; range?: boolean; defaultSince?: string} = {},
                             now: Date = new Date()): Query {
  const q: Query = {};
  if (opts.time !== false) {
    if (opts.range !== false && (argv.from !== undefined || argv.to !== undefined)) {
      if (argv.since !== undefined) throw new Error('Use either --since or --from/--to');
      if (argv.from !== undefined) q.from = parseTime(argv.from, '--from', now).toISOString();
      if (argv.to !== undefined) q.to = parseTime(argv.to, '--to', now).toISOString();
    }
    else {
      if (opts.range === false && argv.from !== undefined)
        throw new Error('A saved job takes --since (--to names the destination)');
      q.since = sinceArg(argv.since ?? opts.defaultSince ?? '24h');
    }
  }
  for (const f of STRING_FILTERS)
    if (hasValue(argv[f])) q[f] = String(argv[f]);
  if (argv.service !== undefined) {
    if (!['server', 'client'].includes(String(argv.service))) throw new Error(`--service is server or client, got '${argv.service}'`);
    q.service = String(argv.service);
  }
  if (argv.route !== undefined) q.route = normalizeRoute(argv.route);
  if (argv.regressed === true) q.regressed = true;
  if (argv['min-users'] !== undefined) q.minUsers = Number(argv['min-users']);
  if (argv['min-count'] !== undefined) q.minCount = Number(argv['min-count']);
  return q;
}

export function byArg(value: any, fallback?: string): string[] {
  const by = listArg(value);
  if (!by.length) return fallback ? [fallback] : [];
  for (const d of by)
    if (!DIMENSIONS.includes(d)) throw new Error(`Unknown --by dimension '${d}'. Valid: ${DIMENSIONS.join(', ')}`);
  if (by.length > 3) throw new Error('--by takes at most three dimensions');
  return by;
}

function trendArg(value: any): string {
  const trend = String(value ?? 'day');
  if (trend !== 'hour' && trend !== 'day') throw new Error(`--trend is hour or day, got '${value}'`);
  return trend;
}

export function occurrenceRow(r: any): Record<string, any> {
  return {
    TIME: fmtTime(r?.time),
    USER: r?.user ?? '',
    SOURCE: r?.service ?? '',
    SIG: shortSig(r?.signature),
    ERROR: truncate(r?.error, 50),
    PACKAGE: r?.package ?? '',
    VERSION: r?.version ?? '',
    ROUTE: r?.route ?? '',
    SERVER: r?.server ?? '',
    REQ: shortRequestId(r?.requestId),
  };
}

function stateText(r: any): string {
  if (r?.state === 'muted')
    return r.stateVersion ? `muted → ${r.stateVersion}` : r.stateUntil ? `muted until ${fmtTime(r.stateUntil)}` : 'muted';
  return r?.state ?? '';
}

/** One aggregate row: signature-level figures when grouped by signature, rollups (SIGNATURES, NEW) otherwise. */
export function aggregateRow(r: any, by: string[], trend: string): Record<string, any> {
  const row: Record<string, any> = {};
  const bySignature = by.includes('signature');
  for (const d of bySignature ? ['signature', ...by.filter((x) => x !== 'signature')] : by)
    row[d === 'signature' ? 'SIG' : d.toUpperCase()] = d === 'signature' ? shortSig(r?.signature) : (r?.[d] ?? '');
  const trendCol = trend === 'hour' ? 'TREND (hourly)' : 'TREND (daily)';
  if (bySignature) {
    row.ERROR = truncate(r?.topError, 50);
    row.USERS = r?.users ?? 0;
    row.COUNT = r?.count ?? 0;
    row['FIRST SEEN'] = [r?.firstVersion, fmtDate(r?.firstSeen)].filter(Boolean).join(' · ');
    row.LAST = fmtTime(r?.lastSeen);
    row[trendCol] = sparkline(r?.trend);
    row.STATE = stateText(r);
  }
  else {
    row.SIGNATURES = r?.signatures ?? 0;
    row.USERS = r?.users ?? 0;
    row.COUNT = r?.count ?? 0;
    row.NEW = r?.newInRange ?? 0;
    row['FIRST SEEN IN'] = r?.firstVersion ?? '';
    row.LAST = fmtTime(r?.lastSeen);
    row[trendCol] = sparkline(r?.trend);
    row['TOP ERROR'] = truncate(r?.topError, 50);
  }
  return row;
}

async function writeExport(connect: Connect, argv: any, q: Query, output: OutputFormat): Promise<boolean> {
  const format = String(argv.format);
  if (!FORMATS.includes(format)) throw new Error(`--format is csv, json or parquet, got '${format}'`);
  const errors = (await connect(singleHost(argv, 'errors --format'))).errors;
  let bytes: Buffer;
  if (format === 'csv') {
    const csv = await errors.query({...q, format: 'csv'});
    bytes = Buffer.from(typeof csv === 'string' ? csv : rowsToCsv(csv ?? []));
  }
  else {
    const rows: any[] = await errors.query({...q, format: 'json'}) ?? [];
    bytes = format === 'json' ? Buffer.from(JSON.stringify(rows, null, 2) + '\n') : toParquet(rows);
  }
  const outFile: string | undefined = argv['output-file'] ?? argv.O;
  if (outFile) {
    fs.writeFileSync(String(outFile), bytes);
    if (output !== 'quiet') console.log(`Wrote ${bytes.length} bytes to ${outFile}`);
  }
  else
    process.stdout.write(bytes);
  return true;
}

/** Parquet is written here from the JSON rows; the server exports CSV and JSON only. */
export function toParquet(rows: any[], load: (name: string) => any = require): Buffer {
  let arrow: any;
  let parquet: any;
  try {
    arrow = load('apache-arrow');
    parquet = load('parquet-wasm');
  }
  catch {
    throw new Error('Parquet output needs the apache-arrow and parquet-wasm packages, which this grok does not have ' +
      'installed; use --format csv or --format json');
  }
  const flat = rows.map((r) => Object.fromEntries(Object.entries(r ?? {})
    .map(([k, v]) => [k, v !== null && typeof v === 'object' ? JSON.stringify(v) : v])));
  const ipc = arrow.tableToIPC(arrow.tableFromJSON(flat), 'stream');
  return Buffer.from(parquet.writeParquet(parquet.Table.fromIPCStream(ipc)));
}

/** `top: <sig> <package> <error> · <n> users`, or `· <after> vs <before>` for a risen signature. */
function topLine(t: any, risen: boolean = false): string {
  if (!t) return '';
  const what = [shortSig(t.signature), t.package, risen ? '' : truncate(t.error, 50)].filter(Boolean).join(' ');
  return `top: ${what}${risen ? ` · ${t.after ?? 0} vs ${t.before ?? 0}` : t.users !== undefined ? ` · ${t.users} users` : ''}`;
}

async function diffWindows(connect: Connect, argv: any, output: OutputFormat): Promise<boolean> {
  if (argv.before === undefined || argv.after === undefined)
    return usage('diff --before <from>..<to> --after <from>..<to> [filters]   or   diff --since 7d --host <a> --host <b>');
  for (const w of [argv.before, argv.after])
    if (!String(w).includes('..')) throw new Error(`A window is <from>..<to>, got '${w}'`);
  const q = {...errorFilters(argv, {time: false}), before: String(argv.before), after: String(argv.after)};
  const dapi = await connect(singleHost(argv, 'errors diff --before/--after'));
  const doc = await dapi.errors.diff(q);
  if (output === 'json') { printOutput(doc, output); return true; }
  if (output === 'csv' || output === 'quiet') { printOutput(doc?.rows ?? [], output); return true; }
  const c = doc?.categories ?? {};
  const lines: [string, any, string, string][] = [
    ['NEW', c.new, `first seen in ${argv.after}`, topLine(c.new?.top)],
    ['GONE', c.gone, `not seen in ${argv.after}`, topLine(c.gone?.top)],
    ['RISEN', c.risen, '≥ 2× occurrences', topLine(c.risen?.top, true)],
    ['REGRESSED', c.regressed, 'muted, back on a newer version', topLine(c.regressed?.top)],
  ];
  const width = Math.max(...lines.map((l) => l[2].length)) + 3;
  for (const [label, cat, text, top] of lines)
    console.log(`${label.padEnd(10)}${String(cat?.count ?? 0).padStart(4)}   ${text.padEnd(width)}${top}`.trimEnd());
  const i = doc?.incidents ?? {};
  console.log(`${'incidents'.padEnd(10)}${i.opened ?? 0} opened · ${i.resolved ?? 0} resolved · ` +
    `median time to resolve ${fmtMinutes(i.medianTtrMinutes)}`);
  const r = doc?.reports ?? {};
  console.log(`${'reports'.padEnd(10)}${r.total ?? 0} (human ${r.human ?? 0}, auto ${r.auto ?? 0}) · ` +
    `${r.linkedToNew ?? 0} linked to a NEW signature`);
  return true;
}

/** Signatures of two deployments over one window, by where they occur. */
export function classifyHosts(a: any[], b: any[]): {onlyA: any[]; onlyB: any[]; both: string[]} {
  const bySig = (rows: any[]) => {
    const m = new Map<string, any>();
    for (const r of rows) {
      const prev = m.get(r?.signature);
      if (!prev || (r?.users ?? 0) > (prev.users ?? 0)) m.set(r?.signature, r);
    }
    return m;
  };
  const ma = bySig(a);
  const mb = bySig(b);
  const rank = (x: any, y: any) => (y.users ?? 0) - (x.users ?? 0) || (y.count ?? 0) - (x.count ?? 0);
  return {
    onlyA: [...ma.values()].filter((r) => !mb.has(r.signature)).sort(rank),
    onlyB: [...mb.values()].filter((r) => !ma.has(r.signature)).sort(rank),
    both: [...ma.keys()].filter((s) => mb.has(s)),
  };
}

async function diffHosts(connect: Connect, argv: any, output: OutputFormat): Promise<boolean> {
  const hosts = hostList(argv.host);
  if (hosts.length !== 2) throw new Error('errors diff compares exactly two --host');
  const q = {...errorFilters(argv), by: 'signature,package', limit: 10000};
  const sides = [];
  for (const host of hosts) {
    const dapi = await connect(host);
    const version = (await dapi.serverInfo()).version;
    sides.push({host, version, rows: (await dapi.errors.query(q) ?? []) as any[]});
  }
  const {onlyA, onlyB, both} = classifyHosts(sides[0].rows, sides[1].rows);
  if (output === 'json') {
    printOutput({hosts: [{host: sides[0].host, version: sides[0].version, only: onlyA},
      {host: sides[1].host, version: sides[1].version, only: onlyB}], both}, output);
    return true;
  }
  const top = (r: any) => topLine(r && {signature: r.signature, package: r.package, error: r.topError, users: r.users ?? 0});
  const labels = sides.map((s) => `ONLY ON ${s.host} (${s.version})`);
  const width = Math.max(...labels.map((l) => l.length), 'ON BOTH'.length) + 2;
  console.log(`${labels[0].padEnd(width)}${String(onlyA.length).padStart(4)}   ${top(onlyA[0])}`.trimEnd());
  console.log(`${labels[1].padEnd(width)}${String(onlyB.length).padStart(4)}   ${top(onlyB[0])}`.trimEnd());
  console.log(`${'ON BOTH'.padEnd(width)}${String(both.length).padStart(4)}`);
  return true;
}

async function saveJob(connect: Connect, rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const name = rest[0] === undefined ? undefined : String(rest[0]);
  if (!name || argv.to === undefined || argv.to === true)
    return usage('save "<name>" [filters] [--by ...] --to "<Share:path/>" [--format csv|json] [--schedule "MON 07:00"]');
  const format = String(argv.format ?? 'csv');
  if (format === 'parquet')
    throw new Error('Saved jobs export csv or json; Parquet is written by the CLI only (errors export --format parquet)');
  if (format !== 'csv' && format !== 'json') throw new Error(`--format is csv or json, got '${format}'`);
  const by = byArg(argv.by);
  const spec: Query = {...errorFilters(argv, {range: false}), by: by.length ? by.join(',') : undefined};
  const to = String(argv.to);
  const path = to.endsWith('/') ? `${to}${slug(name)}-{date}.${format}` : to;
  const cron = argv.schedule === undefined ? undefined : cronFromSchedule(argv.schedule);
  const dapi = await connect(singleHost(argv, 'errors save'));
  const job = await dapi.errors.saveJob({name, spec, format, path, cron});
  if (output === 'json' || output === 'csv') printOutput(job, output);
  else if (output === 'quiet') console.log(job?.id ?? '');
  else console.log(`saved job "${name}" ${cron ?? '(no schedule)'} → ${path}`);
  return true;
}

function printShow(d: any): void {
  const groups: any[] = Array.isArray(d?.groups) ? d.groups.slice(0, 3) : [];
  const reports: any[] = Array.isArray(d?.reports) ? d.reports : [];
  const groupText = groups.map((g) => `${g?.name} ${g?.users}`).join(' · ') || '(none)';
  const alertText = d?.alert ? `alert ${d.alert.kind}, ${d.alert.status} since ${fmtTime(d.alert.openedAt)}` : '';
  const lines: [string, string][] = [
    ['signature', `${shortSig(d?.signature).padEnd(8)}${d?.error ?? ''}`],
    ['package', `${String(d?.package ?? '').padEnd(8)}first seen in ${d?.firstVersion ?? '?'} at ${fmtTime(d?.firstSeen)}` +
      ` · last seen ${fmtTime(d?.lastSeen)}`],
    ['occurrences', `${String(d?.occurrences ?? 0).padEnd(8)}users ${String(d?.users ?? 0).padEnd(7)}groups ${groupText}`],
    ['reports', `${(reports.map((n) => `#${n}`).join(' ') || '(none)').padEnd(21)}${alertText}`],
  ];
  if (d?.change)
    lines.push(['change', `${d.change.type} ${d.change.package} ${d.change.version} by ${d.change.by} at ${fmtTime(d.change.at)}`]);
  if (d?.sessions !== undefined)
    lines.push(['sessions', String(Array.isArray(d.sessions) ? d.sessions.length : d.sessions)]);
  printBlock(lines);
}
