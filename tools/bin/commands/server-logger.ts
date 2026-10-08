/// `grok s observe logger ...` — the server's logging policy: base settings, time-boxed overrides, locks, history
/// (ActionLoggerRouter, `/logging/policy`).
import {NodeDapi} from '../utils/node-dapi';
import {Connect, eachHost, forEachHost, hostList} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {applyListSpec, fmtTime, normalizeFlag, optString, parseDuration, parseTime, printBlock, rows,
  truncate, valueText} from '../utils/obs-format';

export const LOGGER_USAGE = `Usage: grok s observe logger <verb> [server] [options]
  get [server] [--scope user:<login>|group:<name>|session:<id>|package:<name>] [--host a --host b ...]
  set [server] [--print-levels L] [--post-levels L] [--save-levels L] [--debug-flags L] [--print-format text|json]
               [--print-details true|false] [--print-date-time ...] [--print-stack-traces ...] [--print-client-messages ...]
               [--print-synced-messages ...] [--additional-func-logs ...] [--set <path>=<json>]...
               [--scope all|group:<name>|user:<login>|session:<id>|package:<name>] [--for 30m | --until <iso>] [--reason <text>]
  diff [--version <n> | --host <a> --host <b>]
  overrides [--host a --host b ...]
  history [--limit 20]
  revert [<version> | --override <id> | --overrides] [--reason <text>]
L is a comma list: a,b replaces the list, +a,-b adds and removes (--debug-flags +queries, --save-levels -audit).
--for/--until make a time-boxed override; a user, session or package scope without them lasts 1 h;
scope all or group without them changes the base settings.`;

type Kind = 'levels' | 'flags' | 'format' | 'bool';

/** The overridable LoggerSettings properties, by their kebab-case flag. */
const PROPS: Record<string, {prop: string; kind: Kind}> = {
  'print-levels': {prop: 'printLevels', kind: 'levels'},
  'post-levels': {prop: 'postLevels', kind: 'levels'},
  'save-levels': {prop: 'saveLevels', kind: 'levels'},
  'debug-flags': {prop: 'debugFlags', kind: 'flags'},
  'print-format': {prop: 'printFormat', kind: 'format'},
  'print-details': {prop: 'printDetails', kind: 'bool'},
  'print-date-time': {prop: 'printDateTime', kind: 'bool'},
  'print-stack-traces': {prop: 'printStackTraces', kind: 'bool'},
  'print-client-messages': {prop: 'printClientMessages', kind: 'bool'},
  'print-synced-messages': {prop: 'printSyncedMessages', kind: 'bool'},
  'additional-func-logs': {prop: 'additionalFuncLogs', kind: 'bool'},
};

const SCOPES = ['all', 'group', 'user', 'session', 'package'];
const GROUP_PATH = /^userGroupSettings\.([^.]+)\.(.+)$/;
const EXPORT_PATH = /^exportSettings\.([^.]+)\.(.+)$/;

export async function handleLogger(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                   output: OutputFormat): Promise<boolean> {
  if ((verb === 'get' || verb === 'set') && rest[0] !== undefined && String(rest[0]) !== 'server')
    throw new Error(`Only 'server' is supported (got '${rest[0]}')`);
  switch (verb) {
    case 'get': return await loggerGet(connect, argv, output);
    case 'set': return await loggerSet(connect, argv, output);
    case 'diff': return await loggerDiff(connect, argv, output);
    case 'overrides': {
      printOutput(await forEachHost(argv, connect, async (dapi) =>
        rows(await dapi.logging.overrides() ?? [], output, overrideRow), output), output);
      return true;
    }
    case 'history': {
      const history: any[] = await (await connect()).logging.history({limit: argv.limit ?? 20}) ?? [];
      printOutput(output === 'json' ? history : history.map((h) => ({
        VERSION: h?.version ?? '', CHANGED: fmtTime(h?.changedAt), BY: h?.changedBy ?? '', SOURCE: h?.source ?? '',
        REASON: truncate(h?.reason, 40), CHANGES: truncate(h?.changed, 60),
      })), output);
      return true;
    }
    case 'revert': {
      const res = await (await connect()).logging.revert(revertBody(rest, argv));
      if (output === 'table') console.log(`reverted ${res?.reverted ?? ''}`);
      else printOutput(res, output);
      return true;
    }
  }
  printError(new Error(LOGGER_USAGE));
  return false;
}

export function revertBody(rest: string[], argv: any): Record<string, any> {
  const given = [rest[0] !== undefined, argv.override !== undefined, argv.overrides === true].filter(Boolean).length;
  if (given > 1) throw new Error('logger revert takes one of <version>, --override <id>, --overrides');
  const body: Record<string, any> = {reason: optString(argv.reason)};
  if (rest[0] !== undefined) {
    if (!/^\d+$/.test(String(rest[0]))) throw new Error(`logger revert <version>: a history version number, got '${rest[0]}'`);
    body.version = Number(rest[0]);
  }
  if (argv.override !== undefined) body.override = String(argv.override);
  if (argv.overrides === true) body.allOverrides = true;
  return body;
}

export interface Scope { type: string; value?: string }

export function parseScope(value: any): Scope {
  const s = value === undefined || value === true ? 'all' : String(value);
  const colon = s.indexOf(':');
  const type = colon < 0 ? s : s.slice(0, colon);
  const scopeValue = colon < 0 ? undefined : s.slice(colon + 1);
  if (!SCOPES.includes(type) || (type === 'all') !== (scopeValue === undefined) || scopeValue === '')
    throw new Error(`--scope is all, group:<name>, user:<login>, session:<id> or package:<name>, got '${s}'`);
  return {type, value: scopeValue};
}

/** The All Users group of a policy document: its settings are the display path `server.<prop>`. */
function allUsersId(policy: any): string | undefined {
  return Object.entries(policy?.groups ?? {}).find(([, name]) => /^all users$/i.test(String(name)))?.[0];
}

async function groupId(dapi: NodeDapi, policy: any, name: string): Promise<string> {
  const known = Object.entries(policy?.groups ?? {}).find(([, n]) => String(n).toLowerCase() === name.toLowerCase());
  return known ? known[0] : (await dapi.groups.resolve(name)).id;
}

/** The override's `scopeId` and the matching `effective` query parameter. */
async function scopeTarget(dapi: NodeDapi, policy: any, scope: Scope): Promise<{scopeId?: string; param: Record<string, string>}> {
  if (scope.type === 'all') return {param: {}};
  let id = scope.value!;
  if (scope.type === 'group') id = await groupId(dapi, policy, id);
  if (scope.type === 'user') id = (await dapi.users.find(id))?.id ?? id;
  return {scopeId: id, param: {[scope.type]: id}};
}

/** A canonical settings path as operators read it: `server.<prop>`, `server.groups.<name>.<prop>`, `server.export.<label>.<prop>`. */
export function displayPath(path: string, policy: any): string {
  const g = GROUP_PATH.exec(path);
  if (g)
    return g[1] === allUsersId(policy) ? `server.${g[2]}` : `server.groups.${policy?.groups?.[g[1]] ?? g[1]}.${g[2]}`;
  const e = EXPORT_PATH.exec(path);
  if (e)
    return `server.export.${policy?.destinations?.[e[1]] ?? e[1]}.${e[2]}`;
  return `server.${path}`;
}

function displayMap(flat: Record<string, any> | undefined, policy: any): Record<string, any> {
  const out: Record<string, any> = {};
  for (const [path, v] of Object.entries(flat ?? {}))
    out[displayPath(path, policy)] = v;
  return out;
}

export interface DiffRow { change: string; path: string; value: string; scope?: string; reverts?: string }

/** `+` only on the right, `-` only on the left, `~` on both with different values. */
export function diffMaps(left: Record<string, any>, right: Record<string, any>): DiffRow[] {
  const rows: DiffRow[] = [];
  for (const path of [...new Set([...Object.keys(left), ...Object.keys(right)])].sort()) {
    const inLeft = path in left;
    const inRight = path in right;
    if (inLeft && inRight && JSON.stringify(left[path]) === JSON.stringify(right[path])) continue;
    if (!inLeft) rows.push({change: '+', path, value: valueText(right[path])});
    else if (!inRight) rows.push({change: '-', path, value: valueText(left[path])});
    else rows.push({change: '~', path, value: `${valueText(left[path])} → ${valueText(right[path])}`});
  }
  return rows;
}

function overrideScope(o: any): string {
  return o?.scope === 'all' ? 'all' : `${o?.scope}:${o?.scopeName ?? o?.scopeId ?? ''}`;
}

function overrideRow(o: any): Record<string, any> {
  return {
    ID: String(o?.id ?? '').slice(0, 8),
    SCOPE: overrideScope(o),
    CHANGES: Object.entries(o?.changes ?? {}).map(([k, v]) => `${k}=${valueText(v)}`).join('; '),
    REVERTS: fmtTime(o?.expiresAt),
    BY: o?.createdByLogin ?? o?.createdBy ?? '',
    REASON: truncate(o?.reason, 40),
  };
}

function printDiff(rows: DiffRow[], output: OutputFormat): void {
  if (output !== 'table') { printOutput(rows, output); return; }
  if (!rows.length) { console.log('(no differences)'); return; }
  const pw = Math.max(...rows.map((r) => r.path.length));
  const vw = Math.max(...rows.map((r) => r.value.length));
  for (const r of rows) {
    const extra = [r.scope ? `scope ${r.scope}` : '', r.reverts ? `reverts ${r.reverts}` : ''].filter(Boolean).join('   ');
    console.log(`${r.change} ${r.path.padEnd(pw)}  ${r.value.padEnd(vw)}   ${extra}`.trimEnd());
  }
}

const SHORT_NAMES: Record<string, string> = {printLevels: 'print', postLevels: 'post', saveLevels: 'save', debugFlags: 'debug flags'};

const levels = (v: any) => Array.isArray(v) ? (v.length ? v.join(' ') : '(none)') : valueText(v);

function printPolicy(policy: any): void {
  const all = allUsersId(policy);
  const value = (gid: string | undefined, prop: string) => {
    const path = `userGroupSettings.${gid}.${prop}`;
    return policy?.settings?.[path] ?? policy?.defaults?.[path];
  };
  const lines: [string, string][] = [
    ['print', `${levels(value(all, 'printLevels')).padEnd(28)}post  ${levels(value(all, 'postLevels'))}`],
    ['save', `${levels(value(all, 'saveLevels')).padEnd(28)}debug flags  ${levels(value(all, 'debugFlags'))}`],
  ];
  const locked: string[] = Array.isArray(policy?.locked) ? policy.locked : [];
  if (locked.length)
    lines.push(['locked', `${locked.join(', ')}   # from deployment configuration`]);
  const others = new Map<string, string[]>();
  for (const [path, v] of Object.entries(policy?.settings ?? {})) {
    const g = GROUP_PATH.exec(path);
    if (g && g[1] !== all)
      others.set(g[1], [...others.get(g[1]) ?? [], `${SHORT_NAMES[g[2]] ?? g[2]} ${levels(v)}`]);
  }
  for (const [gid, props] of others)
    lines.push([`group ${policy?.groups?.[gid] ?? gid}`, props.join('   ')]);
  for (const o of Array.isArray(policy?.overrides) ? policy.overrides : [])
    lines.push(['override', `${overrideScope(o)}  ${Object.entries(o?.changes ?? {}).map(([k, v]) => `${k} ${valueText(v)}`).join(', ')}` +
      `  reverts ${fmtTime(o?.expiresAt)}  by ${o?.createdByLogin ?? o?.createdBy ?? ''}  ${o?.reason ?? ''}`]);
  printBlock(lines);
}

async function loggerGet(connect: Connect, argv: any, output: OutputFormat): Promise<boolean> {
  const scope = argv.scope === undefined ? undefined : parseScope(argv.scope);
  if (scope?.type === 'all')
    throw new Error('logger get --scope takes user:, group:, session: or package:; without it the base policy is shown');
  const docs: any[] = [];
  await eachHost(argv, connect, async (dapi, host, multi) => {
    const logging = dapi.logging;
    const policy = await logging.policy();
    const doc = scope ? await logging.effective((await scopeTarget(dapi, policy, scope)).param) : policy;
    if (output !== 'table') { docs.push(multi ? {host, ...doc} : doc); return; }
    if (multi) console.log(`== ${host}`);
    if (!scope) printPolicy(doc);
    else printOutput(Object.entries(doc?.settings ?? {})
      .map(([prop, v]) => ({SETTING: prop, VALUE: valueText(v), SOURCE: doc?.sources?.[prop] ?? ''})), 'table');
  });
  if (output !== 'table') printOutput(docs.length === 1 ? docs[0] : docs, output);
  return true;
}

function boolArg(value: any, flag: string): boolean {
  if (value === true || value === 'true') return true;
  if (value === false || value === 'false') return false;
  throw new Error(`--${flag} is true or false, got '${value}'`);
}

export function propChanges(argv: any, current: (prop: string) => any): Record<string, any> {
  const changes: Record<string, any> = {};
  for (const [flag, {prop, kind}] of Object.entries(PROPS)) {
    const v = argv[flag];
    if (v === undefined) continue;
    if (kind === 'levels') changes[prop] = applyListSpec(current(prop), v, (s) => s.toLowerCase(), `--${flag}`);
    else if (kind === 'flags') changes[prop] = applyListSpec(current(prop), v, normalizeFlag, `--${flag}`);
    else if (kind === 'format') {
      if (v !== 'text' && v !== 'json') throw new Error(`--${flag} is text or json, got '${v}'`);
      changes[prop] = v;
    }
    else changes[prop] = boolArg(v, flag);
  }
  return changes;
}

export function setArgs(value: any): Record<string, any> {
  const out: Record<string, any> = {};
  for (const item of value === undefined ? [] : Array.isArray(value) ? value : [value]) {
    const s = String(item);
    const eq = s.indexOf('=');
    if (eq <= 0) throw new Error(`--set expects <path>=<json>, got '${s}'`);
    let v: any = s.slice(eq + 1);
    try { v = JSON.parse(v); } catch { /* a bare string */ }
    out[s.slice(0, eq)] = v;
  }
  return out;
}

export function isOverride(argv: any, scope: Scope): boolean {
  return argv.for !== undefined || argv.until !== undefined || ['user', 'session', 'package'].includes(scope.type);
}

async function loggerSet(connect: Connect, argv: any, output: OutputFormat): Promise<boolean> {
  const scope = parseScope(argv.scope);
  const override = isOverride(argv, scope);
  const paths = setArgs(argv.set);
  if (override && Object.keys(paths).length)
    throw new Error('--set changes the base settings; it cannot be combined with --for, --until or a user, session or package scope');
  if (argv.for !== undefined && argv.until !== undefined)
    throw new Error('Use either --for or --until');
  if (!Object.keys(PROPS).some((f) => argv[f] !== undefined) && !Object.keys(paths).length) {
    printError(new Error(`Nothing to set.\n${LOGGER_USAGE}`));
    return false;
  }
  const reason = optString(argv.reason);
  const dapi = await connect();
  const logging = dapi.logging;
  const policy = await logging.policy();
  if (override) {
    const target = await scopeTarget(dapi, policy, scope);
    const effective = (await logging.effective(target.param))?.settings ?? {};
    const body: Record<string, any> = {scope: scope.type, scopeId: target.scopeId, set: propChanges(argv, (p) => effective[p]), reason};
    if (argv.for !== undefined) body.forMinutes = parseDuration(argv.for, '--for') / 60000;
    if (argv.until !== undefined) body.expiresAt = parseTime(argv.until, '--until').toISOString();
    const created = await logging.addOverride(body);
    if (output !== 'table') { printOutput(created, output); return true; }
    const label = scope.type === 'all' ? 'all' : `${scope.type}:${scope.value}`;
    for (const [prop, v] of Object.entries(body.set))
      console.log(`+ server.${prop}  ${valueText(v)}  scope ${label}  reverts ${fmtTime(created?.expiresAt)}`);
    return true;
  }
  const gid = scope.type === 'all'
    ? allUsersId(policy) ?? (await dapi.groups.resolve('All users')).id
    : await groupId(dapi, policy, scope.value!);
  const all = allUsersId(policy);
  const current = (prop: string) =>
    policy?.settings?.[`userGroupSettings.${gid}.${prop}`] ?? policy?.settings?.[`userGroupSettings.${all}.${prop}`];
  const set: Record<string, any> = {...paths};
  for (const [prop, v] of Object.entries(propChanges(argv, current)))
    set[`userGroupSettings.${gid}.${prop}`] = v;
  const res = await logging.setPolicy({set, reason});
  if (output !== 'table') { printOutput(res, output); return true; }
  if (res?.version == null) { console.log('(no change: the base settings already have these values)'); return true; }
  const named = {...policy, groups: {...policy?.groups, [gid]: policy?.groups?.[gid] ?? scope.value ?? 'All users'}};
  for (const [path, v] of Object.entries(set))
    console.log(`~ ${displayPath(path, named)}  ${valueText(v)}`);
  console.log(`policy version ${res.version}`);
  return true;
}

async function loggerDiff(connect: Connect, argv: any, output: OutputFormat): Promise<boolean> {
  const hosts = hostList(argv.host);
  if (hosts.length > 1) {
    if (hosts.length !== 2 || argv.version !== undefined) throw new Error('logger diff compares exactly two --host, without --version');
    const maps: Record<string, any>[] = [];
    for (const host of hosts) {
      const policy = await (await connect(host)).logging.policy();
      maps.push(displayMap(policy?.settings, policy));
    }
    printDiff(diffMaps(maps[0], maps[1]), output);
    return true;
  }
  const logging = (await connect()).logging;
  if (argv.version !== undefined) {
    if (!/^\d+$/.test(String(argv.version))) throw new Error(`--version is a history version number, got '${argv.version}'`);
    const [then, now] = [await logging.policy({version: String(argv.version)}), await logging.policy()];
    printDiff(diffMaps(displayMap(then?.settings, then), displayMap(now?.settings, now)), output);
    return true;
  }
  const policy = await logging.policy({defaults: true});
  const diff = diffMaps(displayMap(policy?.defaults, policy), displayMap(policy?.settings, policy));
  for (const o of Array.isArray(policy?.overrides) ? policy.overrides : [])
    for (const [prop, v] of Object.entries(o?.changes ?? {}))
      diff.push({change: '+', path: `server.${prop}`, value: valueText(v), scope: overrideScope(o), reverts: fmtTime(o?.expiresAt)});
  printDiff(diff, output);
  return true;
}
