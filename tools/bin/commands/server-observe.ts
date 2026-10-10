/// `grok s observe` (alias `grok s o`): problems, problem rules, logger settings, errors, the timeline and cloud logs,
/// over `/problems`, `/admin/plugins/<key>/settings`, `/alerts/Test_rules`, `/log`, `/admin/metrics`, `/log/timeline`,
/// `/log/cloud/*` and `/log/archive/*`.
import * as fs from 'fs';
import {NodeDapi, buildQuery} from '../utils/node-dapi';
import {printOutput, OutputFormat} from '../utils/server-output';

export const OBSERVE_USAGE = `Usage: grok s observe <command> <verb> [args]      (alias: grok s o)

  problems list [--status <s>] [--state <s>] [--alert-status <s>] [--kind <k>] [--key <k>] [--since 7d] [--limit 100]
                                       Comma lists or all; --alert-status open lists the open alerts
  problems get <id|kind:key>
  problems history <id|kind:key> [--limit 50]
  problems status <id|kind:key> --data '{"status": "muted", "reason": "...", "until": "<ISO>"}' | --json body.json
  problems ack <id|kind:key> --reason <text>
  problems resolve <id|kind:key> --reason <text>
  rules get                            The problem rules (Settings > Alerts), as JSON
  rules put --json rules.json          Replace them; the server checks them first
  rules test                           Evaluate the saved rules once, now
  logger get [<path>]                  The logger settings, or one value (exportSettings.0.types)
  logger set --set <path>=<json> ...   Change values; with --user <login> or --group <name|id>, that
             [--for 30m | --until <ISO>] [--reason <text>]   entry of userGroupSettings, removed when it expires
  logger history [--limit 20]          The log-settings-changed audit records, with each change's diff
  errors top [--since 24h | --date "this week" | --from <ISO> [--to <ISO>]] [--by message|signature]
             [--package <name>] [--user <login>] [--source <s>] [--limit 10]
  timeline --session <id> [--from <ISO|-15m>] [--to <ISO>] [--limit 500]
                                       A session's events and requests in a window: no --from, the last
                                       15 min; --to defaults to --from + 10 min, at most 2 h
  logs cloud [--group <g>] [--since 1h | --from <ISO> [--to <ISO>]] [--filter <pattern>] [--limit 1000]
                                       CloudWatch events; no --group, the instance's own group;
                                       --filter is a CloudWatch filter pattern
  logs archive list [--prefix <p>] [--since 7d] [--limit n]
                                       Objects in the log archive (WORM bucket), by key
  logs archive read <key>              One archived object decoded to events
                                       logs take [--connection <id>]; without it, the instance's AWS role

Durations are 30m (minutes), 24h, 7d, 2w; times are ISO 8601 (UTC without an offset) or
relative, written with = (--from=-2h).`;

const UNIT_MS: Record<string, number> = {m: 60000, h: 3600000, d: 86400000, w: 604800000};

/** `30m`, `2h`, `7d`, `1w` in milliseconds; `m` is minutes, unlike `grok s pull --since`, where it is months. */
export function parseDuration(value: any, flag: string): number {
  const m = /^-?(\d+)([mhdw])$/.exec(String(value ?? '').trim());
  if (!m || Number(m[1]) <= 0)
    throw new Error(`${flag} expects a duration such as 30m, 2h, 7d or 1w, got '${value}'`);
  return Number(m[1]) * UNIT_MS[m[2]];
}

/** ISO 8601 (UTC without an offset), or relative to [now] (`-2h`). */
export function parseTime(value: any, flag: string, now: Date = new Date()): Date {
  const s = String(value ?? '').trim();
  if (s.startsWith('-'))
    return new Date(now.getTime() - parseDuration(s, flag));
  const d = new Date(/^\d{4}-\d{2}-\d{2}$/.test(s) || /(Z|[+-]\d{2}:?\d{2})$/i.test(s) ? s : `${s}Z`);
  if (!/^\d{4}-\d{2}-\d{2}/.test(s) || isNaN(d.getTime()))
    throw new Error(`${flag} expects an ISO time (2026-10-04T06:00) or a relative one (-2h), got '${value}'`);
  return d;
}

/** `a.b.0.c=<json>` applied to [target]; a value that is not JSON is a string. */
export function setPath(target: any, assignment: string): string {
  const at = assignment.indexOf('=');
  if (at <= 0)
    throw new Error(`--set expects <path>=<json>, got '${assignment}'`);
  const keys = assignment.slice(0, at).split('.');
  const text = assignment.slice(at + 1);
  let value: any;
  try {
    value = JSON.parse(text);
  }
  catch {
    value = text;
  }
  let node = target;
  for (const key of keys.slice(0, -1))
    node = node[key] ??= {};
  node[keys[keys.length - 1]] = value;
  return keys[0];
}

function readBody(argv: any): any {
  if (argv.json)
    return JSON.parse(fs.readFileSync(argv.json, 'utf8'));
  if (argv.data !== undefined)
    return JSON.parse(String(argv.data));
  throw new Error('Expected a body: --data \'<json>\' or --json <file>');
}

const iso = (value: any, flag: string) => value === undefined ? undefined : parseTime(value, flag).toISOString();

export async function handleObserve(dapi: NodeDapi, command: string | undefined, verb: string | undefined,
                                    rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  const id = rest[0] === undefined ? undefined : encodeURIComponent(String(rest[0]));
  const needId = () => {
    if (!id)
      throw new Error(`Usage: grok s observe ${command} ${verb} <id|kind:key>`);
  };
  const settings = async (key: string) => (await dapi.raw('GET', `/admin/plugins/${key}/settings`))?.settings ?? {};
  const saveSettings = (key: string, body: any) => dapi.raw('POST', `/admin/plugins/${key}/settings`, body);

  switch (`${command} ${verb}`) {
    case 'problems list':
      printOutput(await dapi.problems.list({status: argv.status, state: argv.state, alertStatus: argv['alert-status'],
        kind: argv.kind, key: argv.key, since: argv.since, limit: argv.limit}), output);
      return true;
    case 'problems get': {
      needId();
      const problem = await dapi.problems.find(String(rest[0]));
      if (!problem)
        throw new Error(`No problem '${rest[0]}'`);
      printOutput(problem, output);
      return true;
    }
    case 'problems history':
      needId();
      printOutput(await dapi.raw('GET', `/problems/${id}/history${buildQuery({limit: argv.limit})}`), output);
      return true;
    case 'problems status':
      needId();
      printOutput(await dapi.raw('POST', `/problems/${id}/status`, readBody(argv)), output);
      return true;
    case 'problems ack':
    case 'problems resolve':
      needId();
      printOutput(await dapi.raw('POST', `/problems/${id}/${verb}`, {reason: argv.reason}), output);
      return true;
    case 'rules get': {
      const rules = (await settings('alerts'))?.problemRules;
      printOutput(rules ? JSON.parse(rules) : [], output === 'table' ? 'json' : output);
      return true;
    }
    case 'rules put': {
      const rules = readBody(argv);
      await saveSettings('alerts', {problemRules: typeof rules === 'string' ? rules : JSON.stringify(rules)});
      if (output !== 'quiet')
        console.log('Saved the problem rules');
      return true;
    }
    case 'rules test':
      console.log(await dapi.raw('POST', '/alerts/Test_rules'));
      return true;
    case 'logger get': {
      let value = await settings('logger');
      for (const key of rest[0] ? String(rest[0]).split('.') : [])
        value = value?.[key];
      printOutput(value, output === 'table' ? 'json' : output);
      return true;
    }
    case 'logger set':
      return loggerSet(dapi, argv, await settings('logger'), (body) => saveSettings('logger', body), output);
    case 'logger history': {
      const text = encodeURIComponent('eventType.name = "log-settings-changed"');
      const events: any[] = await dapi.raw('GET', `/log?text=${text}&limit=${argv.limit ?? 20}&page=1`) ?? [];
      printOutput(events.map((e) => {
        const params: Record<string, any> = {};
        for (const p of e.parameters ?? [])
          params[p.parameter?.name] = p.value;
        const row = {time: e.eventTime, user: params.user, changed: params.changed};
        return output === 'table' ? row : {...row, diff: params.diff ? JSON.parse(params.diff) : null};
      }), output);
      return true;
    }
    case 'errors top': {
      if (argv.since !== undefined && (argv.from !== undefined || argv.date !== undefined))
        throw new Error('--since is not used with --from or --date');
      const metrics = await dapi.raw('GET', `/admin/metrics${buildQuery({
        date_start: argv.since !== undefined ? new Date(Date.now() - parseDuration(argv.since, '--since')).toISOString()
          : iso(argv.from, '--from'),
        date_end: iso(argv.to, '--to'), date: argv.date, limit: argv.limit, errors_by: argv.by,
        package: argv.package, user: argv.user, source: argv.source})}`);
      printOutput(output === 'json' ? metrics?.errors : metrics?.errors?.top ?? [], output);
      return true;
    }
  }
  if (command === 'logs')
    return logs(dapi, verb, rest, argv, output);
  if (command === 'timeline') {
    if (argv.session === undefined)
      throw new Error('Usage: grok s observe timeline --session <id> [--from <ISO|-15m>] [--to <ISO>] [--limit 500]');
    const window = argv.from === undefined ? {from: iso('-15m', '--from'), to: new Date().toISOString()}
      : {from: iso(argv.from, '--from'), to: iso(argv.to, '--to')};
    printOutput(await dapi.raw('GET', `/log/timeline${buildQuery({session: argv.session, ...window,
      limit: argv.limit})}`), output);
    return true;
  }
  console.log(OBSERVE_USAGE);
  return true;
}

async function logs(dapi: NodeDapi, verb: string | undefined, rest: string[], argv: any,
                    output: OutputFormat): Promise<boolean> {
  const connection = argv.connection;
  if (verb === 'cloud') {
    if (argv.since !== undefined && argv.from !== undefined)
      throw new Error('--since is not used with --from');
    const group = argv.group ?? (await dapi.raw('GET', `/log/cloud/groups${buildQuery({connection})}`))?.[0];
    if (!group)
      throw new Error('No log group to read: pass --group');
    const window = argv.from !== undefined ? {start: iso(argv.from, '--from'), end: iso(argv.to, '--to')}
      : {start: new Date(Date.now() - parseDuration(argv.since ?? '1h', '--since')).toISOString()};
    printOutput(await dapi.raw('GET', `/log/cloud/events${buildQuery({connection, group, ...window,
      filter: argv.filter, limit: argv.limit, format: 'json'})}`), output);
    return true;
  }
  if (verb === 'archive' && rest[0] === 'list') {
    let objects: any[] = await dapi.raw('GET', `/log/archive/objects${buildQuery({connection, prefix: argv.prefix,
      limit: argv.limit, format: 'json'})}`) ?? [];
    if (argv.since !== undefined) {
      const from = Date.now() - parseDuration(argv.since, '--since');
      objects = objects.filter((o) => Date.parse(o.modified) >= from);
    }
    printOutput(objects, output);
    return true;
  }
  if (verb === 'archive' && rest[0] === 'read' && rest[1] !== undefined) {
    printOutput(await dapi.raw('GET', `/log/archive/events${buildQuery({connection, key: rest[1],
      format: 'json'})}`), output);
    return true;
  }
  throw new Error('Usage: grok s observe logs cloud [--group <g>] ... | logs archive list [--prefix <p>] | ' +
    'logs archive read <key>');
}

async function loggerSet(dapi: NodeDapi, argv: any, current: any, save: (body: any) => Promise<any>,
                         output: OutputFormat): Promise<boolean> {
  const sets: string[] = argv.set === undefined ? [] : (Array.isArray(argv.set) ? argv.set : [argv.set]).map(String);
  const entry = argv.user !== undefined || argv.group !== undefined;
  if (!sets.length && !entry)
    throw new Error('Usage: grok s observe logger set --set <path>=<json> ... [--user <login> | --group <name|id>]');
  if (!entry && (argv.for !== undefined || argv.until !== undefined || argv.reason !== undefined))
    throw new Error('--for, --until and --reason apply to a --user or --group entry');
  if (!entry) {
    const body: any = {};
    for (const s of sets) {
      const key = setPath(current, s);
      body[key] = current[key];
    }
    await save(body);
  }
  else {
    const group = argv.user !== undefined
      ? await dapi.groups.resolve(String(argv.user), {personalOnly: true})
      : await dapi.groups.resolve(String(argv.group));
    const groups = {...current.userGroupSettings};
    const settings = {...groups[group.id]};
    for (const s of sets)
      setPath(settings, s);
    if (argv.for !== undefined)
      settings.expiresAt = new Date(Date.now() + parseDuration(argv.for, '--for')).toISOString();
    else if (argv.until !== undefined)
      settings.expiresAt = iso(argv.until, '--until');
    if (argv.reason !== undefined)
      settings.reason = String(argv.reason);
    groups[group.id] = settings;
    await save({userGroupSettings: groups});
  }
  if (output !== 'quiet')
    console.log(entry ? `Saved the logger settings of ${argv.user ?? argv.group}` : 'Saved the logger settings');
  return true;
}
