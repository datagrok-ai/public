/// `grok s observe` (alias `grok s o`): problems, problem rules, logger settings, errors, the timeline and cloud logs,
/// over `/problems`, `/admin/plugins/<key>/settings`, `/alerts/Test_rules`, `/log`, `/admin/metrics`, `/log/timeline`,
/// `/log/cloud/*`, `/log/archive/*` and `/docker/containers`.
import * as fs from 'fs';
import {NodeDapi, UUID_RE, buildQuery, tryParseJson} from '../utils/node-dapi';
import {printOutput, OutputFormat} from '../utils/server-output';
import {ALL_USERS} from '../utils/migrate/walker';

export const OBSERVE_USAGE = `Usage: grok s observe <command> <verb> [args]      (alias: grok s o)

  problems list [--status <s>] [--state <s>] [--alert-status <s>] [--kind <k>] [--key <k>] [--since 7d] [--limit 100]
                                       Comma lists or all; --alert-status open lists the open alerts
  problems get <id|kind:key>           One field per line, details by path
  problems history <id|kind:key> [--limit 50]
  problems status <id|kind:key> --data '{"status": "muted", "reason": "...", "until": "<ISO>"}' | --json body.json
  problems severity <id|kind:key> <info|warning|critical|none> [--reason <text>]
                                       Pin the severity its raises use; none unpins it
  problems ack <id|kind:key> --reason <text>
  problems resolve <id|kind:key> --reason <text>
  rules get                            The problem rules (Settings > Alerts), as JSON
  rules put --json rules.json          Replace them; the server checks them first
  rules test                           Evaluate the saved rules once, now
  logger get [<path>]                  The logger settings, or one value (exportSettings.0.types)
  logger set --set <path>=<json> ...   Change values; with --user <login> or --group <name|id>, that
             [--for 30m | --until <ISO>] [--reason <text>]   entry of userGroupSettings (a new one starts from
             [--remove]                                       All users), removed when it expires; --remove deletes it
  logger history [--limit 20]          The log-settings-changed audit records, with each change's diff
  errors top [--since 24h | --date "this week" | --from <ISO> [--to <ISO>]] [--by message|signature]
             [--package <name>] [--user <login>] [--source <s>] [--limit 10]
                                       Ranked by distinct users, so errors nobody saw (scheduler, health
                                       probes) come last: narrow those with --source server
  errors get <signature|problem> [--since 7d] [--limit 50]
                                       Each occurrence of an error (a signature from errors top, or the
                                       id or kind:key of its problem): time, user, function, session
  logs container [<name|id> [--since 1h | --from <ISO> [--to <ISO>]] [--grep <regex>] [--limit 20000]]
                                       A docker container's log lines (its latest --limit, then the
                                       window and --grep); no name, the containers and their status
  timeline --session <id> [--from <ISO|-15m>] [--to <ISO>] [--limit 500]
                                       A session's events and requests in a window: no --from, the last
                                       15 min; --to defaults to --from + 10 min, at most 2 h
  logs cloud [--group <g>] [--since 1h | --from <ISO> [--to <ISO>]] [--filter <pattern>] [--limit 1000]
                                       CloudWatch events; no --group, the instance's own group;
                                       --filter is a CloudWatch filter pattern
  logs archive list [--prefix <p>] [--since 7d] [--limit n]
                                       Objects in the log archive (WORM bucket), by key; --since filters
                                       only the page the server returned, so narrow it with --prefix
  logs archive read <key>              One archived object decoded to events
                                       logs cloud and logs archive take [--connection <id>]; without it,
                                       the instance's AWS role
  logs capture [<login|group> [--for 30m] [--reason <text>]]
                                       Save every API request of a user or a group for a while: logger set
                                       of saveHttpRequests on their expiring entry (an entry that already
                                       expires later keeps its time); no argument, the active captures;
                                       read them with timeline --session <id>
  logs capture stop <login|group>      Remove that expiring entry; one that holds other settings is left
                                       to logger set --remove
  logs request <logId>                 What the requests of one command logged (grok s ... --log prints
                                       its logId), gathered from every instance's memory; once only

Durations are 30m (minutes), 24h, 7d, 2w; times are ISO 8601 (UTC without an offset) or
relative, written with = (--from=-2h).

Any grok s command takes --log [<spec>]: its requests share one logId, every instance keeps what they
log at <spec> (verbose, standard, levels, debug flags or all, comma-separated) in memory, and the
command gathers and prints it to stderr at the end, as logs request does; functions run --debug adds
the debug call's own log.`;

const UNIT_MS: Record<string, number> = {m: 60000, h: 3600000, d: 86400000, w: 604800000};

/** `30m`, `2h`, `7d`, `1w` in milliseconds; `m` is minutes, unlike `grok s pull --since`, where it is months. */
export function parseDuration(value: any, flag: string): number {
  const m = /^(\d+)([mhdw])$/.exec(String(value ?? '').trim());
  if (!m || Number(m[1]) <= 0)
    throw new Error(`${flag} expects a duration such as 30m, 2h, 7d or 1w, got '${value}'`);
  return Number(m[1]) * UNIT_MS[m[2]];
}

/** ISO 8601 (UTC without an offset), or relative to [now] (`-2h`). */
export function parseTime(value: any, flag: string, now: Date = new Date()): Date {
  let s = String(value ?? '').trim();
  if (s.startsWith('-'))
    return new Date(now.getTime() - parseDuration(s.substring(1), flag));
  if (/[T ].*[+-]\d{2}$/.test(s))
    s = `${s}:00`;
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

/** A record as `field, value` rows, nested objects by dotted path, so a table shows every value in full. */
export function fieldRows(record: any, prefix: string = ''): {field: string, value: any}[] {
  const rows: {field: string, value: any}[] = [];
  for (const [key, value] of Object.entries(record ?? {}))
    if (value !== null && typeof value === 'object' && !Array.isArray(value))
      rows.push(...fieldRows(value, `${prefix}${key}.`));
    else
      rows.push({field: `${prefix}${key}`, value});
  return rows;
}

export async function handleObserve(dapi: NodeDapi, command: string | undefined, verb: string | undefined,
                                    rest: string[], argv: any, output: OutputFormat): Promise<boolean> {
  // A `:` cannot pass through a route parameter behind the proxy, so kind:key is looked up by the list's filters.
  let id: string | undefined;
  const needId = async () => {
    const given = rest[0] === undefined ? '' : String(rest[0]);
    const colon = given.indexOf(':');
    if (!given || colon === 0 || colon === given.length - 1)
      throw new Error(`Usage: grok s observe ${command} ${verb} <id|kind:key>`);
    if (colon < 0) {
      id = given;
      return;
    }
    const found = await dapi.problems.list({kind: given.substring(0, colon), key: given.substring(colon + 1), limit: 1});
    if (!found?.length)
      throw new Error(`No problem '${given}'`);
    id = found[0].id;
  };
  const problemPath = () => `/problems/${encodeURIComponent(id!)}`;
  // The route answers JSON as text/plain.
  const settings = async (key: string) => {
    const r = await dapi.raw('GET', `/admin/plugins/${key}/settings`);
    return (typeof r === 'string' ? JSON.parse(r) : r)?.settings ?? {};
  };
  // A save within 5 s of the previous one is refused (429) until that one is applied.
  const saveSettings = async (key: string, body: any): Promise<any> => {
    for (let attempt = 0; ; attempt++) {
      try {
        return await dapi.raw('POST', `/admin/plugins/${key}/settings`, body);
      }
      catch (e: any) {
        if (attempt >= 2 || e?.apiError?.errorCode !== 429)
          throw e;
        await new Promise((resolve) => setTimeout(resolve, 5000));
      }
    }
  };
  const saveLogger = (body: any) => saveSettings('logger', body);

  switch (`${command} ${verb}`) {
    case 'problems list':
      printOutput(await dapi.problems.list({status: argv.status, state: argv.state, alertStatus: argv['alert-status'],
        kind: argv.kind, key: argv.key, since: argv.since, limit: argv.limit}), output);
      return true;
    case 'problems get': {
      await needId();
      const problem = await dapi.problems.find(id!);
      if (!problem)
        throw new Error(`No problem '${rest[0]}'`);
      printOutput(output === 'table' ? fieldRows(problem) : problem, output);
      return true;
    }
    case 'problems history':
      await needId();
      printOutput(await dapi.raw('GET', `${problemPath()}/history${buildQuery({limit: argv.limit})}`), output);
      return true;
    case 'problems status':
      await needId();
      printOutput(await dapi.raw('POST', `${problemPath()}/status`, readBody(argv)), output);
      return true;
    case 'problems severity': {
      await needId();
      const severity = rest[1];
      if (!['info', 'warning', 'critical', 'none'].includes(severity))
        throw new Error('Usage: grok s observe problems severity <id|kind:key> <info|warning|critical|none> [--reason <text>]');
      printOutput(await dapi.raw('POST', `${problemPath()}/severity`,
        {severity: severity === 'none' ? null : severity, reason: argv.reason}), output);
      return true;
    }
    case 'problems ack':
    case 'problems resolve':
      await needId();
      printOutput(await dapi.raw('POST', `${problemPath()}/${verb}`, {reason: argv.reason}), output);
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
      return loggerSet(dapi, argv, await settings('logger'), saveLogger, output);
    case 'logger history': {
      const text = encodeURIComponent('eventType.name = "log-settings-changed"');
      const events: any[] = await dapi.raw('GET', `/log?text=${text}&limit=${argv.limit ?? 20}&page=1`) ?? [];
      printOutput(events.map((e) => {
        const params: Record<string, any> = {};
        for (const p of e.parameters ?? [])
          params[p.parameter?.name] = p.value;
        const row = {time: e.eventTime, user: params.user, changed: params.changed};
        return output === 'table' ? row : {...row, diff: params.diff ? tryParseJson(params.diff) ?? params.diff : null};
      }), output);
      return true;
    }
    case 'errors top': {
      if (argv.since !== undefined && (argv.from !== undefined || argv.date !== undefined))
        throw new Error('--since is not used with --from or --date');
      if (argv.date !== undefined && argv.from !== undefined)
        throw new Error('--date is not used with --from');
      const metrics = await dapi.raw('GET', `/admin/metrics${buildQuery({
        date_start: argv.since !== undefined ? new Date(Date.now() - parseDuration(argv.since, '--since')).toISOString()
          : iso(argv.from, '--from'),
        date_end: iso(argv.to, '--to'), date: argv.date, limit: argv.limit, errors_by: argv.by,
        package: argv.package, user: argv.user, source: argv.source})}`);
      // The server answers the top `limit` of every source, in the global order.
      const top = (metrics?.errors?.top ?? []).slice(0, Number(argv.limit ?? 10));
      printOutput(output === 'json' ? {...metrics?.errors, top} : top, output);
      return true;
    }
    case 'errors get': {
      if (rest[0] === undefined)
        throw new Error('Usage: grok s observe errors get <signature|problem id|kind:key> [--since 7d] [--limit 50]');
      const given = String(rest[0]);
      if (given.includes(':'))
        await needId();
      const problem = await dapi.problems.find(id ?? given);
      const signature = problem ? problem.details?.signatureId : given;
      if (!signature)
        throw new Error(`Problem '${given}' is not about an error: it has no signatureId`);
      // An error with a stack is grouped by its stack hash, one without by its event type (ErrorIncidentProblem).
      const text = `(eventType.id = "${signature}" or eventType.errorStackTraceHash = "${signature}") and ` +
        `eventTime > "${new Date(Date.now() - parseDuration(argv.since ?? '7d', '--since')).toISOString()}"`;
      const events: any[] = await dapi.raw('GET', `/log${buildQuery({text, limit: argv.limit ?? 50, page: 1})}`) ?? [];
      const logins = new Map<string, string>();
      for (const e of events)
        for (const p of e.parameters ?? [])
          if (p.parameter?.name === 'user' && p.value && !logins.has(p.value))
            logins.set(p.value, (await dapi.users.find(p.value))?.login ?? p.value);
      printOutput(events.map((e) => {
        const params: Record<string, any> = {};
        for (const p of e.parameters ?? [])
          params[p.parameter?.name] = p.value;
        const row = {time: e.eventTime, source: params.source, user: logins.get(params.user) ?? null,
          function: params.function ?? null, session: e.session?.id ?? null, message: e.description};
        return output === 'table' ? row : {...row, eventId: e.id, params};
      }), output);
      return true;
    }
  }
  if (command === 'logs' && verb === 'capture')
    return capture(dapi, rest, argv, await settings('logger'), saveLogger, output);
  if (command === 'logs')
    return logs(dapi, verb, rest, argv, output);
  if (command === 'timeline') {
    if (argv.session === undefined)
      throw new Error('Usage: grok s observe timeline --session <id> [--from <ISO|-15m>] [--to <ISO>] [--limit 500]');
    needFrom(argv);
    const window = argv.from === undefined ? {from: iso('-15m', '--from'), to: new Date().toISOString()}
      : {from: iso(argv.from, '--from'), to: iso(argv.to, '--to')};
    printOutput(await dapi.raw('GET', `/log/timeline${buildQuery({session: argv.session, ...window,
      limit: argv.limit})}`), output);
    return true;
  }
  throw new Error(`Unknown command '${[command, verb].filter(Boolean).join(' ')}'\n${OBSERVE_USAGE}`);
}

const needFrom = (argv: any) => {
  if (argv.to !== undefined && argv.from === undefined)
    throw new Error('--to needs --from');
};

async function logs(dapi: NodeDapi, verb: string | undefined, rest: string[], argv: any,
                    output: OutputFormat): Promise<boolean> {
  const connection = argv.connection;
  if (verb === 'cloud') {
    if (argv.since !== undefined && argv.from !== undefined)
      throw new Error('--since is not used with --from');
    needFrom(argv);
    const group = argv.group ?? (await dapi.raw('GET', `/log/cloud/groups${buildQuery({connection})}`))?.[0];
    if (!group)
      throw new Error('No log group to read: pass --group');
    const window = argv.from !== undefined ? {start: iso(argv.from, '--from'), end: iso(argv.to, '--to')}
      : {start: new Date(Date.now() - parseDuration(argv.since ?? '1h', '--since')).toISOString()};
    printOutput(await dapi.raw('GET', `/log/cloud/events${buildQuery({connection, group, ...window,
      filter: argv.filter, limit: argv.limit, format: 'json'})}`), output);
    return true;
  }
  if (verb === 'container')
    return containerLogs(dapi, rest[0], argv, output);
  if (verb === 'request' && rest[0] !== undefined) {
    const rows = await readRequestLog(dapi, String(rest[0]));
    printOutput(output === 'table' ? rows.map((r) => ({time: r.time, level: r.level, flag: r.flag ?? '',
      source: r.source ?? '', message: r.message})) : rows, output);
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
  throw new Error('Usage: grok s observe logs cloud [--group <g>] ... | logs container [<name|id>] | ' +
    'logs archive list [--prefix <p>] | logs archive read <key> | logs request <logId> | logs capture [<login|group>] | ' +
    'logs capture stop <login|group>');
}

/// `GET /log/requests/<logId>`: what every isolate of every instance kept in memory for the id, as JSON lines,
/// oldest first. Gathering hands the records over, so a second read returns none.
async function readRequestLog(dapi: NodeDapi, logId: string): Promise<any[]> {
  return await dapi.raw('GET', `/log/requests/${encodeURIComponent(logId)}`) ?? [];
}

/** After a command run with `--log`: its log on stderr, so stdout keeps the result; JSON lines under `--output json`. */
export async function printRequestLog(dapi: NodeDapi, logId: string, output: OutputFormat): Promise<void> {
  let rows: any[];
  try {
    rows = await readRequestLog(dapi, logId);
  }
  catch (e: any) {
    process.stderr.write(`--log: could not read the log of ${logId}: ${e?.message ?? e}\n`);
    return;
  }
  if (output === 'json') {
    for (const r of rows)
      process.stderr.write(JSON.stringify({logId, ...r}) + '\n');
    return;
  }
  process.stderr.write(`── log ${logId}: ${rows.length} record${rows.length === 1 ? '' : 's'} (grok s o logs request ${logId})\n`);
  for (const r of rows) {
    const tags = [r.flag, r.source].filter(Boolean).join(' ');
    process.stderr.write(`${String(r.time ?? '').substring(11, 23)} ${String(r.level ?? '').toUpperCase().padEnd(7)} ` +
      `${tags ? `[${tags}] ` : ''}${r.message ?? ''}\n`);
    if (r.stackTrace)
      process.stderr.write(String(r.stackTrace).replace(/^/gm, '    ') + '\n');
  }
}

/// `logs container`: the server only cuts a container's log by line count, so the window and --grep apply here.
async function containerLogs(dapi: NodeDapi, name: string | undefined, argv: any,
                             output: OutputFormat): Promise<boolean> {
  if (name === undefined) {
    const all: any[] = await dapi.raw('GET', `/docker/containers${buildQuery({limit: 1000})}`) ?? [];
    printOutput(all.map((c) => ({id: c.id, name: c.name, status: c.status})), output);
    return true;
  }
  if (argv.since !== undefined && argv.from !== undefined)
    throw new Error('--since is not used with --from');
  needFrom(argv);
  let id = String(name);
  if (!UUID_RE.test(id)) {
    const found: any[] = await dapi.raw('GET', `/docker/containers${buildQuery({text: id, limit: 100})}`) ?? [];
    const container = found.find((c) => c.name === id);
    if (!container)
      throw new Error(`No container '${id}'${found.length ? `; similar: ${found.map((c) => c.name).join(', ')}` : ''}`);
    id = container.id;
  }
  const from = argv.from !== undefined ? parseTime(argv.from, '--from').getTime()
    : argv.since !== undefined ? Date.now() - parseDuration(argv.since, '--since') : -Infinity;
  const to = argv.to !== undefined ? parseTime(argv.to, '--to').getTime() : Infinity;
  const grep = argv.grep === undefined ? null : new RegExp(String(argv.grep), 'i');
  const text = await dapi.raw('GET', `/docker/containers/${id}/logs${buildQuery({limit: argv.limit ?? 20000})}`);
  // A line without its own timestamp (a stack frame) belongs to the window of the line before it.
  let inWindow = from === -Infinity && to === Infinity;
  const lines: string[] = [];
  for (const line of String(text ?? '').split('\n')) {
    const stamp = /^(\d{4}-\d{2}-\d{2}T[\d:]{8})(?:\.\d+)?(Z|[+-]\d{2}:\d{2})/.exec(line);
    if (stamp) {
      const t = Date.parse(`${stamp[1]}${stamp[2]}`);
      inWindow = t >= from && t <= to;
    }
    if (inWindow && line && (!grep || grep.test(line)))
      lines.push(line);
  }
  if (output === 'json')
    printOutput(lines, output);
  else
    for (const line of lines)
      console.log(line);
  return true;
}

async function loggerSet(dapi: NodeDapi, argv: any, current: any, save: (body: any) => Promise<any>,
                         output: OutputFormat): Promise<boolean> {
  const sets: string[] = argv.set === undefined ? [] : (Array.isArray(argv.set) ? argv.set : [argv.set]).map(String);
  const entry = argv.user !== undefined || argv.group !== undefined || argv.target !== undefined;
  if (!sets.length && !entry)
    throw new Error('Usage: grok s observe logger set --set <path>=<json> ... [--user <login> | --group <name|id>]');
  if (!entry && (argv.for !== undefined || argv.until !== undefined || argv.reason !== undefined || argv.remove))
    throw new Error('--for, --until, --reason and --remove apply to a --user or --group entry');
  const name = argv.user ?? argv.group ?? argv.name;
  if (!entry) {
    const body: any = {};
    for (const s of sets) {
      const key = setPath(current, s);
      body[key] = current[key];
    }
    await save(body);
  }
  else {
    const group = argv.target ?? (argv.user !== undefined
      ? await dapi.groups.resolve(String(argv.user), {personalOnly: true})
      : await dapi.groups.resolve(String(argv.group)));
    const groups = {...current.userGroupSettings};
    const existing = groups[group.id];
    if (argv.remove) {
      if (!existing)
        throw new Error(`${name} has no logger entry`);
      delete groups[group.id];
    }
    else {
      if (existing && !existing.expiresAt && (argv.for !== undefined || argv.until !== undefined))
        throw new Error(`${name} has a logger entry that does not expire; change it without --for or --until`);
      // A user's settings come from one entry, so a new one starts from All users instead of the defaults.
      const settings = existing ? {...existing} : {...groups[ALL_USERS]};
      if (!existing)
        delete settings.reason;
      for (const s of sets)
        setPath(settings, s);
      if (argv.for !== undefined)
        settings.expiresAt = new Date(Date.now() + parseDuration(argv.for, '--for')).toISOString();
      else if (argv.until !== undefined) {
        settings.expiresAt = iso(argv.until, '--until');
        if (Date.parse(settings.expiresAt) <= Date.now())
          throw new Error(`--until ${argv.until} is in the past`);
      }
      if (argv.reason !== undefined)
        settings.reason = String(argv.reason);
      groups[group.id] = settings;
    }
    await save({userGroupSettings: groups});
  }
  if (output !== 'quiet')
    console.log(!entry ? 'Saved the logger settings' : argv.remove ? `Removed the logger entry of ${name}`
      : `Saved the logger settings of ${name}`);
  return true;
}

/// `logs capture`: logger set of saveHttpRequests, and the usage level it is saved at, on an expiring entry.
async function capture(dapi: NodeDapi, rest: string[], argv: any, current: any,
                       save: (body: any) => Promise<any>, output: OutputFormat): Promise<boolean> {
  const groups = current.userGroupSettings ?? {};
  const stop = rest[0] === 'stop';
  const name = stop ? rest[1] : rest[0];
  if (name === undefined) {
    if (stop)
      throw new Error('Usage: grok s observe logs capture stop <login|group>');
    const rows = [];
    for (const [id, s] of Object.entries<any>(groups))
      if (s?.saveHttpRequests === true && s.expiresAt) {
        const group = await dapi.groups.find(id);
        rows.push({id, name: group?.friendlyName ?? group?.name ?? id, expiresAt: s.expiresAt, reason: s.reason ?? ''});
      }
    printOutput(rows, output);
    return true;
  }
  let group: any;
  try {
    group = await dapi.groups.resolve(String(name), {personalOnly: true});
  }
  catch {
    group = await dapi.groups.resolve(String(name));
  }
  const existing = groups[group.id];
  if (stop) {
    if (!existing?.saveHttpRequests)
      throw new Error(`No capture for ${name}`);
    if (!existing.expiresAt)
      throw new Error(`The logger entry of ${name} does not expire; change it with logger set`);
    const own = ['saveHttpRequests', 'postLevels', 'saveLevels', 'expiresAt', 'reason'];
    const others = Object.keys(existing).filter((k) => !own.includes(k) &&
      JSON.stringify(existing[k]) !== JSON.stringify(groups[ALL_USERS]?.[k]));
    if (others.length)
      throw new Error(`The logger entry of ${name} also holds ${others.join(', ')}; remove it with ` +
        `grok s o logger set --${group.personal ? 'user' : 'group'} ${name} --remove`);
    return loggerSet(dapi, {target: group, name, remove: true}, current, save, output);
  }
  const base = existing ?? groups[ALL_USERS] ?? {};
  const sets = ['saveHttpRequests=true', ...['postLevels', 'saveLevels']
    .filter((l) => Array.isArray(base[l]) && !base[l].includes('usage'))
    .map((l) => `${l}=${JSON.stringify([...base[l], 'usage'])}`)];
  const end = Math.max(Date.now() + parseDuration(argv.for ?? '30m', '--for'), Date.parse(existing?.expiresAt) || 0);
  return loggerSet(dapi, {target: group, name, set: sets, reason: argv.reason, until: new Date(end).toISOString()},
    current, save, output);
}
