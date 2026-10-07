/// `grok s observe alerts ...` — the alerts the deployment's problems raised (AlertsRouter, `/alerts`), and
/// `grok s observe problems ...` — what the deployment detects and what people decided about it (`/problems`).
import {Query} from '../utils/node-observability';
import {Connect, eachHost, forEachHost, hostList} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {fmtTime, fmtDateTime, hasValue, optString, parseDuration, parseTime, printBlock, rows, sinceArg,
  truncate} from '../utils/obs-format';

export const ALERTS_USAGE = `Usage: grok s observe alerts <verb> [args]
  list [--status open,acknowledged|resolved|all] [--kind <k>] [--since 24h] [--limit n] [--host a --host b ...]
  get <id|kind:key>
  ack <id|kind:key> [--reason <text>]
  resolve <id|kind:key> [--reason <text>]   a condition that still holds raises a new alert on its next check
  mute <id|kind:key> --reason <text> [--for 2h | --until <iso|HH:MM> | --until-version <v>]
                                           stops the alerts of the alert's problem and resolves this one;
                                           no end: until unmuted
  unmute <id|kind:key> [--reason <text>]   the alert's problem alerts again
  detection [--all] [--host a --host b ...]      live servers and those stopped or last seen within 1 h;
                                                --all lists every server; hosts on one database print once
An alert is one message that a problem went wrong. It stays open until a person resolves it; CLEARED says its
condition ended. Resolving leaves the problem active, so it alerts again: to stop its alerts, mute it.
Ids: a UUID, a unique UUID prefix (6+ characters) or kind:key (connection:ELN:Prod), which names an open or
acknowledged alert (its key exactly, else a unique key prefix); unmute also finds a resolved one.`;

export const PROBLEMS_USAGE = `Usage: grok s observe problems <verb> [args]
  list [--status active,muted,not-a-problem,fixed|all] [--state ongoing|cleared] [--kind <k>] [--since 7d]
       [--limit n] [--host a --host b ...]
  get <id|kind:key>
  history <id|kind:key> [--limit n]   the problem's records: alerts opened, escalated, acknowledged, cleared,
                                      resolved, and status changes
  mute <id|kind:key> --reason <text> [--for 2h | --until <iso|HH:MM> | --until-version <v>]
                                        stops its alerts and resolves the open one; no end: until activated
  dismiss <id|kind:key> --reason <text>   not a problem: never alerts again
  fix <id|kind:key> [--reason <text>]     fixed: alerts again, as a regression, if it comes back
  activate <id|kind:key> [--reason <text>]   alerts again (unmutes)
A problem is what is wrong, kept for good; its status decides whether it alerts. Muting a problem is how
you stop its alerts; activating it lets it alert again. Active problems raise an alert when they start;
muted ones and those that are not a problem never do. Every status but active resolves the open alert. Ids: a UUID, a unique UUID prefix (6+ characters) or kind:key (connection:ELN:Prod).`;

const MUTE_USAGE = 'mute <id|kind:key> --reason <text> [--for 2h | --until <iso|HH:MM> | --until-version <v>]';
const ALERT_PAST: Record<string, string> = {ack: 'acknowledged', unmute: 'unmuted', resolve: 'resolved'};
const STATUS: Record<string, string> = {dismiss: 'not-a-problem', fix: 'fixed', activate: 'active'};
const PROBLEM_PAST: Record<string, string> = {dismiss: 'dismissed as not a problem', fix: 'marked fixed', activate: 'made active'};

export async function handleAlerts(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                   output: OutputFormat): Promise<boolean> {
  const id = optString(rest[0]);
  const alert = async (status: string = 'open,acknowledged') => {
    const alerts = (await connect()).alerts;
    const all = status === 'all';
    const list = (q: Query) => alerts.list({...q, status, limit: all ? 1 : undefined});
    return {alerts, id: await target(id!, list, all ? 'alert' : 'open alert')};
  };
  switch (verb) {
    case 'list': {
      const since = argv.since === undefined ? undefined : sinceArg(argv.since);
      const q = {status: optString(argv.status), kind: optString(argv.kind), since, limit: argv.limit};
      printOutput(await forEachHost(argv, connect, async (dapi) =>
        rows(await dapi.alerts.list(q) ?? [], output, alertRow), output), output);
      return true;
    }
    case 'detection': {
      const deployments: {hosts: string[]; ids: string[]; lease: any}[] = [];
      await eachHost(argv, connect, async (dapi, host) => {
        const lease = await dapi.alerts.detection();
        const servers: any[] = lease?.servers ?? [];
        const ids = servers.map((s) => String(s?.id));
        const same = deployments.find((d) => d.ids.some((x) => ids.includes(x)));
        if (same)
          same.hosts.push(host);
        else
          deployments.push({hosts: [host], ids, lease: {...lease, servers: argv.all === true ? servers
            : servers.filter((s) => s?.live || Date.now() - Date.parse(s?.stoppedAt ?? s?.lastSeen) <= 3600000)}});
      });
      if (hostList(argv.host).length < 2)
        printOutput(output === 'json' ? deployments[0].lease : detectionRows(deployments[0].lease), output);
      else if (output === 'json')
        printOutput(deployments.map((d) => ({hosts: d.hosts, ...d.lease})), output);
      else
        printOutput(deployments.flatMap((d) => detectionRows(d.lease).map((r) => ({HOST: d.hosts.join(', '), ...r}))), output);
      return true;
    }
    case 'get': {
      if (!id) return usage('alerts', 'get <id|kind:key>');
      const t = await alert();
      const a = await t.alerts.get(t.id);
      if (output === 'table') printAlert(a);
      else printOutput(a, output);
      return true;
    }
    case 'ack':
    case 'unmute':
    case 'resolve': {
      if (!id) return usage('alerts', `${verb} <id|kind:key> [--reason <text>]`);
      const reason = optString(argv.reason);
      const t = await alert(verb === 'unmute' ? 'all' : undefined);
      const a = await t.alerts.transition(t.id, verb, {reason});
      report(a, `${ALERT_PAST[verb]} ${identity(a, id)}${reason ? ` — ${reason}` : ''}`, output);
      return true;
    }
    case 'mute': {
      if (!id) return usage('alerts', MUTE_USAGE);
      const {body, until} = muteBody(argv);
      const t = await alert();
      const a = await t.alerts.transition(t.id, 'mute', body);
      report(a, `muted ${identity(a, id)} ${until} — ${body.reason}`, output);
      return true;
    }
  }
  printError(new Error(ALERTS_USAGE));
  return false;
}

export async function handleProblems(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                     output: OutputFormat): Promise<boolean> {
  const id = optString(rest[0]);
  const problem = async () => {
    const alerts = (await connect()).alerts;
    return {alerts, id: await target(id!, (q) => alerts.problems({...q, status: 'all'}), 'problem')};
  };
  switch (verb) {
    case 'list': {
      const since = argv.since === undefined ? undefined : sinceArg(argv.since);
      const q = {status: optString(argv.status), state: optString(argv.state), kind: optString(argv.kind), since, limit: argv.limit};
      printOutput(await forEachHost(argv, connect, async (dapi) =>
        rows(await dapi.alerts.problems(q) ?? [], output, problemRow), output), output);
      return true;
    }
    case 'get': {
      if (!id) return usage('problems', 'get <id|kind:key>');
      const t = await problem();
      const p = await t.alerts.problem(t.id);
      if (output === 'table') printProblem(p);
      else printOutput(p, output);
      return true;
    }
    case 'history': {
      if (!id) return usage('problems', 'history <id|kind:key> [--limit n]');
      const t = await problem();
      printOutput(rows(await t.alerts.history(t.id, {limit: argv.limit}) ?? [], output, historyRow), output);
      return true;
    }
    case 'mute': {
      if (!id) return usage('problems', MUTE_USAGE);
      const {body, until} = muteBody(argv);
      const t = await problem();
      const p = await t.alerts.setStatus(t.id, {status: 'muted', ...body});
      report(p, `muted ${identity(p, id)} ${until} — ${body.reason}`, output);
      return true;
    }
    case 'dismiss':
    case 'fix':
    case 'activate': {
      const reason = optString(argv.reason);
      if (!id || (verb === 'dismiss' && !reason))
        return usage('problems', `${verb} <id|kind:key> ${verb === 'dismiss' ? '--reason <text>' : '[--reason <text>]'}`);
      const t = await problem();
      const p = await t.alerts.setStatus(t.id, {status: STATUS[verb], reason});
      report(p, `${PROBLEM_PAST[verb]}: ${identity(p, id)}${reason ? ` — ${reason}` : ''}`, output);
      return true;
    }
  }
  printError(new Error(PROBLEMS_USAGE));
  return false;
}

async function target(id: string, list: (q: Query) => Promise<any[]>, noun: string): Promise<string> {
  const colon = id.indexOf(':');
  if (colon < 0) return id;
  const matches: any[] = await list({kind: id.slice(0, colon), key: id.slice(colon + 1)}) ?? [];
  if (!matches.length)
    throw new Error(`No ${noun} ${id}`);
  if (matches.length > 1)
    throw new Error(`Several ${noun}s match ${id}; pass an id:\n` +
      matches.map((m) => `  ${m?.id}  ${m?.kind}:${m?.key}  ${m?.status}  ${truncate(m?.summary, 60)}`).join('\n'));
  return String(matches[0].id);
}

function usage(command: string, line: string): boolean {
  printError(new Error(`Usage: grok s observe ${command} ${line}`));
  return false;
}

function identity(x: any, fallback: string): string {
  return x?.kind ? `${x.kind}:${x.key}` : fallback;
}

function report(x: any, line: string, output: OutputFormat): void {
  output === 'table' ? console.log(line) : printOutput(x, output);
}

/** A reason and at most one of `--for`, `--until`, `--until-version`; none (or `--forever`) mutes until a person lifts it. */
export function muteBody(argv: any, now: Date = new Date()): {body: Record<string, any>; until: string} {
  const given = ['for', 'until', 'until-version'].filter((k) => k in argv);
  if (given.length > 1 || given.some((k) => !hasValue(argv[k])) || (given.length && argv.forever === true))
    throw new Error('mute takes at most one of --for <duration>, --until <iso|HH:MM>, --until-version <v>');
  const reason = optString(argv.reason);
  if (!reason)
    throw new Error('mute needs --reason <text>');
  if (!given.length)
    return {body: {reason}, until: 'until lifted'};
  if (hasValue(argv['until-version']))
    return {body: {reason, untilVersion: String(argv['until-version'])}, until: `until ${argv['until-version']}`};
  const until = hasValue(argv.for)
    ? new Date(now.getTime() + parseDuration(argv.for, '--for'))
    : parseTime(argv.until, '--until', now);
  if (until.getTime() <= now.getTime())
    throw new Error(`--until ${argv.until} is in the past`);
  return {body: {reason, until: until.toISOString()}, until: `until ${fmtTime(until, now)}`};
}

export function alertRow(a: any): Record<string, any> {
  return {
    KIND: a?.kind ?? '',
    KEY: a?.key ?? '',
    SEV: a?.severity ?? '',
    AUDIENCE: a?.audience ?? '',
    STATUS: a?.status ?? '',
    OPENED: fmtTime(a?.openedAt),
    CLEARED: fmtTime(a?.clearedAt),
    SUMMARY: truncate(a?.summary, 60),
  };
}

export function problemRow(p: any): Record<string, any> {
  return {
    KIND: p?.kind ?? '',
    KEY: p?.key ?? '',
    SEV: p?.severity ?? '',
    STATUS: p?.status ?? '',
    STATE: p?.state ?? '',
    EPISODES: p?.episodes ?? 0,
    'LAST SEEN': fmtTime(p?.lastSeen),
    SUMMARY: truncate(p?.summary, 60),
  };
}

export function historyRow(r: any): Record<string, any> {
  return {
    TIME: fmtDateTime(r?.time),
    RECORD: r?.type ?? '',
    ALERT: truncate(r?.alertId, 8),
    SUMMARY: truncate(r?.summary, 80),
  };
}

export function detectionRows(detection: any): Record<string, any>[] {
  return (detection?.servers ?? []).map((s: any) => ({
    SERVER: s?.name ?? '',
    'HOST NAME': s?.host ?? '',
    VERSION: s?.version ?? '',
    'LAST SEEN': fmtTime(s?.lastSeen),
    LIVE: s?.live ? 'yes' : 'no',
  }));
}

function printAlert(a: any): void {
  const lines: [string, string][] = [
    ['alert', `${a?.kind}:${a?.key}  ${a?.alertname ?? ''}`],
    ['status', `${a?.status ?? ''}  ${a?.severity ?? ''}  audience ${a?.audience ?? ''}`],
    ['summary', a?.summary ?? ''],
    ['opened', `${fmtDateTime(a?.openedAt)}` +
      `  · last seen ${fmtDateTime(a?.lastSeen)}  · ${a?.occurrences ?? 0} occurrences`],
  ];
  if (a?.clearedAt) lines.push(['cleared', `${fmtDateTime(a.clearedAt)} — the condition ended; open until resolved`]);
  if (a?.ackedAt) lines.push(['acknowledged', fmtDateTime(a.ackedAt)]);
  if (a?.resolvedAt) lines.push(['resolved', `${fmtDateTime(a.resolvedAt)}${a?.resolveReason ? ` — ${a.resolveReason}` : ''}`]);
  if (a?.url) lines.push(['url', a.url]);
  if (a?.details) lines.push(['details', JSON.stringify(a.details)]);
  lines.push(['id', a?.id ?? '']);
  if (a?.problemId) lines.push(['problem', a.problemId]);
  printBlock(lines);
}

function printProblem(p: any): void {
  const lines: [string, string][] = [
    ['problem', `${p?.kind}:${p?.key}  ${p?.name ?? ''}`],
    ['status', `${p?.status ?? ''}${p?.statusReason ? ` — ${p.statusReason}` : ''}`],
    ['state', `${p?.state ?? ''}  ${p?.severity ?? ''}  audience ${p?.audience ?? ''}`],
    ['summary', p?.summary ?? ''],
    ['seen', `first ${fmtDateTime(p?.firstSeen)}  · last ${fmtDateTime(p?.lastSeen)}  · ${p?.episodes ?? 0} episodes, ` +
      `${p?.occurrences ?? 0} occurrences`],
  ];
  if (p?.mutedUntil) lines.push(['muted until', fmtDateTime(p.mutedUntil)]);
  if (p?.mutedUntilVersion) lines.push(['muted until', `version ${p.mutedUntilVersion}`]);
  if (p?.url) lines.push(['url', p.url]);
  if (p?.details) lines.push(['details', JSON.stringify(p.details)]);
  lines.push(['id', p?.id ?? '']);
  printBlock(lines);
}
