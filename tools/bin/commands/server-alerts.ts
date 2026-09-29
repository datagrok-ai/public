/// `grok s alerts ...` — the deployment's alert state (AlertsRouter, `/alerts`).
import {NodeAlertsClient} from '../utils/node-observability';
import {Connect, eachHost, forEachHost, hostList, singleHost} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {fmtTime, fmtDateTime, hasValue, parseDuration, parseTime, printBlock, sinceArg, truncate} from '../utils/obs-format';

export const ALERTS_USAGE = `Usage: grok s alerts <verb> [args]
  list [--status open,acknowledged|muted|resolved|all] [--kind <k>] [--since 24h] [--limit n] [--host a --host b ...]
  get <id|kind:key>
  ack <id|kind:key> [--reason <text>]
  mute <id|kind:key> --reason <text> (--for 2h | --until <iso|HH:MM> | --until-version <v> | --forever)
  unmute <id|kind:key> [--reason <text>]
  resolve <id|kind:key> [--reason <text>]
  detection [--all] [--host a --host b ...]      live servers and those stopped or last seen within 1 h;
                                                --all lists every server; hosts on one database print once
Ids: a UUID, a unique UUID prefix (6+ characters) or kind:key (connection:ELN:Prod), which names an open,
acknowledged or muted alert (its key exactly, else a unique key prefix).`;

const PAST: Record<string, string> = {ack: 'acknowledged', unmute: 'unmuted', resolve: 'resolved'};

export async function handleAlerts(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                   output: OutputFormat): Promise<boolean> {
  const id = rest[0] === undefined ? undefined : String(rest[0]);
  switch (verb) {
    case 'list': {
      const since = argv.since === undefined ? undefined : sinceArg(argv.since);
      const q = {status: optString(argv.status), kind: optString(argv.kind), since, limit: argv.limit};
      const rows = await forEachHost(argv, connect, async (dapi) => {
        const alerts: any[] = await dapi.alerts.list(q) ?? [];
        return output === 'json' || output === 'quiet' ? alerts : alerts.map(alertRow);
      }, output);
      printOutput(rows, output);
      return true;
    }
    case 'detection': {
      const deployments: {hosts: string[]; ids: string[]; lease: any}[] = [];
      await eachHost(argv, connect, async (dapi, host) => {
        const lease = await dapi.alerts.detection();
        const ids = serverIds(lease);
        const same = deployments.find((d) => d.ids.some((id) => ids.includes(id)));
        if (same) same.hosts.push(host);
        else deployments.push({hosts: [host], ids, lease: {...lease, servers: recentServers(lease, argv.all === true)}});
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
      if (!id) return usage('get <id|kind:key>');
      const t = await target(connect, argv, verb, id);
      const alert = await t.alerts.get(t.id);
      if (output === 'table') printAlert(alert);
      else printOutput(alert, output);
      return true;
    }
    case 'ack':
    case 'unmute':
    case 'resolve': {
      if (!id) return usage(`${verb} <id|kind:key> [--reason <text>]`);
      const reason = optString(argv.reason);
      const t = await target(connect, argv, verb, id);
      const alert = await t.alerts.transition(t.id, verb, {reason});
      report(alert, `${PAST[verb]} ${identity(alert, id)}${reason ? ` — ${reason}` : ''}`, output);
      return true;
    }
    case 'mute': {
      if (!id) return usage('mute <id|kind:key> --reason <text> (--for 2h | --until <iso|HH:MM> | --until-version <v> | --forever)');
      const {body, until} = muteBody(argv);
      const t = await target(connect, argv, verb, id);
      const alert = await t.alerts.transition(t.id, 'mute', body);
      report(alert, `muted ${identity(alert, id)} ${until} — ${body.reason}`, output);
      return true;
    }
  }
  printError(new Error(ALERTS_USAGE));
  return false;
}

/** The alert an id names; `kind:key` cannot travel in the path, so it is looked up among the open alerts first. */
async function target(connect: Connect, argv: any, verb: string, id: string): Promise<{alerts: NodeAlertsClient; id: string}> {
  const alerts = (await connect(singleHost(argv, `alerts ${verb}`))).alerts;
  return {alerts, id: await alertId(alerts, id)};
}

export async function alertId(alerts: NodeAlertsClient, id: string): Promise<string> {
  const colon = id.indexOf(':');
  if (colon < 0) return id;
  const matches: any[] = await alerts.list({kind: id.slice(0, colon), key: id.slice(colon + 1), status: 'open,acknowledged,muted'}) ?? [];
  if (!matches.length)
    throw new Error(`No open alert ${id}`);
  if (matches.length > 1)
    throw new Error(`Several open alerts match ${id}; pass an id:\n` +
      matches.map((a) => `  ${a?.id}  ${a?.kind}:${a?.key}  ${a?.status}  ${truncate(a?.summary, 60)}`).join('\n'));
  return String(matches[0].id);
}

function usage(line: string): boolean {
  printError(new Error(`Usage: grok s alerts ${line}`));
  return false;
}

function optString(v: any): string | undefined {
  return v === undefined || v === null || v === true ? undefined : String(v);
}

function identity(alert: any, fallback: string): string {
  return alert?.kind ? `${alert.kind}:${alert.key}` : fallback;
}

function report(alert: any, line: string, output: OutputFormat): void {
  if (output === 'json' || output === 'csv') printOutput(alert, output);
  else if (output === 'quiet') console.log(alert?.id ?? '');
  else console.log(line);
}

/** Exactly one of `--for`, `--until`, `--until-version`, `--forever`, and a reason. */
export function muteBody(argv: any, now: Date = new Date()): {body: Record<string, any>; until: string} {
  const given = ['for', 'until', 'until-version'].filter((k) => hasValue(argv[k])).concat(argv.forever === true ? ['forever'] : []);
  if (given.length !== 1)
    throw new Error('alerts mute takes exactly one of --for <duration>, --until <iso|HH:MM>, --until-version <v>, --forever');
  const reason = optString(argv.reason);
  if (!reason)
    throw new Error('alerts mute needs --reason <text>');
  if (argv.forever === true)
    return {body: {reason, forever: true}, until: 'forever'};
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
    BY: a?.openedOnServerName ?? '',
    SUMMARY: truncate(a?.summary, 60),
  };
}

/**
 * Two aliases on one database answer from the same `servers` table, so leases that share a server id
 * are one deployment (the lease holder alone is empty while nothing holds it). The ids are taken before
 * the recent-server filter, which hides different rows as time passes between the calls.
 */
export function serverIds(lease: any): string[] {
  const servers: any[] = Array.isArray(lease?.servers) ? lease.servers : [];
  return servers.map((s) => String(s?.id));
}

/** Live servers, and those that stopped or were last seen within the hour; every row with [all]. */
export function recentServers(lease: any, all: boolean, now: number = Date.now()): any[] {
  const servers: any[] = Array.isArray(lease?.servers) ? lease.servers : [];
  return all ? servers : servers.filter((s) => s?.live || now - Date.parse(s?.stoppedAt ?? s?.lastSeen) <= 3600000);
}

export function detectionRows(lease: any): Record<string, any>[] {
  const servers: any[] = Array.isArray(lease?.servers) ? lease.servers : [];
  return servers.map((s) => ({
    SERVER: s?.name ?? '',
    'HOST NAME': s?.host ?? '',
    VERSION: s?.version ?? '',
    'LAST SEEN': fmtTime(s?.lastSeen),
    LIVE: s?.live ? 'yes' : 'no',
    ELIGIBLE: s?.eligible ? 'yes' : 'no',
    OWNER: s?.id && s.id === lease?.holder ? '*' : '',
  }));
}

function printAlert(a: any): void {
  const lines: [string, string][] = [
    ['alert', `${a?.kind}:${a?.key}  ${a?.alertname ?? ''}`],
    ['status', `${a?.status ?? ''}  ${a?.severity ?? ''}  audience ${a?.audience ?? ''}`],
    ['summary', a?.summary ?? ''],
    ['opened', `${fmtDateTime(a?.openedAt)}${a?.openedOnServerName ? ` on ${a.openedOnServerName}` : ''}` +
      `  · last seen ${fmtDateTime(a?.lastSeen)}  · ${a?.occurrences ?? 0} occurrences`],
  ];
  if (a?.ackedAt) lines.push(['acknowledged', fmtDateTime(a.ackedAt)]);
  if (a?.resolvedAt) lines.push(['resolved', `${fmtDateTime(a.resolvedAt)}${a?.resolveReason ? ` — ${a.resolveReason}` : ''}`]);
  if (a?.url) lines.push(['url', a.url]);
  if (a?.details) lines.push(['details', JSON.stringify(a.details)]);
  lines.push(['id', a?.id ?? '']);
  printBlock(lines);
}
