/// `grok s o problems ...` — what the deployment detects and what people decided about it (`/problems`).
import {NodeAlertsClient} from '../utils/node-observability';
import {Connect, forEachHost, singleHost} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {fmtTime, fmtDateTime, printBlock, sinceArg, truncate} from '../utils/obs-format';
import {alertRow, muteBody} from './server-alerts';

export const PROBLEMS_USAGE = `Usage: grok s o problems <verb> [args]
  list [--status active,muted,not-a-problem,fixed|all] [--state ongoing|cleared] [--kind <k>] [--since 7d]
       [--limit n] [--host a --host b ...]
  get <id|kind:key>
  alerts <id|kind:key> [--status open,acknowledged,resolved|all]   the alerts the problem raised
  mute <id|kind:key> --reason <text> [--for 2h | --until <iso|HH:MM> | --until-version <v>]
                                        no end: muted until a person makes it active
  dismiss <id|kind:key> --reason <text>   not a problem: never alerts again
  fix <id|kind:key> [--reason <text>]     fixed: alerts again, as a regression, if it comes back
  activate <id|kind:key> [--reason <text>]
Statuses: active problems raise an alert when they start; muted ones and those that are not a problem never
do. Every status but active resolves the open alert. Ids: a UUID, a unique UUID prefix (6+ characters) or
kind:key (connection:ELN:Prod).`;

const STATUS: Record<string, string> = {dismiss: 'not-a-problem', fix: 'fixed', activate: 'active'};
const PAST: Record<string, string> = {dismiss: 'dismissed as not a problem', fix: 'marked fixed', activate: 'made active'};

export async function handleProblems(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                     output: OutputFormat): Promise<boolean> {
  const id = rest[0] === undefined ? undefined : String(rest[0]);
  switch (verb) {
    case 'list': {
      const since = argv.since === undefined ? undefined : sinceArg(argv.since);
      const q = {status: optString(argv.status), state: optString(argv.state), kind: optString(argv.kind), since, limit: argv.limit};
      const rows = await forEachHost(argv, connect, async (dapi) => {
        const problems: any[] = await dapi.alerts.problems(q) ?? [];
        return output === 'json' || output === 'quiet' ? problems : problems.map(problemRow);
      }, output);
      printOutput(rows, output);
      return true;
    }
    case 'get': {
      if (!id) return usage('get <id|kind:key>');
      const t = await target(connect, argv, verb, id);
      const problem = await t.alerts.problem(t.id);
      if (output === 'table') printProblem(problem);
      else printOutput(problem, output);
      return true;
    }
    case 'alerts': {
      if (!id) return usage('alerts <id|kind:key>');
      const t = await target(connect, argv, verb, id);
      const alerts: any[] = await t.alerts.list({problem: t.id, status: optString(argv.status) ?? 'all'}) ?? [];
      printOutput(output === 'json' || output === 'quiet' ? alerts : alerts.map(alertRow), output);
      return true;
    }
    case 'mute': {
      if (!id) return usage('mute <id|kind:key> --reason <text> [--for 2h | --until <iso|HH:MM> | --until-version <v>]');
      const {body, until} = muteBody(argv);
      const t = await target(connect, argv, verb, id);
      const problem = await t.alerts.setStatus(t.id, {status: 'muted', ...body});
      report(problem, `muted ${identity(problem, id)} ${until} — ${body.reason}`, output);
      return true;
    }
    case 'dismiss':
    case 'fix':
    case 'activate': {
      const reason = optString(argv.reason);
      if (!id || (verb === 'dismiss' && !reason))
        return usage(`${verb} <id|kind:key> ${verb === 'dismiss' ? '--reason <text>' : '[--reason <text>]'}`);
      const t = await target(connect, argv, verb, id);
      const problem = await t.alerts.setStatus(t.id, {status: STATUS[verb], reason});
      report(problem, `${PAST[verb]}: ${identity(problem, id)}${reason ? ` — ${reason}` : ''}`, output);
      return true;
    }
  }
  printError(new Error(PROBLEMS_USAGE));
  return false;
}

async function target(connect: Connect, argv: any, verb: string, id: string): Promise<{alerts: NodeAlertsClient; id: string}> {
  const alerts = (await connect(singleHost(argv, `problems ${verb}`))).alerts;
  return {alerts, id: await problemId(alerts, id)};
}

/** A problem keeps its identity for good, so `kind:key` is looked up whatever its status. */
export async function problemId(alerts: NodeAlertsClient, id: string): Promise<string> {
  const colon = id.indexOf(':');
  if (colon < 0) return id;
  const matches: any[] = await alerts.problems({kind: id.slice(0, colon), key: id.slice(colon + 1), status: 'all'}) ?? [];
  if (!matches.length)
    throw new Error(`No problem ${id}`);
  if (matches.length > 1)
    throw new Error(`Several problems match ${id}; pass an id:\n` +
      matches.map((p) => `  ${p?.id}  ${p?.kind}:${p?.key}  ${p?.status}  ${truncate(p?.summary, 60)}`).join('\n'));
  return String(matches[0].id);
}

function usage(line: string): boolean {
  printError(new Error(`Usage: grok s o problems ${line}`));
  return false;
}

function optString(v: any): string | undefined {
  return v === undefined || v === null || v === true ? undefined : String(v);
}

function identity(problem: any, fallback: string): string {
  return problem?.kind ? `${problem.kind}:${problem.key}` : fallback;
}

function report(problem: any, line: string, output: OutputFormat): void {
  if (output === 'json' || output === 'csv') printOutput(problem, output);
  else if (output === 'quiet') console.log(problem?.id ?? '');
  else console.log(line);
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
