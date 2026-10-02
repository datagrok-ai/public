/// `grok s o capture ...` (capture rules, LoggingRouter `/logging/capture`) and `grok s o timeline`
/// (one action, request, session, report or rule in time order, `/log/timeline`).
import {Query} from '../utils/node-observability';
import {Connect, singleHost} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {fmtClock, fmtDateTime, fmtSpan, hasValue, listArg, normalizeFlag, normalizeLevel, parseDuration, parseTime,
  printBlock, shortRequestId, sinceArg, truncate} from '../utils/obs-format';

export const CAPTURE_USAGE = `Usage: grok s o capture <verb> [args]
  add (--user <login> | --group <name> | --package <name> | --everyone)
      [--view <name> | --element "<path>" | --function <nqName> | --error <signature>]
      --capture <items> (--for 2d | --until <iso>) --reason <text>
      [--limit <max events>] [--max-sessions n] [--window 10m] [--anonymous] [--name <n>]
  list [--all] [--since 90d]
  show <cap-N|id> [--timeline] [--output csv|json]
  stop <cap-N|id> [--reason <text>]
Capture items: clicks, inputs, requests, calls, errors, and last server:<level>=<flag>[,<flag>...]
  (--capture clicks,inputs,requests,calls,errors,server:debug=queries,files)
--anonymous applies to --group and --everyone only; --everyone needs a scope.`;

export const TIMELINE_USAGE = `Usage: grok s o timeline (--action <id> | --request <id> | --session <id> | --report <n> | --rule <cap-N>)
                       [--from <iso|-1h>] [--to <iso>] [--limit n]`;

const CAPTURE_ITEMS = ['clicks', 'inputs', 'requests', 'calls', 'errors'];
const SUBJECTS = ['user', 'group', 'package', 'everyone'];
const RULE_SCOPES = ['view', 'element', 'function', 'error'];
const TIMELINE_KEYS = ['action', 'request', 'session', 'report', 'rule'];

export interface CaptureSpec {
  clicks: boolean; inputs: boolean; requests: boolean; calls: boolean; errors: boolean;
  serverLevel: string | null; debugFlags: string[];
}

/** `clicks,requests,server:debug=queries,files`: everything after `server:<level>=` is debug flags. */
export function parseCapture(value: any): CaptureSpec {
  const raw = value === undefined || value === true ? '' : String(value);
  const at = raw.indexOf('server:');
  const plain = listArg(at < 0 ? raw : raw.slice(0, at));
  const spec: CaptureSpec = {clicks: false, inputs: false, requests: false, calls: false, errors: false, serverLevel: null, debugFlags: []};
  if (!plain.length && at < 0)
    throw new Error(`--capture needs items: ${CAPTURE_ITEMS.join(', ')}, server:<level>=<flags>`);
  for (const item of plain) {
    if (!CAPTURE_ITEMS.includes(item))
      throw new Error(`Unknown capture item '${item}'. Valid: ${CAPTURE_ITEMS.join(', ')}, server:<level>=<flags> (last)`);
    (spec as any)[item] = true;
  }
  if (at >= 0) {
    const server = raw.slice(at + 'server:'.length);
    const eq = server.indexOf('=');
    spec.serverLevel = normalizeLevel(eq < 0 ? server : server.slice(0, eq));
    spec.debugFlags = eq < 0 ? [] : listArg(server.slice(eq + 1)).map(normalizeFlag);
    if (spec.debugFlags.includes('credentials'))
      throw new Error('A capture rule never turns on the credentials flag');
  }
  return spec;
}

function one(argv: any, keys: string[], what: string, required: boolean): string | undefined {
  const given = keys.filter((k) => argv[k] !== undefined && argv[k] !== false);
  if (given.length > 1 || (required && !given.length))
    throw new Error(`capture add takes ${required ? 'exactly' : 'at most'} one ${what}: ${keys.map((k) => `--${k}`).join(', ')}`);
  return given[0];
}

/** The POST /logging/capture body; every refusal the server makes about the shape is made here first. */
export function captureBody(argv: any, now: Date = new Date()): Record<string, any> {
  const subject = one(argv, SUBJECTS, 'subject', true)!;
  const scope = one(argv, RULE_SCOPES, 'scope', false);
  const reason = argv.reason === undefined || argv.reason === true ? undefined : String(argv.reason);
  if (!reason) throw new Error('A capture rule needs --reason <text>');
  if (argv.for === undefined && argv.until === undefined) throw new Error('A capture rule needs --for <duration> or --until <iso>');
  if (argv.for !== undefined && argv.until !== undefined) throw new Error('Use either --for or --until');
  if (argv.anonymous === true && subject !== 'group' && subject !== 'everyone')
    throw new Error('--anonymous applies to --group and --everyone rules only');
  if (subject === 'everyone' && !scope) throw new Error('An --everyone rule needs a scope: --view, --element, --function or --error');
  if (subject !== 'everyone' && (argv[subject] === true || argv[subject] === ''))
    throw new Error(`--${subject} needs a value`);
  const body: Record<string, any> = {
    name: argv.name === undefined ? undefined : String(argv.name),
    subject: {type: subject, value: subject === 'everyone' ? undefined : String(argv[subject])},
    scope: scope ? {type: scope, value: String(argv[scope])} : undefined,
    capture: parseCapture(argv.capture),
    anonymous: argv.anonymous === true,
    reason,
  };
  if (argv.window !== undefined) body.windowMinutes = parseDuration(argv.window, '--window') / 60000;
  if (argv['max-sessions'] !== undefined) body.maxSessions = Number(argv['max-sessions']);
  if (argv.limit !== undefined) body.maxEvents = Number(argv.limit);
  if (argv.for !== undefined) body.forMinutes = parseDuration(argv.for, '--for') / 60000;
  else {
    const until = parseTime(argv.until, '--until', now);
    if (until.getTime() <= now.getTime()) throw new Error(`--until ${argv.until} is in the past`);
    body.expiresAt = until.toISOString();
  }
  return body;
}

const ruleId = (r: any) => r?.number !== undefined ? `cap-${r.number}` : String(r?.id ?? '');

function subjectText(r: any): string {
  const s = r?.subject ?? {};
  return s.type === 'everyone' ? 'everyone' : `${s.type ?? ''} ${s.value ?? ''}`.trim();
}

function scopeText(r: any): string {
  if (r?.scope?.type) return `${r.scope.type} ${r.scope.value ?? ''}`.trim();
  return r?.capture?.serverLevel ? `server ${r.capture.serverLevel}` : 'all activity';
}

/** How long the rule runs (or ran), with how it ended. */
function activeText(r: any): string {
  const start = Date.parse(r?.createdAt);
  const end = Date.parse(r?.status === 'active' ? r?.expiresAt : (r?.endedAt ?? r?.expiresAt));
  const span = isNaN(start) || isNaN(end) ? '' : fmtSpan(end - start);
  return r?.status && r.status !== 'active' ? `${span} (${r.status})`.trim() : span;
}

export function ruleRow(r: any): Record<string, any> {
  return {
    RULE: ruleId(r),
    AUTHOR: r?.author ?? '',
    SUBJECT: subjectText(r),
    SCOPE: scopeText(r),
    REASON: truncate(r?.reason, 40),
    ACTIVE: activeText(r),
    EVENTS: r?.events ?? 0,
  };
}

/** `rule cap-17  active until 2026-09-30 10:20 · 1 user · 1 view · 0/2000 events` */
export function ruleSummary(r: any): string {
  const who = r?.subject?.type === 'everyone' ? 'everyone' : `1 ${r?.subject?.type ?? 'subject'}`;
  const where = r?.scope?.type ? `1 ${r.scope.type}` : 'all activity';
  const status = r?.status && r.status !== 'active' ? r.status : 'active until';
  const until = status === 'active until' ? ` ${fmtDateTime(r?.expiresAt)}` : '';
  return `rule ${ruleId(r)}  ${status}${until} · ${who} · ${where} · ${r?.events ?? 0}/${r?.maxEvents ?? '?'} events`;
}

function captureItems(c: any): string {
  const items = CAPTURE_ITEMS.filter((i) => c?.[i]);
  if (c?.serverLevel) items.push(`server:${c.serverLevel}${c.debugFlags?.length ? `=${c.debugFlags.join(',')}` : ''}`);
  return items.join(',');
}

export function timelineRow(e: any): Record<string, any> {
  return {
    TIME: fmtClock(e?.time),
    SOURCE: e?.source ?? '',
    SERVER: e?.server ?? '',
    KIND: e?.kind ?? '',
    SUMMARY: truncate(e?.summary, 70),
    STATUS: e?.status ?? '',
    MS: e?.ms ?? '',
    REQ: shortRequestId(e?.requestId),
  };
}

function printTimeline(events: any[], output: OutputFormat): void {
  printOutput(output === 'json' ? events : events.map(timelineRow), output);
}

export async function handleCapture(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                    output: OutputFormat): Promise<boolean> {
  const id = rest[0] === undefined ? undefined : String(rest[0]);
  const logging = async () => (await connect(singleHost(argv, `capture ${verb}`))).logging;
  switch (verb) {
    case 'add': {
      const rule = await (await logging()).addCaptureRule(captureBody(argv));
      if (output === 'table') console.log(ruleSummary(rule));
      else if (output === 'quiet') console.log(ruleId(rule));
      else printOutput(rule, output);
      return true;
    }
    case 'list': {
      const since = argv.since === undefined ? undefined : sinceArg(argv.since);
      const rules: any[] = await (await logging()).captureRules({all: argv.all === true ? true : undefined, since}) ?? [];
      if (output === 'quiet') for (const r of rules) console.log(ruleId(r));
      else printOutput(output === 'json' ? rules : rules.map(ruleRow), output);
      return true;
    }
    case 'show': {
      if (!id) return usage('show <cap-N|id> [--timeline] [--output csv|json]');
      const client = await logging();
      if (argv.timeline === true) {
        printTimeline(await client.timeline({rule: id, limit: argv.limit}) ?? [], output);
        return true;
      }
      const rule = await client.captureRule(id);
      if (output !== 'table') { printOutput(rule, output); return true; }
      printBlock([
        ['rule', `${ruleId(rule)}${rule?.name ? `  ${rule.name}` : ''}`],
        ['author', rule?.author ?? ''],
        ['subject', `${subjectText(rule)}${rule?.anonymous ? '  (anonymous)' : ''}`],
        ['scope', scopeText(rule)],
        ['capture', captureItems(rule?.capture)],
        ['reason', rule?.reason ?? ''],
        ['status', `${rule?.status ?? ''}  ${activeText(rule)}  · created ${fmtDateTime(rule?.createdAt)}` +
          ` · expires ${fmtDateTime(rule?.expiresAt)}`],
        ['events', `${rule?.events ?? 0}/${rule?.maxEvents ?? '?'}  · window ${rule?.windowMinutes ?? '?'} min` +
          ` · max sessions ${rule?.maxSessions ?? '?'}`],
      ]);
      const activations: any[] = Array.isArray(rule?.activations) ? rule.activations : [];
      if (activations.length) {
        console.log('');
        printOutput(activations.map((a) => ({
          SESSION: String(a?.sessionId ?? '').slice(0, 8), USER: a?.user ?? '', ACTIVATED: fmtDateTime(a?.activatedAt),
          UNTIL: fmtDateTime(a?.until), TRIGGER: truncate(a?.triggerDetail, 50),
        })), 'table');
      }
      return true;
    }
    case 'stop': {
      if (!id) return usage('stop <cap-N|id> [--reason <text>]');
      const reason = argv.reason === undefined ? undefined : String(argv.reason);
      const rule = await (await logging()).stopCaptureRule(id, reason);
      if (output === 'table') console.log(`stopped ${ruleId(rule) || id}${reason ? ` — ${reason}` : ''}`);
      else printOutput(rule, output);
      return true;
    }
  }
  printError(new Error(CAPTURE_USAGE));
  return false;
}

function usage(line: string): boolean {
  printError(new Error(`Usage: grok s o capture ${line}`));
  return false;
}

export function timelineQuery(argv: any, now: Date = new Date()): Query {
  const keys = TIMELINE_KEYS.filter((k) => hasValue(argv[k]));
  if (keys.length !== 1)
    throw new Error(`timeline takes exactly one of ${TIMELINE_KEYS.map((k) => `--${k}`).join(', ')}`);
  const q: Query = {[keys[0]]: String(argv[keys[0]]), limit: argv.limit};
  if (argv.from !== undefined) q.from = parseTime(argv.from, '--from', now).toISOString();
  if (argv.to !== undefined) q.to = parseTime(argv.to, '--to', now).toISOString();
  return q;
}

export async function handleTimeline(connect: Connect, _verb: string | undefined, _rest: string[], argv: any,
                                     output: OutputFormat): Promise<boolean> {
  let q: Query;
  try {
    q = timelineQuery(argv);
  }
  catch (err: any) {
    printError(new Error(`${err.message}\n${TIMELINE_USAGE}`));
    return false;
  }
  const dapi = await connect(singleHost(argv, 'timeline'));
  printTimeline(await dapi.logging.timeline(q) ?? [], output);
  return true;
}
