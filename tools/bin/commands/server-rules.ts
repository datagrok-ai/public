/// `grok s observe rules ...` — problem rules: your own problem types, conditions over the log (`/problems/rules`).
import * as fs from 'fs';
import {Connect} from '../utils/server-client';
import {printOutput, printError, OutputFormat} from '../utils/server-output';
import {fmtDateTime, fmtTime, optString, printBlock, rows, truncate} from '../utils/obs-format';

export const RULES_USAGE = `Usage: grok s observe rules <verb> [args]
  list                                   every rule: yours and the deployment's (read-only)
  get <name>                             the rule and its definition
  add [<name>] --json <rule.json|'{...}'>    <name> sets the rule's name
  edit <name> --json <rule.json|'{...}'>     replaces the definition; enabled stays as it is unless given
  enable <name> | disable <name> | delete <name>
  test (<name> | --json <rule.json|'{...}'>) [--hours 24]
                                         what holds now and what would have raised in the last hours
                                         (1 to 24); raises nothing
A rule: {"name", "match": {...}, "groupBy", "window", "when": {...}, "severity", "summary"}; rules the
deployment defines (GROK_PARAMETERS problemRules) are read-only. Needs ManageAlerts.
See https://datagrok.ai/help/govern/audit/problem-rules`;

const PAST: Record<string, string> = {enable: 'enabled', disable: 'disabled'};

/** The rule from `--json`: a file path, or the JSON itself when it starts with `{`; [name] sets its name. */
export function ruleBody(argv: any, name?: string): Record<string, any> {
  const src = optString(argv.json);
  if (!src)
    throw new Error('Pass the rule as --json rule.json or --json \'{...}\'');
  const inline = src.trim().startsWith('{');
  let rule: any;
  try {
    rule = JSON.parse(inline ? src : fs.readFileSync(src, 'utf8'));
  }
  catch (err: any) {
    throw new Error(`Invalid rule JSON in ${inline ? '--json' : src}: ${err.message}`);
  }
  if (rule === null || typeof rule !== 'object' || Array.isArray(rule))
    throw new Error(`Invalid rule JSON in ${inline ? '--json' : src}: a rule is one JSON object`);
  return name ? {...rule, name} : rule;
}

export function ruleRow(r: any): Record<string, any> {
  const d = r?.definition ?? {};
  return {
    NAME: r?.name ?? '',
    SOURCE: r?.source ?? '',
    ON: r?.enabled ? 'yes' : 'no',
    WHEN: Object.keys(d.when ?? {count: 1}).join('+'),
    GROUP: [].concat(d.groupBy ?? []).join(','),
    SEV: d.severity ?? 'warning',
    ONGOING: r?.problems?.ongoing ?? 0,
    ALERTING: r?.problems?.alerting ?? 0,
    ERROR: truncate(r?.error, 60),
  };
}

export function testRows(r: any): {now: Record<string, any>[]; raised: Record<string, any>[]} {
  return {
    now: (r?.now ?? []).map((x: any) => ({GROUP: x?.group ?? x?.key ?? '', SUMMARY: truncate(x?.summary, 80), COUNT: x?.count ?? ''})),
    raised: (r?.raised ?? []).map((x: any) => ({FIRST: fmtTime(x?.firstAt), GROUP: x?.group ?? x?.key ?? '',
      SUMMARY: truncate(x?.summary, 80)})),
  };
}

export async function handleRules(connect: Connect, verb: string | undefined, rest: string[], argv: any,
                                  output: OutputFormat): Promise<boolean> {
  const name = optString(rest[0]);
  const report = (x: any, line: string) => output === 'table' ? console.log(line) : printOutput(x, output);
  switch (verb) {
    case 'list':
      printOutput(rows(await (await connect()).alerts.rules() ?? [], output, ruleRow), output);
      return true;
    case 'get': {
      if (!name) return usage('get <name>');
      const r = await (await connect()).alerts.rule(name);
      if (output === 'table') printRule(r);
      else printOutput(r, output);
      return true;
    }
    case 'add': {
      const body = ruleBody(argv, name);
      const r = await (await connect()).alerts.addRule(body);
      report(r, `added rule ${r?.name ?? name} (${r?.enabled ? 'enabled' : 'disabled'})`);
      return true;
    }
    case 'edit': {
      if (!name) return usage('edit <name> --json <rule.json|\'{...}\'>');
      const body = ruleBody(argv);
      const r = await (await connect()).alerts.editRule(name, body);
      report(r, `replaced the definition of ${name}`);
      return true;
    }
    case 'enable':
    case 'disable': {
      if (!name) return usage(`${verb} <name>`);
      report(await (await connect()).alerts.enableRule(name, verb === 'enable'), `${PAST[verb]} ${name}`);
      return true;
    }
    case 'delete': {
      if (!name) return usage('delete <name>');
      report(await (await connect()).alerts.deleteRule(name), `deleted ${name}`);
      return true;
    }
    case 'test': {
      if (!name && !optString(argv.json)) return usage('test (<name> | --json <rule.json|\'{...}\'>) [--hours 24]');
      const body = optString(argv.json) ? ruleBody(argv, name) : {name};
      const r = await (await connect()).alerts.testRule(body, {hours: argv.hours});
      if (output === 'table') printTest(r, argv.hours ?? 24);
      else printOutput(r, output);
      return true;
    }
  }
  printError(new Error(RULES_USAGE));
  return false;
}

function usage(line: string): boolean {
  printError(new Error(`Usage: grok s observe rules ${line}`));
  return false;
}

function printRule(r: any): void {
  const lines: [string, string][] = [
    ['rule', `${r?.name ?? ''}  ${r?.kind ?? ''}`],
    ['source', `${r?.source ?? ''}${r?.readOnly ? '  (read-only: change it in GROK_PARAMETERS problemRules)' : ''}`],
    ['enabled', `${r?.enabled ? 'yes' : 'no'}${r?.shadowed ? '  (shadowed by the deployment\'s rule)' : ''}`],
    ['problems', `${r?.problems?.ongoing ?? 0} ongoing  · ${r?.problems?.alerting ?? 0} alerting`],
  ];
  if (r?.error) lines.push(['error', r.error]);
  if (r?.createdAt) lines.push(['created', `${fmtDateTime(r.createdAt)}  ${r?.createdBy ?? ''}`]);
  if (r?.updatedAt) lines.push(['updated', `${fmtDateTime(r.updatedAt)}  ${r?.updatedBy ?? ''}`]);
  printBlock(lines);
  console.log('');
  console.log(JSON.stringify(r?.definition ?? {}, null, 2));
}

function printTest(r: any, hours: any): void {
  const t = testRows(r);
  console.log(`Matched ${r?.matched ?? 0} events in ${hours} h (${r?.ms ?? '?'} ms)${r?.truncated ? '  (stopped early: budget)' : ''}`);
  console.log('');
  console.log('Holds now');
  printOutput(t.now.length ? t.now : null, 'table');
  console.log('');
  console.log(`Would have raised in the last ${hours} h`);
  printOutput(t.raised.length ? t.raised : null, 'table');
}
