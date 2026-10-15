/* The automatic audits the harness runs after every passing scenario. On a u2 page, on the
   scenario's final state and before the shell reset: the layout audit (layout.ts) and the console-error
   floor. After the reset: the leak check, which asks `DG.U2.debug` what was created since the baseline
   taken when the scenario started and is still live (scopes, timers, global listeners). Tracking is on
   only in pages whose project audits leaks: an init script sets `__U2_DEBUG__` before the platform
   loads. Modes per project (`audits` in bdd.config.json), env overrides first; `report` warns and
   annotates the test (`bdd-audit`, counted by tool/audit-summary.mjs), `fail` fails the scenario. */
import type {Page, TestInfo} from '@playwright/test';
import {type AuditConfig, type AuditMode, CONFIG_FILE, readConfig} from '../project.js';
import {auditLayout, type LayoutRule} from './layout.js';

export type AuditName = 'layout' | 'errors' | 'leaks';
export type AuditModes = Record<AuditName, AuditMode>;
export type AuditRule = LayoutRule | 'console-error' | 'leak';
export type AuditAllow = Required<NonNullable<AuditConfig['allow']>>;

export interface AuditFinding {audit: AuditName; rule: AuditRule; detail: string}

/** A console error or uncaught exception as the harness collected it: the text and where it came from. */
export interface LoggedError {message: string; url: string}

/** What a scenario's audits did: the audits that ran (a page without u2 runs none) and what they found. */
export interface AuditOutcome {ran: AuditName[]; findings: AuditFinding[]}

/** One annotation per finding; one more naming the audits that ran on the scenario. */
export const AUDIT_ANNOTATION = 'bdd-audit';
export const AUDITS_RAN_ANNOTATION = 'bdd-audits';

const AUDITS: AuditName[] = ['layout', 'errors', 'leaks'];
const MODES: AuditMode[] = ['fail', 'report', 'off'];
const DEFAULT_MODES: AuditModes = {layout: 'report', errors: 'report', leaks: 'off'};

/** Findings every project accepts — u2's own page-wide singletons, platform noise — each with its reason. */
export const BUILTIN_ALLOW: AuditAllow = {layout: [], errors: [], leaks: []};

function modeOf(value: unknown, source: string): AuditMode | undefined {
  if (value === undefined || value === '')
    return undefined;
  if (!MODES.includes(value as AuditMode))
    throw new Error(`${source}: "${value}" is not an audit mode (${MODES.join(', ')})`);
  return value as AuditMode;
}

/** `BDD_AUDIT_LAYOUT|ERRORS|LEAKS` > `BDD_AUDIT` > the project's `audits` > layout and errors `report`, leaks `off`. */
export function resolveModes(config: AuditConfig | undefined, env: Record<string, string | undefined> = process.env): AuditModes {
  const all = modeOf(env.BDD_AUDIT, 'BDD_AUDIT');
  const modes = {} as AuditModes;
  for (const name of AUDITS) {
    const key = `BDD_AUDIT_${name.toUpperCase()}`;
    modes[name] = modeOf(env[key], key) ?? all ?? modeOf(config?.[name], `${CONFIG_FILE} audits.${name}`) ?? DEFAULT_MODES[name];
  }
  return modes;
}

/** The project's allowlists after the library's; an entry without a reason, or with a pattern that is not a
 * regular expression, is refused. */
export function resolveAllow(config: AuditConfig | undefined): AuditAllow {
  const allow: AuditAllow = {
    layout: [...BUILTIN_ALLOW.layout, ...config?.allow?.layout ?? []],
    errors: [...BUILTIN_ALLOW.errors, ...config?.allow?.errors ?? []],
    leaks: [...BUILTIN_ALLOW.leaks, ...config?.allow?.leaks ?? []],
  };
  for (const name of AUDITS) {
    for (const entry of allow[name] as {reason?: string}[]) {
      const where = `${CONFIG_FILE} audits.allow.${name} ${JSON.stringify(entry)}`;
      if (typeof entry.reason !== 'string' || entry.reason.trim() === '')
        throw new Error(`${where}: an allowed finding says why (reason)`);
      for (const key of ['message', 'url', 'label', 'stack'] as const) {
        const source = (entry as Record<string, unknown>)[key];
        if (source === undefined)
          continue;
        try {
          new RegExp(String(source));
        }
        catch (e) {
          throw new Error(`${where}: ${key} is not a regular expression: ${(e as Error).message}`);
        }
      }
    }
  }
  return allow;
}

interface AuditSettings {modes: AuditModes; allow: AuditAllow}
const settings = new Map<string, AuditSettings>();

/** The modes and allowlists of the project at `root`, read once per worker. */
export function auditSettings(root = process.env.BDD_ROOT): AuditSettings {
  const key = root ?? '';
  let found = settings.get(key);
  if (!found) {
    const config = root ? readConfig(root).audits : undefined;
    found = {modes: resolveModes(config), allow: resolveAllow(config)};
    settings.set(key, found);
  }
  return found;
}

export const auditModes = (root = process.env.BDD_ROOT): AuditModes => auditSettings(root).modes;

/** Audits judge a scenario that passed as expected: a failure has its own report, a skip did nothing, and a
 * `@known-failure` scenario stops midway on purpose. */
export function shouldAudit(info: Pick<TestInfo, 'status' | 'expectedStatus' | 'tags'>): boolean {
  return info.status === 'passed' && info.expectedStatus === 'passed' &&
    !info.tags.some((t) => t === '@known-failure' || t.startsWith('@known-failure:'));
}

export function errorFindings(errors: LoggedError[], allow: AuditAllow['errors']): AuditFinding[] {
  return errors
    .filter((e) => !allow.some((a) => new RegExp(a.message).test(e.message) &&
      (a.url === undefined || new RegExp(a.url).test(e.url))))
    .map((e) => ({audit: 'errors', rule: 'console-error',
      detail: `${e.message.split('\n')[0].slice(0, 300)}${e.url ? ` (${e.url})` : ''}`}));
}

/** What the leak check reads back: the records the baseline did not have, singletons left out. */
export interface LeakRead {
  scopes: {id: number; label: string; name?: string; stack?: string}[];
  timers: {kind: string; ms?: number; stack?: string}[];
  listeners: {target: string; type: string; capture: boolean; stack?: string}[];
  singletons: number;
}

/** The first three frames of a creation stack, each `<fn> (<file>:<line>:<col>)` without the URL's path. */
function createdAt(stack: string | undefined): string {
  const frames = (stack ?? '').split('\n').map((l) => l.trim()).filter((l) => l.startsWith('at '))
    .slice(0, 3).map((l) => l.slice(3).replace(/[a-z][\w+.-]*:\/\/[^\s)]*\//gi, ''));
  return frames.length > 0 ? ` — created at ${frames.join(' < ')}` : '';
}

export function leakFindings(read: LeakRead, allow: AuditAllow['leaks']): AuditFinding[] {
  const records = [
    ...read.scopes.map((s) => ({label: s.label, stack: s.stack,
      what: `scope #${s.id} ${s.label}${s.name ? ` (${s.name})` : ''} still live`})),
    ...read.timers.map((t) => ({label: `${t.kind} timer`, stack: t.stack,
      what: `${t.kind}${t.ms !== undefined ? ` ${t.ms} ms` : ''} timer still pending`})),
    ...read.listeners.map((l) => ({label: `${l.target} ${l.type} listener`, stack: l.stack,
      what: `${l.target} ${l.type}${l.capture ? ' capture' : ''} listener still attached`})),
  ];
  return records
    .filter((r) => !allow.some((a) => (a.label !== undefined || a.stack !== undefined) &&
      (a.label === undefined || new RegExp(a.label).test(r.label)) &&
      (a.stack === undefined || new RegExp(a.stack).test(r.stack ?? ''))))
    .map((r) => ({audit: 'leaks', rule: 'leak', detail: `${r.what} after the view closed${createdAt(r.stack)}`}));
}

const tracked = new WeakSet<Page>();

/** On a new page of a project that audits leaks: u2's tracking switched on before the platform loads. */
export async function trackLeaks(page: Page): Promise<void> {
  if (auditSettings().modes.leaks === 'off' || tracked.has(page))
    return;
  tracked.add(page);
  await page.addInitScript(() => { (globalThis as any).__U2_DEBUG__ = true; });
}

let leakCheckUnavailable = false;

/** What is live now, kept in the page as the scenario's leak baseline (a navigation drops it with the
 * document, so a check never compares two documents). A page not in the shell has none. */
export async function markLeakBaseline(page: Page): Promise<void> {
  if (auditSettings().modes.leaks === 'off')
    return;
  const missing: string = await page.evaluate(() => {
    const w = window as any;
    delete w.__bddLeakBaseline;
    if (typeof w.grok?.shell?.closeAll !== 'function')
      return '';
    const debug = w.DG?.U2?.debug;
    if (!debug)
      return 'the platform has no DG.U2.debug (a js-api older than the u2 engine)';
    if (!debug.enabled)
      return 'DG.U2.debug is not tracking (__U2_DEBUG__ was not set before the platform loaded)';
    w.__bddLeakBaseline = debug.snapshot();
    return '';
  }).catch(() => '');
  if (missing && !leakCheckUnavailable) {
    leakCheckUnavailable = true;
    console.warn(`bdd: leaks are audited, but ${missing}: no leak check on this worker`);
  }
}

/** After the shell reset: two frames, the orphan sweep, then what was created since the baseline and is
 * still live, read again every 100 ms for up to 2 s while anything is (a debounce clears itself).
 * `undefined` — no baseline, so no check. */
export async function readLeaks(page: Page): Promise<LeakRead | undefined> {
  const read = await page.evaluate(async () => {
    const w = window as any;
    const since = w.__bddLeakBaseline;
    delete w.__bddLeakBaseline;
    const debug = w.DG?.U2?.debug;
    if (!since || !debug?.enabled)
      return null;
    await new Promise((resolve) => {
      requestAnimationFrame(() => requestAnimationFrame(resolve));
      setTimeout(resolve, 200);
    });
    const deadline = performance.now() + 2000;
    debug.sweep();
    let report = debug.diff(since);
    while (report.total > 0 && performance.now() < deadline) {
      await new Promise((resolve) => setTimeout(resolve, 100));
      debug.sweep();
      report = debug.diff(since);
    }
    const all: {label?: string}[] = [...report.scopes, ...report.timers, ...report.listeners];
    return {
      singletons: all.filter((r) => r.label === 'singleton').length,
      scopes: report.scopes.filter((r: any) => r.label !== 'singleton')
        .map((r: any) => ({id: r.id, label: r.label, name: r.name, stack: r.stack?.split('\n').slice(0, 8).join('\n')})),
      timers: report.timers.filter((r: any) => r.label !== 'singleton')
        .map((r: any) => ({kind: r.kind, ms: r.ms, stack: r.stack?.split('\n').slice(0, 8).join('\n')})),
      listeners: report.listeners.filter((r: any) => r.label !== 'singleton')
        .map((r: any) => ({target: r.target, type: r.type, capture: r.capture,
          stack: r.stack?.split('\n').slice(0, 8).join('\n')})),
    };
  }).catch(() => null);
  return read ?? undefined;
}

/** The leak check of a scenario whose view the shell reset has closed. */
export async function leakCheck(page: Page): Promise<AuditOutcome> {
  const {modes, allow} = auditSettings();
  const read = modes.leaks === 'off' ? undefined : await readLeaks(page);
  return read ? {ran: ['leaks'], findings: leakFindings(read, allow.leaks)} : {ran: [], findings: []};
}

/** The layout audit and the console-error floor on the scenario's final state, when the current view or an
 * overlay over it shows u2. `errors` takes the errors logged since the scenario's floor. */
export async function auditScenario(page: Page, errors: () => LoggedError[]): Promise<AuditOutcome> {
  const {modes, allow} = auditSettings();
  const outcome: AuditOutcome = {ran: [], findings: []};
  if (modes.layout === 'off' && modes.errors === 'off')
    return outcome;
  // one cheap read first: most scenarios of a platform suite show no u2, and the audit waits two frames
  let u2 = await page.evaluate(() => {
    const w = window as any;
    let root: unknown = null;
    try {
      root = w.grok?.shell?.v?.root;
    }
    catch {
      root = null;
    }
    const scope = root instanceof Element && root.isConnected ? root : document.body;
    return scope.querySelector('[data-u2]') !== null ||
      document.querySelector('[data-u2="dialog"], [data-u2-part="popup"], [data-u2="menu"]') !== null;
  }).catch(() => false);
  if (u2 && modes.layout !== 'off') {
    const layout = await auditLayout(page, {allow: allow.layout});
    u2 = layout.u2;
    if (u2) {
      outcome.ran.push('layout');
      outcome.findings.push(...layout.findings.map((f): AuditFinding => ({audit: 'layout', rule: f.rule, detail: f.detail})));
    }
  }
  if (u2 && modes.errors !== 'off') {
    outcome.ran.push('errors');
    outcome.findings.push(...errorFindings(errors(), allow.errors));
  }
  return outcome;
}

const listed = (findings: AuditFinding[]): string => findings.map((f) => `  ${f.rule}: ${f.detail}`).join('\n');

/** Every finding becomes a `bdd-audit` annotation (with the journey scenario it belongs to); the ones in
 * `report` mode are printed, and those in `fail` mode come back as the error that fails the scenario. */
export function settleAudit(outcome: AuditOutcome, title: string, info: Pick<TestInfo, 'annotations'>,
  scenario?: string): Error | undefined {
  const {modes} = auditSettings();
  const where = scenario ? ` (scenario "${scenario}")` : '';
  if (outcome.ran.length > 0)
    info.annotations.push({type: AUDITS_RAN_ANNOTATION, description: `${outcome.ran.join(', ')}${where}`});
  for (const f of outcome.findings)
    info.annotations.push({type: AUDIT_ANNOTATION, description: `${f.rule}: ${f.detail}${where}`});
  const reported = outcome.findings.filter((f) => modes[f.audit] === 'report');
  if (reported.length > 0)
    console.warn(`bdd: audit (report) ${title}: ${reported.length} finding(s)\n${listed(reported)}`);
  const failing = outcome.findings.filter((f) => modes[f.audit] === 'fail');
  return failing.length > 0 ? new Error(`audit: ${failing.length} finding(s) after "${title}"\n${listed(failing)}`) : undefined;
}
