import {registry} from './debug.js';

/** How u2 answers a misuse it detects: `off` says nothing, `report` logs it once per code and call site and keeps it
 * in {@link Usage.reports}, `throw` throws a {@link U2UsageError}. */
export type UsageMode = 'off' | 'report' | 'throw';

/** A misuse of u2 the engine detected. `code` names the rule (`owner-missing`, `owner-disposed`, …); the message
 * names the fix. */
export class U2UsageError extends Error {
  /** The rule broken, as {@link UsageReport.code}. */
  readonly code: string;

  constructor(code: string, message: string) {
    super(`u2 usage: ${code}: ${message}`);
    this.name = 'U2UsageError';
    this.code = code;
  }
}

/** One misuse as {@link Usage.reports} lists it. */
export interface UsageReport {
  /** The rule broken: `owner-missing`, `owner-disposed`, `reserved-key`, `source-missing`. */
  code: string;
  /** What happened and how to fix it. */
  message: string;
  /** Where to fix it: the first stack frame outside the engine and u2 — in a bundle, where u2's frames cannot be told
   * apart, the first three outside the engine. */
  site: string;
  /** How many times it happened at that site. */
  count: number;
}

interface UsageState {
  mode?: UsageMode;
  setting?: {value: unknown};
  reports: Map<string, UsageReport>;
}

const g = globalThis as Record<PropertyKey, any>;
const state: UsageState = g[Symbol.for('u2.usage')] ??= {reports: new Map()};

const ENGINE = /u2core[\\/]|[\\/]js-api[\\/](dist|src)[\\/]|[\\/]js-api(\.min)?\.js|datagrok-api[\\/]|[\\/]datagrok\.js/;
const LIBRARY = /[\\/]libraries[\\/]u2[\\/](src|dist|\.test-build)[\\/]|@datagrok-libraries[\\/]u2[\\/]/;
const INTERNAL = /\(node:|^at node:|<anonymous>\)?$/;

/** What the engine does about a misuse of u2 — the one place every usage rule reports through. The mode is
 * `report` unless `__U2_STRICT__`, `localStorage['u2.strict']` or `U2_STRICT` says `off`, `report` or `throw` (`1`);
 * with u2 debug on in the browser (`localStorage['u2.debug'] = '1'`, dev) it is `throw`. The test runners set
 * `throw`: u2's node runner, bdd and `grok test`. A report is a console warning; in Node a line on stderr, and
 * `U2_USAGE_REPORT=<file>` appends each new one to that file as a JSON line.
 * @example
 * DG.U2.Usage.mode = 'throw';
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-usage} */
export class Usage {
  /** The mode in force: set here, else read from the environment, else `report`. */
  static get mode(): UsageMode {
    return state.mode ??= Usage._modeOf(Usage._setting()) ?? 'report';
  }

  static set mode(x: UsageMode) {
    state.mode = x;
  }

  /** Strict mode: Node, a page whose `__U2_STRICT__` or `u2.strict` is set, or u2 debug tracking on (dev, bdd).
   * Checks that cost a walk per class run only there. */
  static get strict(): boolean {
    const setting = Usage._setting();
    return g.process?.versions?.node !== undefined || registry.enabled ||
      setting !== undefined && Usage._modeOf(setting) !== 'off';
  }

  /** Reports the misuse `code`: in `throw` mode throws a {@link U2UsageError}, in `report` mode logs it the first
   * time it happens at a call site, in `off` mode does nothing. */
  static report(code: string, message: string): void {
    const mode = Usage.mode;
    if (mode === 'off')
      return;
    const fresh = Usage._record(code, message);
    if (mode === 'throw')
      throw new U2UsageError(code, message);
    if (fresh !== undefined)
      Usage._log(fresh);
  }

  /** Throws a {@link U2UsageError} whatever the mode — for a misuse u2 cannot go on from — and records it. */
  static fail(code: string, message: string): never {
    Usage._record(code, message);
    throw new U2UsageError(code, message);
  }

  /** Every misuse reported on this page so far, one entry per code, message and site. */
  static reports(): UsageReport[] {
    return [...state.reports.values()].map((r) => ({...r}));
  }

  /** Forgets the reports, so each misuse is logged again. */
  static clear(): void {
    state.reports.clear();
  }

  private static _record(code: string, message: string): UsageReport | undefined {
    const site = Usage._site();
    const key = `${code}\n${message}\n${site}`;
    const known = state.reports.get(key);
    if (known !== undefined) {
      known.count++;
      return undefined;
    }
    const report = {code, message, site, count: 1};
    state.reports.set(key, report);
    return report;
  }

  private static _log(report: UsageReport): void {
    const text = `u2 usage: ${report.code}: ${report.message}${report.site ? ` (at ${report.site})` : ''}`;
    const p = g.process;
    if (p?.versions?.node === undefined) {
      console.warn(text);
      return;
    }
    p.stderr?.write?.(`${text}\n`);
    const file = p.env?.U2_USAGE_REPORT;
    if (file)
      p.getBuiltinModule?.('node:fs')?.appendFileSync(file, JSON.stringify({...report, test: p.argv?.[1]}) + '\n');
  }

  private static _site(): string {
    const E = Error as {stackTraceLimit?: number};
    const limit = E.stackTraceLimit;
    E.stackTraceLimit = 40;
    const frames = (new Error().stack ?? '').split('\n').slice(1).map((l) => l.trim())
      .filter((l) => l !== '' && !ENGINE.test(l) && !INTERNAL.test(l));
    E.stackTraceLimit = limit;
    // where u2's frames can be told apart (Node), the first frame outside it; in a bundle, the first three
    const outside = frames.filter((l) => !LIBRARY.test(l));
    const picked = outside.length < frames.length ? outside.slice(0, 1) : frames.slice(0, 3);
    return picked.map((frame) => frame.replace(/^at /, '').replace(/file:\/\//g, '')
      .replace(/[^\s(]*[\\/]public[\\/]/g, '')).join(' < ');
  }

  private static _setting(): unknown {
    return (state.setting ??= {value: Usage._read()}).value;
  }

  private static _read(): unknown {
    if (g.__U2_STRICT__ !== undefined)
      return g.__U2_STRICT__;
    try {
      const stored = g.localStorage?.getItem('u2.strict');
      if (stored != null)
        return stored;
      if (g.localStorage?.getItem('u2.debug') === '1')
        return 'throw';
    }
    catch (_) {
    }
    return g.process?.env?.U2_STRICT;
  }

  private static _modeOf(x: unknown): UsageMode | undefined {
    if (x === true || x === 1 || x === '1' || x === 'true' || x === 'throw')
      return 'throw';
    if (x === 'report')
      return 'report';
    return x === false || x === 0 || x === '0' || x === 'false' || x === 'off' ? 'off' : undefined;
  }
}
