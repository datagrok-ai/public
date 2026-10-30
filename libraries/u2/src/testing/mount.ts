import {Scope, Usage, debug} from 'datagrok-api/u2core';
import type {DebugSnapshot, LeakReport} from 'datagrok-api/u2core';

/** What {@link mount} attaches: a control (its root) or an element. */
export type Mountable = HTMLElement | {root: HTMLElement};

/** Options of {@link mount}. */
export interface MountOptions {
  /** Where the root goes: `document.body` by default. */
  into?: Element;
}

interface Leaks {
  scopes: number;
  report: LeakReport;
}

const SETTLE_MS = 2000;

/** A test body as any runner takes it (`node:test`, `@datagrok-libraries/test`): `body` runs under a root scope of
 * the test's own, so what it builds is owned and goes when it ends, as does what {@link mount} attached. The test
 * fails when, once its scope is gone, a scope, a timer or a global listener it made is still live; when a u2 usage
 * error is reported while it runs; or when an error reaches `Scope.onUnhandledError` (a handler, an effect or a
 * cleanup that threw).
 * @example
 * test('the pencil edits the name', u2test(() => {
 *   const name = mount(() => new TextInput({label: 'Name', value: 'KHL-3', readOnly: true}));
 *   user.click(byName(name, 'Edit Name'));
 * })); */
export function u2test(body: () => unknown): () => Promise<void> {
  return () => TestRun.run(body);
}

/** Builds `build`'s control (or element) under the running {@link u2test}'s scope and attaches its root to the page;
 * it is disposed and detached when the test ends.
 * @example
 * const input = mount(() => new TextInput({label: 'Name'})); */
export function mount<T extends Mountable>(build: () => T, options?: MountOptions): T {
  const run = TestRun.current('mount');
  const built = run.root.build(build);
  const el = built instanceof HTMLElement ? built : built.root;
  (options?.into ?? document.body).append(el);
  run.mounted.push(el);
  return built;
}

class TestRun {
  private static _current: TestRun | undefined;
  readonly root = Scope.root('u2test');
  readonly mounted: Element[] = [];
  readonly errors: unknown[] = [];

  static current(what: string): TestRun {
    if (TestRun._current === undefined)
      throw new Error(`${what} runs inside a u2test body: test(name, u2test(() => …))`);
    return TestRun._current;
  }

  static async run(body: () => unknown): Promise<void> {
    if (TestRun._current !== undefined)
      throw new Error('u2test: another u2test body is running; u2 tests run one at a time');
    const live = Scope.liveCount;
    const tracked = debug.enabled;
    if (!tracked)
      debug.enable({orphanSeconds: 0});
    const since = debug.snapshot();
    const reported = TestRun._usage();
    const unhandled = Scope.onUnhandledError;
    const run = TestRun._current = new TestRun();
    Scope.onUnhandledError = (error, scope) => {
      run.errors.push(error);
      unhandled(error, scope);
    };
    let leaks: Leaks;
    try {
      try {
        await run.root.build(body);
      } finally {
        TestRun._current = undefined;
        run.root.dispose();
        for (const el of run.mounted)
          el.remove();
        const active = document.activeElement as HTMLElement | null;
        if (active !== null && !active.isConnected)
          active.blur();
        Scope.onUnhandledError = unhandled;
      }
      leaks = await TestRun._settled(live, since);
    } finally {
      if (!tracked)
        debug.disable();
    }
    const problems = [...TestRun._leaks(leaks), ...TestRun._reports(reported), ...run.errors.map((e) =>
      `an error reached Scope.onUnhandledError (a handler, an effect or a cleanup threw): ${TestRun._message(e)}`)];
    if (problems.length > 0)
      throw new Error(`u2test: ${problems.join('\n  ')}`, {cause: run.errors[0]});
  }

  private static async _settled(live: number, since: DebugSnapshot): Promise<Leaks> {
    const deadline = Date.now() + SETTLE_MS;
    for (;;) {
      await new Promise((resolve) => setTimeout(resolve, 1));
      const found = {scopes: Scope.liveCount - live, report: debug.diff(since)};
      if (found.scopes <= 0 && found.report.total === 0 || Date.now() > deadline)
        return found;
      await new Promise((resolve) => setTimeout(resolve, 50));
    }
  }

  private static _leaks({scopes, report}: Leaks): string[] {
    const counted = <T extends {label?: string}>(records: T[]) => records.filter((r) => r.label !== 'singleton');
    const [s, t, l] = [counted(report.scopes), counted(report.timers), counted(report.listeners)];
    if (scopes <= 0 && s.length + t.length + l.length === 0)
      return [];
    const items = [
      ...s.map((r) => `scope ${r.label}${r.name ? ` "${r.name}"` : ''}${TestRun._at(r.stack)}`),
      ...t.map((r) => `${r.kind}${r.ms === undefined ? '' : ` ${r.ms} ms`}${TestRun._at(r.stack)}`),
      ...l.map((r) => `${r.target} '${r.type}' listener${TestRun._at(r.stack)}`),
    ];
    return [`the test left ${Math.max(scopes, s.length)} scope(s), ${t.length} timer(s) and ${l.length} global ` +
      `listener(s) live once its own scope was disposed${items.length > 0 ? `: ${items.join('; ')}` : ''}`];
  }

  private static _at(stack: string | undefined): string {
    const frame = stack?.split('\n').map((l) => l.trim())
      .find((l) => l.startsWith('at ') && !/[\\/]testing[\\/]/.test(l));
    return frame === undefined ? '' : ` (${frame})`;
  }

  private static _usage(): Map<string, number> {
    return new Map(Usage.reports().map((r) => [`${r.code}\n${r.message}\n${r.site}`, r.count]));
  }

  private static _reports(before: Map<string, number>): string[] {
    return Usage.reports().filter((r) => r.count > (before.get(`${r.code}\n${r.message}\n${r.site}`) ?? 0))
      .map((r) => `u2 usage: ${r.code}: ${r.message}${r.site ? ` (at ${r.site})` : ''}`);
  }

  private static _message(e: unknown): string {
    return e instanceof Error ? e.message : String(e);
  }
}
