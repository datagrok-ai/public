import {batch, rawEffect, untracked} from './signals.js';
import {asHelper, captureStack, forgetListener, forgetScope, forgetTimer, registry, trackListener, trackScope,
  trackTimer} from './debug.js';
import {Usage} from './usage.js';
import type {ObservableLike} from './protocol.js';

/** Receives an error a scope caught: from an effect, a bound handler, a listener or a cleanup. */
export type ErrorHandler = (error: unknown, scope: Scope) => void;

/** What the function {@link Scope.task} runs gets: the run's abort signal, and the way back into its owner after an
 * `await`. */
export interface TaskContext {
  /** Aborted when the owner is disposed or the run cancelled: hand it to `fetch`. */
  readonly signal: AbortSignal;
  /** The scope the task runs for. */
  readonly owner: Scope;
  /** Runs `fn` with the owner ambient, so what it builds is the owner's; once the owner is gone or the run
   * cancelled, skips it and answers undefined. */
  build<T>(fn: () => T): T | undefined;
}

/** A run {@link Scope.task} started. */
export interface TaskRun {
  /** Aborts the run's signal; a run still waiting behind a serial one is dropped. */
  cancel(): void;
  /** Settles once the run is over (finished, dropped or cancelled). A throw or a rejection goes to the owner's
   * error handler, never here. */
  readonly done: Promise<void>;
}

/** How {@link Scope.task} runs. */
export interface TaskOptions {
  /** A queue name: one run of the queue in flight on this owner. A run asked for meanwhile waits for it, and a
   * newer one replaces the waiting one — three saves while one runs make that one and the last. */
  serial?: string;
}

interface TaskQueue {
  waiting?: {start: () => void, replace: (by: Promise<void>) => void};
}

const LOOP_RUNS = 25;
const activationHooks = new WeakMap<Scope, () => void>();

/** @internal Runs `hook` once, when the lazy `scope` starts owning something. */
export function whenActive(scope: Scope, hook: () => void): void {
  activationHooks.set(scope, hook);
}

/** Owner of effects, cleanups, listeners, timers and child scopes. Everything registered on a scope
 * is released by a single `dispose()` call. Errors raised in what it owns are reported to
 * {@link onError}, then up the parent chain, then to {@link Scope.onUnhandledError} — never thrown
 * at whoever happened to write a signal. */
export class Scope {
  private static _stack: (Scope | undefined)[] = [];

  /** Where errors no scope handles go. Default: `console.error('u2:', error, '(<label>)')`. */
  static onUnhandledError: ErrorHandler = (error, scope) => console.error('u2:', error, `(${scope.label})`);

  /** The scope this one is disposed with, if any. */
  readonly parent: Scope | undefined;
  /** What created it, for leak reports: `scope`, `effect`, `row`, `popup`, `dialog`, … */
  readonly label: string;
  /** Receives errors raised in what this scope owns; unset, they go to the parent. */
  onError: ErrorHandler | undefined;

  private _disposers?: Set<() => void>;
  private _disposed = false;
  private _counted = false;
  private _unlink?: () => void;
  private _abort?: AbortController;
  private _queues?: Map<string, TaskQueue>;

  /** `lazy`: counted, tracked and linked to `parent` only once it owns something — an effect run
   * that builds nothing costs nothing. */
  constructor(parent?: Scope, options?: {label?: string, lazy?: boolean}) {
    this.parent = parent;
    this.label = options?.label ?? 'scope';
    if (!options?.lazy)
      this._activate();
  }

  /** Scopes created and not yet disposed — a counter, so production holds no references. */
  static get liveCount(): number {
    return registry.live;
  }

  /** The scope that owns whatever is being built right now, if any. Element factories bind
   * signals to it; components created under it are disposed with it. */
  static get ambient(): Scope | undefined {
    return Scope._stack[Scope._stack.length - 1];
  }

  /** Runs `fn` with `scope` ambient; `undefined` runs it with NO owner.
   * @deprecated Use {@link Scope.build} (`owner.build(fn)`); for what must outlive whoever is building, a new tree
   * under {@link Scope.root}. */
  static runWith<T>(scope: Scope | undefined, fn: () => T): T {
    return Scope._with(scope, fn);
  }

  /** A new tree on purpose — what outlives whoever is building (a platform view opened from a row): a scope with no
   * owner, counted and named `label` in leak reports. Whoever holds it disposes it.
   * @example
   * const session = DG.U2.Scope.root('campaign session');
   * session.build(() => new DG.U2.Control());
   * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-task} */
  static root(label: string): Scope {
    return new Scope(undefined, {label});
  }

  private static _with<T>(scope: Scope | undefined, fn: () => T): T {
    Scope._stack.push(scope);
    try {
      return fn();
    }
    finally {
      Scope._stack.pop();
    }
  }

  /** Runs `fn` now with this scope ambient, so what it builds is owned here; synchronous, and a throw propagates.
   * @example
   * const control = owner.build(() => new DG.U2.Control());
   * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-task} */
  build<T>(fn: () => T): T {
    return Scope._with(this, fn);
  }

  /** `x` as a promise when it is one — a native promise, or an object or function with a callable `then` — else
   * undefined. Reading `then` never throws here: a proxy that answers every name (a `DG.Row`) may throw on it, and
   * such an `x` is no promise. The one test of "did a handler answer a promise" for u2 and its hosts. */
  static thenable(x: unknown): PromiseLike<unknown> | undefined {
    if (x instanceof Promise)
      return x;
    if (x === null || typeof x !== 'object' && typeof x !== 'function')
      return undefined;
    try {
      return typeof (x as PromiseLike<unknown>).then === 'function' ? x as PromiseLike<unknown> : undefined;
    }
    catch (_) {
      return undefined;
    }
  }

  /** Whether {@link dispose} ran. */
  get isDisposed(): boolean {
    return this._disposed;
  }

  /** A scope disposed with this one, or earlier on its own. */
  child(label?: string): Scope {
    return new Scope(this, {label});
  }

  /** Runs `fn` now and on every change of the signals it reads. Each run gets its own child scope,
   * ambient while `fn` runs and disposed before the next run and with this scope: what a run builds
   * dies with it. A returned function is that run's cleanup. A throw is reported, not propagated. */
  effect(fn: () => void): this {
    if (this._disposed) {
      this._reportDisposed('effect');
      return this;
    }
    let run: Scope | undefined;
    let runs = 0;
    let warned = false;
    const stack = registry.enabled ? captureStack() : undefined;
    const stop = rawEffect(() => {
      if (run !== undefined) {
        const previous = run;
        run = undefined;
        untracked(() => previous.dispose());
      }
      if (this._disposed)
        return;
      if (registry.enabled && ++runs === 1)
        queueMicrotask(() => runs = 0);
      if (runs === LOOP_RUNS && !warned) {
        warned = true;
        console.warn(`u2: effect re-ran ${LOOP_RUNS} times in one update — a loop? (${this.label})`, stack);
      }
      const current = run = new Scope(this, {label: 'effect', lazy: true});
      try {
        const cleanup: unknown = current.build(fn);
        if (typeof cleanup === 'function')
          current.own(cleanup as () => void);
      }
      catch (e) {
        this.reportError(e);
      }
    });
    return this.own(() => {
      stop();
      run?.dispose();
      run = undefined;
    });
  }

  /** Registers a cleanup to run on {@link dispose}. Runs immediately if already disposed. */
  own(dispose: () => void): this {
    if (!this._disposed)
      this._activate();
    if (this._disposed)
      dispose();
    else
      (this._disposers ??= new Set()).add(dispose);
    return this;
  }

  /** Drops a cleanup registered by {@link own} — for things that come and go many times over
   * the owner's life (an overlay reopened on every click) and would otherwise leave a dead
   * closure behind on every round. */
  disown(dispose: () => void): this {
    this._disposers?.delete(dispose);
    return this;
  }

  /** Wraps a handler the library or the platform calls later: it runs with this scope ambient,
   * inside one `batch`, does nothing once the scope is disposed, and reports a throw instead of
   * propagating it (answering undefined) — and the rejection of a promise it answers. */
  bind<A extends unknown[], R>(fn: (...args: A) => R): (...args: A) => R | undefined {
    return (...args: A) => {
      if (this._disposed)
        return undefined;
      try {
        const result = batch(() => this.build(() => fn(...args)));
        Scope.thenable(result)?.then(undefined, (e: unknown) => this.reportError(e));
        return result;
      }
      catch (e) {
        this.reportError(e);
        return undefined;
      }
    };
  }

  /** Runs `fn` now the way {@link bind} runs a handler. */
  run<T>(fn: () => T): T | undefined {
    return this.bind(fn)();
  }

  /** Subscribes `handler` to `source` — an rxjs observable, a platform event stream, anything with a `subscribe`
   * method — until this scope is disposed. The handler runs as {@link bind} runs one. With u2 debug on the
   * subscription is a live `subscription` scope with its creation stack. Answers the unsubscribe.
   * @example
   * owner.subscribe(grok.events.onViewRemoved, (view) => grok.shell.info(view.name));
   * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-task} */
  subscribe<T>(source: ObservableLike<T>, handler: (value: T) => void): () => void {
    if (this._disposed)
      return this._reportDisposed('subscribe');
    if (typeof (source as Partial<ObservableLike<T>> | null | undefined)?.subscribe !== 'function') {
      Usage.report('source-missing', `owner.subscribe(source, handler) needs an observable; got ${source === null ? 'null' :
        typeof source === 'object' ? 'an object with no subscribe method' : String(source)} — an event this host ` +
        'does not have (an older platform)? Check for it first; nothing was subscribed.');
      return () => {};
    }
    const scope = this.child('subscription');
    // a signal's `subscribe` (plain JS hands one in) answers the unsubscribe itself, not a subscription
    const subscription: {unsubscribe(): void} | (() => void) = source.subscribe(scope.bind(handler));
    scope.own(typeof subscription === 'function' ? subscription : () => subscription.unsubscribe());
    return () => scope.dispose();
  }

  /** Runs async work for this scope: `fn` starts now with the scope ambient, and after each `await` gets back in
   * with `t.build(…)`, which does nothing once the scope is disposed or the run cancelled; `t.signal` aborts then.
   * With `serial`, one run of that queue is in flight: a run asked for meanwhile waits, and a newer one replaces
   * it. A throw or a rejection is reported as {@link bind} reports one.
   * @example
   * owner.task(async (t) => {
   *   const df = await grok.data.getDemoTable('demog.csv');
   *   t.build(() => host.append(ui.divText(`${df.rowCount} rows`)));
   * }, {serial: 'load'});
   * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-task} */
  task<R>(fn: (t: TaskContext) => R | PromiseLike<R>, options?: TaskOptions): TaskRun {
    let settle!: () => void;
    const done = new Promise<void>((resolve) => settle = resolve);
    const controller = new AbortController();
    if (this._disposed) {
      this._reportDisposed('task');
      controller.abort();
      settle();
      return {cancel: () => {}, done};
    }
    const context: TaskContext = {signal: controller.signal, owner: this,
      build: (f) => this._disposed || controller.signal.aborted ? undefined : this.build(f)};
    const key = options?.serial;
    const queue = key === undefined ? undefined : (this._queues ??= new Map()).get(key);
    let waiting: TaskQueue['waiting'];
    const cancel = () => {
      controller.abort();
      this.disown(cancel);
      if (waiting !== undefined && queue?.waiting === waiting) {
        queue.waiting = undefined;
        settle();
      }
    };
    this.own(cancel);
    const start = () => {
      const result = controller.signal.aborted ? undefined : this.bind(() => fn(context))();
      const pending = Scope.thenable(result);
      const over = () => {
        this.disown(cancel);
        settle();
        if (key !== undefined)
          this._nextInQueue(key);
      };
      if (pending === undefined)
        over();
      else
        pending.then(over, over);
    };
    if (key === undefined)
      start();
    else if (queue === undefined) {
      this._queues!.set(key, {});
      start();
    }
    else {
      queue.waiting?.replace(done);
      waiting = queue.waiting = {start, replace: (by) => {
        controller.abort();
        this.disown(cancel);
        void by.then(settle);
      }};
    }
    return {cancel, done};
  }

  private _nextInQueue(key: string): void {
    const queue = this._queues?.get(key);
    const next = queue?.waiting;
    if (next === undefined) {
      this._queues?.delete(key);
      return;
    }
    queue!.waiting = undefined;
    next.start();
  }

  /** Aborted on {@link dispose} — hand it to `fetch`, `addEventListener` or anything cancellable. */
  get abortSignal(): AbortSignal {
    this._abort ??= new AbortController();
    if (this._disposed)
      this._abort.abort();
    return this._abort.signal;
  }

  /** Adds a listener removed on {@link dispose}; the handler runs as {@link bind} runs it. Answers
   * the remover. */
  listen(target: EventTarget, type: string, handler: (e: Event) => void,
    options?: AddEventListenerOptions): () => void {
    if (this._disposed)
      return this._reportDisposed('listen');
    this._activate();
    const bound = this.bind(handler);
    const h = options?.once ? (e: Event) => {
      remove();
      bound(e);
    } : bound;
    const id = trackListener(target, type, options, this);
    const remove = () => {
      target.removeEventListener(type, h, options);
      this.disown(remove);
      forgetListener(id);
    };
    asHelper(() => target.addEventListener(type, h, options));
    this.own(remove);
    return remove;
  }

  /** A timeout cleared on {@link dispose}; `fn` runs as {@link bind} runs it. Answers the canceller. */
  setTimeout(fn: () => void, ms?: number): () => void {
    if (this._disposed)
      return this._reportDisposed('setTimeout');
    this._activate();
    const bound = this.bind(fn);
    const id = trackTimer('timeout', ms, this);
    const handle = globalThis.setTimeout(() => {
      this.disown(cancel);
      forgetTimer(id);
      bound();
    }, ms);
    const cancel = () => {
      globalThis.clearTimeout(handle);
      this.disown(cancel);
      forgetTimer(id);
    };
    this.own(cancel);
    return cancel;
  }

  /** An interval cleared on {@link dispose}; `fn` runs as {@link bind} runs it. Answers the canceller. */
  setInterval(fn: () => void, ms: number): () => void {
    if (this._disposed)
      return this._reportDisposed('setInterval');
    this._activate();
    const bound = this.bind(fn);
    const id = trackTimer('interval', ms, this);
    const handle = asHelper(() => globalThis.setInterval(bound, ms));
    const cancel = () => {
      asHelper(() => globalThis.clearInterval(handle));
      this.disown(cancel);
      forgetTimer(id);
    };
    this.own(cancel);
    return cancel;
  }

  /** An animation frame cancelled on {@link dispose}; `fn` runs as {@link bind} runs it. Answers the
   * canceller. Falls back to a 16 ms timeout where there are no frames (Node). */
  requestFrame(fn: (time: number) => void): () => void {
    if (this._disposed)
      return this._reportDisposed('requestFrame');
    this._activate();
    const bound = this.bind(fn);
    const id = trackTimer('raf', undefined, this);
    const fire = (time: number) => {
      this.disown(cancel);
      forgetTimer(id);
      bound(time);
    };
    const raf = globalThis.requestAnimationFrame;
    const handle = raf === undefined ? globalThis.setTimeout(() => fire(performance.now()), 16) :
      asHelper(() => raf(fire));
    const cancel = () => {
      if (raf === undefined)
        globalThis.clearTimeout(handle as ReturnType<typeof setTimeout>);
      else
        asHelper(() => globalThis.cancelAnimationFrame(handle as number));
      this.disown(cancel);
      forgetTimer(id);
    };
    this.own(cancel);
    return cancel;
  }

  /** Hands an error to {@link onError}, else up the parent chain, else to {@link Scope.onUnhandledError}. */
  reportError(error: unknown): void {
    if (this.onError)
      this.onError(error, this);
    else if (this.parent)
      this.parent.reportError(error);
    else
      Scope.onUnhandledError(error, this);
  }

  /** Runs every cleanup (untracked, each on its own — one throwing does not stop the rest; the
   * errors are reported after all ran), aborts {@link abortSignal} and leaves the parent. */
  dispose(): void {
    if (this._disposed)
      return;
    this._disposed = true;
    this._abort?.abort();
    const disposers = this._disposers;
    this._disposers = undefined;
    let errors: unknown[] | undefined;
    if (disposers !== undefined) {
      untracked(() => {
        for (const d of disposers) {
          try {
            d();
          }
          catch (e) {
            (errors ??= []).push(e);
          }
        }
      });
    }
    if (this._counted) {
      registry.live--;
      forgetScope(this);
    }
    if (this._unlink !== undefined) {
      this.parent?.disown(this._unlink);
      this._unlink = undefined;
    }
    for (const e of errors ?? []) {
      if (this.parent)
        this.parent.reportError(e);
      else
        Scope.onUnhandledError(e, this);
    }
  }

  private _activate(): void {
    if (this._counted)
      return;
    this._counted = true;
    registry.live++;
    trackScope(this, this.label, this.parent);
    if (this.parent !== undefined) {
      const unlink = this._unlink = () => this.dispose();
      this.parent.own(unlink);
    }
    const hook = activationHooks.get(this);
    if (hook !== undefined) {
      activationHooks.delete(this);
      hook();
    }
  }

  private _reportDisposed(what: string): () => void {
    Usage.report('owner-disposed', `${what}() on a scope already disposed (${this.label}) does nothing — a handler ` +
      'that closed its own owner and went on? Register it on a live owner; after an `await`, use ' +
      '`owner.task(async (t) => … t.build(…))`, which skips the work once the owner is gone.');
    return () => {};
  }
}

registry.ambients.push(() => Scope.ambient);
