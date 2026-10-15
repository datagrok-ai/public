import {batch, rawEffect, untracked} from './signals.js';
import {asHelper, captureStack, forgetListener, forgetScope, forgetTimer, registry, trackListener, trackScope,
  trackTimer} from './debug.js';

/** Receives an error a scope caught: from an effect, a bound handler, a listener or a cleanup. */
export type ErrorHandler = (error: unknown, scope: Scope) => void;

const LOOP_RUNS = 25;

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

  /** Runs `fn` with `scope` ambient; `undefined` runs it with NO owner — for what must outlive
   * whoever is building (a platform view opened from a row). */
  static runWith<T>(scope: Scope | undefined, fn: () => T): T {
    Scope._stack.push(scope);
    try {
      return fn();
    }
    finally {
      Scope._stack.pop();
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
        const cleanup: unknown = Scope.runWith(current, fn);
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
        const result = batch(() => Scope.runWith(this, () => fn(...args)));
        if (typeof (result as PromiseLike<unknown> | null | undefined)?.then === 'function')
          (result as PromiseLike<unknown>).then(undefined, (e: unknown) => this.reportError(e));
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
      return () => {};
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
      return () => {};
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
      return () => {};
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
      return () => {};
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
  }
}

registry.ambients.push(() => Scope.ambient);
