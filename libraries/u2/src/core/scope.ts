import {rawEffect} from './signals.js';

/** Owner of effects and cleanups. Everything registered on a scope is released by a single
 * `dispose()` call; the live-scope count feeds leak detection. */
export class Scope {
  private static _live = new Set<Scope>();
  private static _stack: Scope[] = [];
  private _disposers: (() => void)[] = [];
  private _disposed = false;

  constructor() {
    Scope._live.add(this);
  }

  /** Number of scopes created and not yet disposed — the leak detector's ground truth. */
  static get liveCount(): number {
    return Scope._live.size;
  }

  /** The scope that owns whatever is being built right now, if any. Element factories bind
   * signals to it; components created under it are disposed with it. */
  static get ambient(): Scope | undefined {
    return Scope._stack[Scope._stack.length - 1];
  }

  static runWith<T>(scope: Scope, fn: () => T): T {
    Scope._stack.push(scope);
    try {
      return fn();
    } finally {
      Scope._stack.pop();
    }
  }

  get isDisposed(): boolean {
    return this._disposed;
  }

  /** Runs `fn` now and on every change of the signals it reads. Disposed with the scope. */
  effect(fn: () => void): this {
    return this.own(rawEffect(fn));
  }

  /** Registers a cleanup to run on {@link dispose}. Runs immediately if already disposed. */
  own(dispose: () => void): this {
    if (this._disposed)
      dispose();
    else
      this._disposers.push(dispose);
    return this;
  }

  dispose(): void {
    if (this._disposed)
      return;
    this._disposed = true;
    Scope._live.delete(this);
    for (const d of this._disposers)
      d();
    this._disposers.length = 0;
  }
}
