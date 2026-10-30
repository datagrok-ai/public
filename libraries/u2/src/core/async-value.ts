import {signal, computed, untracked, Signal, ReadonlySignal} from './signals.js';
import {Scope} from './scope.js';
import {Ownership} from './owner.js';
import type {Owner} from './owner.js';

export type AsyncValueState<T> =
  {kind: 'idle'} |
  {kind: 'loading'} |
  {kind: 'ready', value: T} |
  {kind: 'error', message: string};

export interface AsyncValueOptions {
  /** Signals read here re-trigger the run (debounced) — how a bound param makes a source reactive. */
  deps?: () => void;
  /** Coalesces the runs `deps` trigger, the first one included (300ms by default). A value without `deps` loads
   * once, at once: there is nothing to coalesce. */
  debounceMs?: number;
  /** Run at construction and on deps changes (default true); false = only {@link AsyncValue.refresh} runs. */
  auto?: boolean;
  /** Disposes the value with it: whoever is building by default; after an `await` nobody is — pass it. */
  owner?: Owner;
}

/** The single-value twin of {@link AsyncSource} (DD8/C2): one abortable run at a time, stale
 * responses dropped, a re-trigger while loading cancels and restarts, dep changes debounced.
 * Consumers never see this object — a source projects its outputs as computeds over `state`. */
export class AsyncValue<T> {
  readonly state: ReadonlySignal<AsyncValueState<T>>;
  /** The loaded value's signal, for wiring; undefined until a run is ready. */
  readonly signal: ReadonlySignal<T | undefined>;
  private readonly _state: Signal<AsyncValueState<T>>;
  private readonly _fetch: (abort: AbortSignal) => Promise<T>;
  private readonly _debounceMs: number;
  private readonly _scope = new Scope();
  private readonly _release: (() => void) | undefined;
  private _timer: ReturnType<typeof setTimeout> | undefined;
  private _abort: AbortController | undefined;
  private _disposed = false;

  constructor(fetch: (abort: AbortSignal) => Promise<T>, options?: AsyncValueOptions) {
    this._fetch = fetch;
    this._debounceMs = options?.debounceMs ?? 300;
    this._state = signal<AsyncValueState<T>>({kind: 'idle'});
    this.state = this._state;
    this.signal = computed(() => {
      const s = this._state.value;
      return s.kind === 'ready' ? s.value : undefined;
    });
    this._release = Ownership.adopt(options?.owner, () => this.dispose());
    if (options?.auto ?? true) {
      const deps = options?.deps;
      if (deps) {
        this._scope.effect(() => {
          deps();
          this._trigger();
        });
      } else
        untracked(() => void this.refresh());
    }
  }

  /** The loaded value, followed when read inside an `effect` or a `computed`; undefined until a run is ready. */
  get value(): T | undefined {
    return this.signal.value;
  }

  /** The loaded value, not followed. */
  peek(): T | undefined {
    return this.signal.peek();
  }

  /** Manual run; retry after an error is the same call. It starts NOW: the debounce exists to
   * coalesce dependency changes, and a user who asked for a run must not wait on a timer that
   * makes the ask look ignored. Answers when the run has landed, so the caller can report it. */
  refresh(): Promise<void> {
    this._cancel();
    if (this._disposed)
      return Promise.resolve();
    this._state.value = {kind: 'loading'};
    return this._run();
  }

  dispose(): void {
    this._disposed = true;
    this._cancel();
    this._scope.dispose();
    this._release?.();
  }

  private _trigger(): void {
    this._cancel();
    if (this._disposed)
      return;
    this._state.value = {kind: 'loading'};
    this._timer = setTimeout(() => this._run(), this._debounceMs);
  }

  private _cancel(): void {
    clearTimeout(this._timer);
    this._timer = undefined;
    if (this._abort)
      this._abort.abort();
    this._abort = undefined;
  }

  private async _run(): Promise<void> {
    const abort = new AbortController();
    this._abort = abort;
    try {
      const value = await this._fetch(abort.signal);
      if (abort.signal.aborted || this._disposed)
        return;
      this._state.value = {kind: 'ready', value};
    } catch (e) {
      if (abort.signal.aborted || this._disposed)
        return;
      this._state.value = {kind: 'error', message: e instanceof Error ? e.message : String(e)};
    }
  }
}

export function asyncValue<T>(fetch: (abort: AbortSignal) => Promise<T>,
  options?: AsyncValueOptions): AsyncValue<T> {
  return new AsyncValue<T>(fetch, options);
}
