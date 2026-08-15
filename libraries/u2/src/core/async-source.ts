import {signal, Signal, ReadonlySignal} from './signals.js';

export type AsyncState<T> =
  {kind: 'idle'} |
  {kind: 'loading'} |
  {kind: 'ready', items: T[]} |
  {kind: 'empty'} |
  {kind: 'error', message: string};

export type AsyncFetch<T> = (query: string, abort: AbortSignal) => Promise<T[]>;

/** The async-data contract every u2 async control renders against: debounce, per-keystroke
 * cancellation, and the standardized idle/loading/ready/empty/error-with-retry states —
 * written once so no control hand-rolls spinners or races. */
export class AsyncSource<T> {
  readonly state: ReadonlySignal<AsyncState<T>>;
  private readonly _state: Signal<AsyncState<T>>;
  private readonly _fetch: AsyncFetch<T>;
  private readonly _debounceMs: number;
  private _query = '';
  private _timer: ReturnType<typeof setTimeout> | undefined;
  private _abort: AbortController | undefined;
  private _disposed = false;

  constructor(fetch: AsyncFetch<T>, options?: {debounceMs?: number}) {
    this._fetch = fetch;
    this._debounceMs = options?.debounceMs ?? 150;
    this._state = signal<AsyncState<T>>({kind: 'idle'});
    this.state = this._state;
  }

  query(q: string): void {
    this._query = q;
    this._cancel();
    if (this._disposed)
      return;
    this._state.value = {kind: 'loading'};
    this._timer = setTimeout(() => this._run(), this._debounceMs);
  }

  retry(): void {
    this.query(this._query);
  }

  dispose(): void {
    this._disposed = true;
    this._cancel();
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
      const items = await this._fetch(this._query, abort.signal);
      if (abort.signal.aborted || this._disposed)
        return;
      this._state.value = items.length > 0 ? {kind: 'ready', items} : {kind: 'empty'};
    } catch (e) {
      if (abort.signal.aborted || this._disposed)
        return;
      this._state.value = {kind: 'error', message: e instanceof Error ? e.message : String(e)};
    }
  }
}
