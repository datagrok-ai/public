import {signal, Signal, ReadonlySignal, batch} from './signals.js';

export type PagerState = 'idle' | 'loading' | 'done' | 'error';

export interface AsyncPagerOptions<T> {
  /** `page` is 0-based; a page shorter than `pageSize` is the last one. */
  fetchPage: (page: number, abort: AbortSignal) => Promise<T[]>;
  count?: (abort: AbortSignal) => Promise<number>;
  pageSize?: number;
}

/** The paged twin of `AsyncSource`: pages accumulate into one array signal, a short page ends
 * the collection, and responses a `reset()` outran are dropped — written once so no browse-style
 * list hand-rolls page counters, end-of-data checks or racing. `reset()` also starts a collection:
 * it is what reads the total and loads the first page. */
export class AsyncPager<T> {
  readonly items: ReadonlySignal<T[]>;
  readonly state: ReadonlySignal<PagerState>;
  readonly total: ReadonlySignal<number | null>;

  private readonly _items: Signal<T[]>;
  private readonly _state: Signal<PagerState>;
  private readonly _total: Signal<number | null>;
  private readonly _fetchPage: (page: number, abort: AbortSignal) => Promise<T[]>;
  private readonly _count: ((abort: AbortSignal) => Promise<number>) | undefined;
  private readonly _pageSize: number;
  private _page = 0;
  private _abort: AbortController | undefined;
  private _disposed = false;

  constructor(options: AsyncPagerOptions<T>) {
    this._fetchPage = options.fetchPage;
    this._count = options.count;
    this._pageSize = options.pageSize ?? 20;
    this._items = signal<T[]>([]);
    this._state = signal<PagerState>('idle');
    this._total = signal<number | null>(null);
    this.items = this._items;
    this.state = this._state;
    this.total = this._total;
  }

  /** Loads the next page; a no-op while one is in flight and once the last one arrived. After an
   * error it retries the page that failed, so the accumulated pages survive a hiccup. */
  loadMore(): void {
    const state = this._state.peek();
    if (this._disposed || state === 'loading' || state === 'done')
      return;
    this._state.value = 'loading';
    void this._load(this._abortSignal());
  }

  reset(): void {
    this._cancel();
    if (this._disposed)
      return;
    this._page = 0;
    batch(() => {
      this._items.value = [];
      this._total.value = null;
      this._state.value = 'idle';
    });
    if (this._count)
      void this._recount(this._abortSignal());
    this.loadMore();
  }

  dispose(): void {
    this._disposed = true;
    this._cancel();
  }

  /** One controller per generation: everything in flight is invalidated together by `reset()`. */
  private _abortSignal(): AbortSignal {
    if (!this._abort)
      this._abort = new AbortController();
    return this._abort.signal;
  }

  private _cancel(): void {
    if (this._abort)
      this._abort.abort();
    this._abort = undefined;
  }

  private async _load(abort: AbortSignal): Promise<void> {
    try {
      const page = await this._fetchPage(this._page, abort);
      if (abort.aborted || this._disposed)
        return;
      this._page++;
      batch(() => {
        this._items.value = this._items.peek().concat(page);
        const done = page.length < this._pageSize;
        this._state.value = done ? 'done' : 'idle';
        if (done && this._total.peek() === null)
          this._total.value = this._items.peek().length;
      });
    } catch (_) {
      if (abort.aborted || this._disposed)
        return;
      this._state.value = 'error';
    }
  }

  private async _recount(abort: AbortSignal): Promise<void> {
    try {
      const total = await this._count!(abort);
      if (!abort.aborted && !this._disposed)
        this._total.value = total;
    } catch (_) { /* the total is a nicety — a failed count leaves it unknown, the pages stand */ }
  }
}
