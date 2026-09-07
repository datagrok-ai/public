/* Server collections as async-control sources. Structural on purpose: nothing here imports
   datagrok-api, so a real `grok.dapi.*` source and a headless fake are equally acceptable. */
import {AsyncFetch} from '../../core/async-source.js';
import {AsyncPager} from '../../core/pager.js';

/** What the adapter needs from `HttpDataSource<T>`; every `grok.dapi.*` source satisfies it. */
export interface DapiSourceLike<T> {
  list(options?: {pageSize?: number, filter?: string}): Promise<T[]>;
}

/** What paging needs on top: a page number, a total, and the fluent setters `count()` reads
 * (it takes no options, so a filtered total has to be set on the source first). */
export interface DapiPagerSourceLike<T> extends DapiSourceLike<T> {
  list(options?: {pageSize?: number, pageNumber?: number, filter?: string}): Promise<T[]>;
  count(): Promise<number>;
  filter(query: string): DapiPagerSourceLike<T>;
  order(field: string, desc?: boolean): DapiPagerSourceLike<T>;
}

export interface DapiSourceOptions {
  pageSize?: number;
  /** Query → smart-filter string (`undefined` lists unfiltered). Default: the trimmed query
   * itself, free-text smart search. */
  filter?: (query: string) => string | undefined;
}

/** An `AsyncFetch` over a dapi collection, for `TypeAhead`/`Combobox`/`AsyncSource`:
 * `dapiSource(() => grok.dapi.users)`. Takes a factory, not an instance — `HttpDataSource`'s
 * fluent methods mutate the instance, so a shared one would accumulate filters across queries;
 * every `grok.dapi.X` getter builds a fresh source. The server offers no cancellation, so the
 * abort signal's only job is done upstream: `AsyncSource` drops stale responses. */
export function dapiSource<T>(factory: () => DapiSourceLike<T>, options: DapiSourceOptions = {}): AsyncFetch<T> {
  const toFilter = options.filter ?? ((q: string) => q === '' ? undefined : q);
  return (query: string) => {
    const filter = toFilter(sanitizeFilterValue(query));
    return factory().list({pageSize: options.pageSize ?? 20, ...(filter === undefined ? {} : {filter})});
  };
}

export interface DapiPagerOptions {
  pageSize?: number;
  /** A thunk is read on every request, so one pager serves a query box for its whole lifetime:
   * change the query, call `reset()`, and the next generation carries the new filter. */
  filter?: string | (() => string | undefined);
  order?: string;
  desc?: boolean;
}

/** An `AsyncPager` over a dapi collection: `dapiPager(() => grok.dapi.reports, {pageSize: 30})`.
 * Takes a factory for the same reason {@link dapiSource} does — the fluent setters mutate the
 * source — so the total and every page are ordered and filtered on a source of their own. */
export function dapiPager<T>(factory: () => DapiPagerSourceLike<T>,
  options: DapiPagerOptions = {}): AsyncPager<T> {
  const pageSize = options.pageSize ?? 20;
  const currentFilter = () => typeof options.filter === 'function' ? options.filter() : options.filter;
  const source = () => options.order === undefined ? factory() : factory().order(options.order, options.desc);
  return new AsyncPager<T>({
    pageSize,
    fetchPage: (page) => {
      const filter = currentFilter();
      return source().list({pageSize, pageNumber: page, ...(filter === undefined ? {} : {filter})});
    },
    count: () => {
      const filter = currentFilter();
      return (filter === undefined ? source() : source().filter(filter)).count();
    },
  });
}

/** Quotes and backslashes break the smart-filter grammar; strip them from user-typed values
 * before they are spliced into a `like "…"` clause. */
export function sanitizeFilterValue(query: string): string {
  return query.replace(/["\\]/g, '').trim();
}
