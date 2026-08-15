/* Async-state visuals: the standardized idle/loading/ready/empty/error-with-retry rendering
   over an AsyncSource, so no control hand-rolls a spinner. Combobox renders these states into
   its popup; AsyncView renders them into a content area. */
import {AsyncSource, AsyncFetch, AsyncState} from '../core/async-source.js';
import {Component} from '../core/component.js';
import {Scope} from '../core/scope.js';
import {button} from '../core/elements.js';

export interface AsyncViewOptions {
  empty?: string;
  skeleton?: boolean;
}

export function loader(text?: string): HTMLElement {
  const el = document.createElement('div');
  el.className = 'u2-loader';
  el.setAttribute('role', 'status');
  const spinner = document.createElement('span');
  spinner.className = 'u2-loader-spinner';
  el.append(spinner);
  if (text !== undefined) {
    const label = document.createElement('span');
    label.className = 'u2-loader-text';
    label.textContent = text;
    el.append(label);
  }
  return el;
}

export function skeleton(lines = 3): HTMLElement {
  const el = document.createElement('div');
  el.className = 'u2-skeleton';
  el.setAttribute('aria-hidden', 'true');
  for (let i = 0; i < lines; i++) {
    const line = document.createElement('div');
    line.className = 'u2-skeleton-line';
    el.append(line);
  }
  return el;
}

export class AsyncView<T> extends Component {
  private readonly _source: AsyncSource<T>;
  private readonly _render: (items: T[]) => HTMLElement;
  private readonly _empty: string;
  private readonly _skeleton: boolean;
  private readonly _retry: HTMLElement;
  private _query = '';
  private _content: Scope | undefined;

  constructor(source: AsyncSource<T>, render: (items: T[]) => HTMLElement, options?: AsyncViewOptions) {
    super();
    this._source = source;
    this._render = render;
    this._empty = options?.empty ?? 'No data';
    this._skeleton = options?.skeleton ?? false;
    this.root.classList.add('u2-async-view');
    this.root.dataset.u2 = 'async-view';
    this._retry = this.run(() => button('Retry', () => source.retry()));
    this.own(() => this._releaseContent());
    this.effect(() => this._apply(source.state.value));
  }

  /** Same view over a source it creates and disposes — for the common case where nobody
   * else shares the data. */
  static owned<T>(fetch: AsyncFetch<T>, render: (items: T[]) => HTMLElement,
    options?: AsyncViewOptions): AsyncView<T> {
    const source = new AsyncSource<T>(fetch);
    const view = new AsyncView<T>(source, render, options);
    view.own(() => source.dispose());
    return view;
  }

  refresh(query?: string): void {
    this._query = query ?? this._query;
    this._source.query(this._query);
  }

  private _apply(state: AsyncState<T>): void {
    this._releaseContent();
    const scope = new Scope();
    this._content = scope;
    const content = Scope.runWith(scope, () => this._build(state));
    this.root.setAttribute('aria-busy', String(state.kind === 'loading'));
    this.root.replaceChildren(...(content ? [content] : []));
  }

  private _build(state: AsyncState<T>): HTMLElement | undefined {
    if (state.kind === 'loading')
      return this._skeleton ? skeleton() : loader();
    if (state.kind === 'ready')
      return this._render(state.items);
    if (state.kind === 'empty')
      return AsyncView._div('u2-async-empty', this._empty);
    if (state.kind === 'error') {
      const row = AsyncView._div('u2-async-error');
      row.append(AsyncView._div('u2-async-error-message', state.message), this._retry);
      return row;
    }
    return undefined;
  }

  private _releaseContent(): void {
    if (this._content)
      this._content.dispose();
    this._content = undefined;
  }

  private static _div(cls: string, text?: string): HTMLElement {
    const el = document.createElement('div');
    el.className = cls;
    if (text !== undefined)
      el.textContent = text;
    return el;
  }
}
