/* Saved runs of one function, newest first, paged off `grok.dapi.functions.calls` and rendered
   through the FuncCall ObjectHandler (FuncCallMeta.renderListItem). `functionName` drives it:
   a change resets the pager and reloads; selection goes into `grok.shell.o` like the
   FunctionsBrowser's does. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {computed, signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {button, div, divV, span} from '../../core/elements.js';
import type {BindProp, BindSource} from '../../core/widget-like.js';
import type {LiveOption} from '../../core/input-base.js';
import {AsyncPager} from '../../core/pager.js';
import {VirtualList} from '../../components/collections/list.js';
import {handlerRenderer} from './entity.js';
import {sanitizeFilterValue} from './dapi-source.js';

export interface FuncCallHistoryBrowserOptions {
  /** Whose runs to list (qualified or short name). Live (UB-10, follow-only): a change resets
   * the pager and reloads. */
  functionName?: LiveOption<string>;
  pageSize?: number;
  itemHeight?: number;
  /** Selection makes the call the shell's current object; default true. */
  setCurrentObject?: boolean;
  /** Row override; `row` is the pooled list row (per-row attributes land there). */
  render?: (call: DG.FuncCall, row: HTMLElement) => HTMLElement;
  onChanged?: (call: DG.FuncCall | null) => void;
  /** Double-click or Enter. */
  onActivate?: (call: DG.FuncCall) => void;
}

/** Everything the Dart row renderer and the populate path need in one fetch. */
const INCLUDE = 'session.user,func.package,func.params,inputs,options';

export class FuncCallHistoryBrowser extends Control {
  /** The control's own writable signal; a bound option is followed into it, never adopted. */
  readonly functionName: Signal<string>;
  readonly selected: ReadonlySignal<DG.FuncCall | null>;

  private readonly options: FuncCallHistoryBrowserOptions;
  private readonly _pager: AsyncPager<DG.FuncCall>;
  private readonly _list: VirtualList<DG.FuncCall>;
  private readonly _renderer = handlerRenderer<DG.FuncCall>();
  private readonly _selected = signal<DG.FuncCall | null>(null);
  private _lastItems: DG.FuncCall[] = [];
  private _scrollPending = false;
  private _selectedId?: ReadonlySignal<string | null>;

  constructor(options: FuncCallHistoryBrowserOptions = {}) {
    super();
    this.options = options;
    this.root.classList.add('u2-fch');
    this.root.dataset.u2 = 'func-call-history-browser';
    this.functionName = this._live(options.functionName, '');

    const pageSize = options.pageSize ?? 20;
    this._pager = new AsyncPager<DG.FuncCall>({
      pageSize,
      fetchPage: (page) => this._fetchPage(page, pageSize),
    });
    this.own(() => this._pager.dispose());

    const itemHeight = options.itemHeight ?? 36;
    this._list = this.runInScope(() => new VirtualList<DG.FuncCall>({
      itemHeight,
      keyOf: (call) => call.id,
      render: (call, _index, row) => this._renderRow(call, row),
    }));
    this._list.root.setAttribute('data-u2', 'fch-list');
    this._list.setItems(this._pager.items);
    const onScroll = () => {
      const root = this._list.root;
      if (root.scrollTop + root.clientHeight > root.scrollHeight - 5 * itemHeight)
        this._pager.loadMore();
    };
    this._list.root.addEventListener('scroll', onScroll);
    this.own(() => this._list.root.removeEventListener('scroll', onScroll));

    const emptyMessage = span('', 'u2-fch-empty-message');
    const retry = this.runInScope(() => button('Retry', () => this._pager.loadMore()));
    retry.dataset.u2 = 'fch-retry';
    const stateArea = divV([emptyMessage, retry], 'u2-async-empty u2-fch-state');
    stateArea.dataset.u2 = 'fch-state';
    this.root.append(divV([this._list, stateArea], 'u2-fch-main'));
    this.effect(() => {
      const state = this._pager.state.value;
      const empty = this._pager.items.value.length === 0;
      const name = this.functionName.value;
      retry.hidden = state !== 'error';
      emptyMessage.textContent =
        state === 'error' ? 'Could not load the runs.' :
          empty && state === 'loading' ? 'Loading runs…' :
            empty && name === '' ? 'Pick a function to see its runs.' :
              empty && state === 'done' ? 'No saved runs.' : '';
      stateArea.hidden = !(empty || state === 'error');
      this._list.root.hidden = empty;
    });

    this.selected = this._selected;
    // never a naive (items, selectedIndex) pairing: on an items change the keyed reselect lives
    // in the list's own effect with no ordering guarantee against this one, so for one flush the
    // OLD index points into the NEW array (functions-browser.ts precedent, keyed by call id)
    this.effect(() => {
      const items = this._pager.items.value;
      const index = this._list.selectedIndex.value;
      const call = items[index] ?? null;
      const key = call?.id ?? null;
      const prevKey = this._selected.peek()?.id ?? null;
      const itemsChanged = items !== this._lastItems;
      this._lastItems = items;
      if (itemsChanged && key !== prevKey) {
        if (prevKey !== null && items.some((c) => c.id === prevKey)) {
          this._scrollPending = true;
          return;
        }
        this._selected.value = null;
        this._scrollPending = false;
        return;
      }
      this._selected.value = call;
      if ((itemsChanged || this._scrollPending) && index >= 0)
        this._list.scrollToIndex(index);
      this._scrollPending = false;
    });

    this.effect(() => {
      this.functionName.value;
      this._pager.reset();
    });

    let initial = true;
    this.effect(() => {
      const call = this.selected.value;
      if (initial) {
        initial = false;
        return;
      }
      if (call !== null && this.options.setCurrentObject !== false)
        grok.shell.o = call;
      this.options.onChanged?.(call);
      this.fireEvent('change', call);
    });

    this._listen(this._list.root, 'dblclick', (e) => {
      const target = e.target as Element | null;
      const row = target?.closest('.u2-list-row') as HTMLElement | null;
      if (row)
        this._activate(this._pager.items.peek()[Number(row.dataset.index)]);
    });
    this._listen(this._list.root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter')
        this._activate(this.selected.peek() ?? undefined);
    });
  }

  refresh(): void {
    this._pager.reset();
  }

  /** A loaded call by id — the popup commit seam (the row carries the id it renders). */
  call(id: string): DG.FuncCall | undefined {
    return this._pager.items.peek().find((c) => c.id === id);
  }

  /** `$.fch.functionName` follows a FunctionsBrowser's `selected` step; `$.fch.selected`
   * answers the picked call's id. */
  bindStep(name: string): Signal<unknown> | BindSource | null {
    if (name === 'functionName')
      return this.functionName as unknown as Signal<unknown>;
    if (name === 'selected') {
      this._selectedId ??= computed(() => this.selected.value?.id ?? null);
      return this._selectedId as unknown as Signal<unknown>;
    }
    return super.bindStep(name);
  }

  bindProps(): BindProp[] {
    const props = super.bindProps();
    if (!props.some((p) => p.name === 'functionName'))
      props.push({name: 'functionName', type: 'string', writable: true});
    props.push({name: 'selected', type: 'string', writable: false});
    return props;
  }

  /** The browser in an anchored popup — the one shape every history picker shares: a row click
   * or activation picks a LOADED call and closes; Esc, an outside click or the anchor leaving
   * the DOM close without a pick. Disposing the returned scope closes the popup; hang
   * close-state resets on it with `scope.own`. */
  static popup(anchor: HTMLElement, popupClass: string, options: FuncCallHistoryBrowserOptions,
    onPick: (call: DG.FuncCall) => void): Scope {
    const scope = new Scope();
    const pick = (call: DG.FuncCall | undefined): void => {
      if (call === undefined)
        return;
      scope.dispose();
      onPick(call);
    };
    Scope.runWith(scope, () => {
      const browser = new FuncCallHistoryBrowser({setCurrentObject: false, ...options,
        onActivate: pick});
      const popup = div([browser.root], popupClass);
      // one click picks: the row carries the id it renders, so the commit never depends on the
      // browser's selection effects having flushed
      popup.addEventListener('click', (e) => {
        const id = (e.target as Element | null)?.closest('[data-u2-funccall]')
          ?.getAttribute('data-u2-funccall');
        if (id)
          pick(browser.call(id));
      });
      popup.addEventListener(OVERLAY_CLOSE_EVENT, () => scope.dispose());
      Overlay.show(anchor, popup, scope);
    });
    return scope;
  }

  // the proven filter term is func.name="<short>" (compute-utils history-utils precedent);
  // cross-package same-name collisions accepted, as in the Dart history panes
  private _fetchPage(page: number, pageSize: number): Promise<DG.FuncCall[]> {
    const name = this.functionName.peek();
    if (name === '')
      return Promise.resolve([]);
    const short = sanitizeFilterValue(name.split(':').pop() ?? name);
    return grok.dapi.functions.calls.allPackageVersions()
      .filter(`func.name="${short}"`)
      .include(INCLUDE)
      .order('started', true)
      .list({pageSize, pageNumber: page});
  }

  private _live(value: LiveOption<string> | undefined, fallback: string): Signal<string> {
    if (!(value instanceof Signal))
      return signal((value as string | undefined) ?? fallback);
    const own = signal((value.peek() as string | undefined) ?? fallback);
    this.effect(() => {
      const x = (value as ReadonlySignal<string>).value;
      if (x !== undefined)
        own.value = x; // a forward-ref proxy starts undefined; skip until linked
    });
    return own;
  }

  private _activate(call: DG.FuncCall | undefined): void {
    if (call === undefined)
      return;
    this.options.onActivate?.(call);
    this.fireEvent('activate', call);
  }

  private _renderRow(call: DG.FuncCall, row: HTMLElement): HTMLElement {
    row.setAttribute('data-u2-funccall', call.id);
    row.setAttribute('aria-label', this._renderer.caption(call));
    const custom = this.options.render;
    if (custom)
      return custom(call, row);
    return div([this._renderer.listItem(call)], 'u2-fch-row');
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}

export function funcCallHistoryBrowser(options: FuncCallHistoryBrowserOptions = {}): FuncCallHistoryBrowser {
  return new FuncCallHistoryBrowser(options);
}
