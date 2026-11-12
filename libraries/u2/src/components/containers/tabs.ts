/* Tab strip in the platform's `ui.tabControl` look (accent underline, no stripe chrome), with the
   IDE document-tab skin as a variant. Hand-written (no VS Code port — its tab bar is welded to the
   editor-group model and the widget layer we skipped). The overflow dropdown is a plain
   absolutely-positioned div rather than the Overlay layer: it is anchored inside the header and
   must not outlive it. */
import {signal, untracked, Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {iconOf} from '../display/icon.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface TabOptions {
  id: string;
  /** An empty title with an icon makes an icon-only tab (the platform's `d4-tab-icons` mode). */
  title: string;
  /** A Font Awesome name (rendered through `icon`) or a ready element. */
  icon?: string | HTMLElement;
  /** Header tooltip; defaults to the title, so an icon-only tab should set it. */
  tooltip?: string;
  closable?: boolean;
  /** An element or a control, or a function building one; either is placed in the tab's panel on its
   * first activation, unless the strip is {@link TabStripOptions.eager}. */
  content: HTMLElement | Control | (() => HTMLElement | Control);
  /** Why the tab's content needs the user, followed live: the tab is marked while it says something
   * ({@link TabStrip.setInvalid}). Unset and the content is a `Form`, the form's shown error with its field
   * ({@link Form.shownProblem}) — nothing on a fresh form, the first failing field after a Save attempt. */
  validity?: ReadonlySignal<string | null>;
}

export type TabOrientation = 'horizontal' | 'vertical';
export type TabVariant = 'platform' | 'document';

export interface TabStripOptions {
  tabs?: TabOptions[];
  orientation?: TabOrientation;
  /** `platform` (default) is `ui.tabControl`; `document` is the IDE document-tab skin. */
  variant?: TabVariant;
  /** Also the `tabClosed` event (`strip.on('tabClosed', …)`): one call per close, a repeated id included. */
  onTabClosed?: (id: string, strip: TabStrip) => void;
  /** Places a tab's element or control content in its hidden panel when the tab is added, not on its
   * first activation: a form on a tab not yet visited then measures its captions — the tabs of an editor
   * sharing a `LabelColumn` line up before the user switches. A function content stays lazy. */
  eager?: boolean;
}

interface Tab {
  options: TabOptions;
  header: HTMLElement;
  panel: HTMLElement;
  mark: HTMLElement;
  built: boolean;
  /** Why the tab's content needs attention, while it does ({@link TabStrip.setInvalid}). */
  invalid?: string | null;
  /** Owns the tab's validity watch; disposed when the tab closes. */
  scope: Scope;
}

let stripCount = 0;

/** Events: `tabClosed` (the tab id). */
export class TabStrip extends Control {
  readonly activeTab: Signal<string | null> = signal<string | null>(null);
  readonly orientation: TabOrientation;

  private readonly _idPrefix = `u2-tabs-${++stripCount}`;
  private readonly _header = document.createElement('div');
  private readonly _scroller = document.createElement('div');
  private readonly _overflow = document.createElement('button');
  private readonly _menu = document.createElement('div');
  private readonly _content = document.createElement('div');
  private readonly _tabs = signal<Tab[]>([]);
  private readonly _menuOpen = signal(false);
  private readonly _eager: boolean;

  constructor(options?: TabStripOptions) {
    super();
    this.orientation = options?.orientation ?? 'horizontal';
    this._eager = options?.eager === true;

    this.root.classList.add('u2-tabs');
    this.root.classList.toggle('u2-tabs-vertical', this.orientation === 'vertical');
    this.root.classList.toggle('u2-tabs-document', options?.variant === 'document');
    this.root.dataset.u2 = 'tabs';

    this._header.className = 'u2-tabs-header';
    this._header.setAttribute('role', 'tablist');
    this._header.setAttribute('aria-orientation', this.orientation);
    this._scroller.className = 'u2-tabs-scroller';
    this._scroller.setAttribute('role', 'presentation');
    this._overflow.className = 'u2-tabs-overflow';
    this._overflow.type = 'button';
    this._overflow.tabIndex = -1;
    this._overflow.textContent = '⌄';
    this._overflow.title = 'Show all tabs';
    this._overflow.setAttribute('aria-haspopup', 'menu');
    this._menu.className = 'u2-tabs-overflow-menu';
    this._menu.setAttribute('role', 'menu');
    this._header.append(this._scroller, this._overflow, this._menu);
    this._content.className = 'u2-tabs-content';
    this.root.append(this._header, this._content);

    const scope = this.scope;
    scope.listen(this._scroller, 'pointerdown', (e) => this._onPointerDown(e as PointerEvent));
    scope.listen(this._scroller, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    scope.listen(this._overflow, 'click', () => this._menuOpen.value = !this._menuOpen.peek());
    scope.listen(this._menu, 'click', (e) => this._onMenuClick(e));
    const onTabClosed = options?.onTabClosed;
    if (onTabClosed)
      this.on('tabClosed', (id, strip) => onTabClosed(id as string, strip));

    const resize = new ResizeObserver(() => this._updateOverflow());
    resize.observe(this._header);
    this.own(() => resize.disconnect());

    this.effect(() => this._applyActive());
    this.effect(() => this._renderMenu());

    for (const tab of options?.tabs ?? [])
      this.addTab(tab);
  }

  addTab(tab: TabOptions): void {
    const header = document.createElement('div');
    header.className = 'u2-tabs-tab';
    header.id = `${this._idPrefix}-tab-${tab.id}`;
    header.dataset.id = tab.id;
    header.dataset.u2Item = tab.id;
    header.tabIndex = -1;
    header.title = tab.tooltip ?? tab.title;
    header.setAttribute('role', 'tab');
    header.setAttribute('aria-selected', 'false');
    header.setAttribute('aria-controls', `${this._idPrefix}-panel-${tab.id}`);

    const label = document.createElement('span');
    label.className = 'u2-tabs-label';
    label.dataset.u2Part = 'item-label';
    label.textContent = tab.title;
    const mark = document.createElement('span');
    mark.className = tab.closable ? 'u2-tabs-close' : 'u2-tabs-mark';
    if (tab.closable) {
      mark.title = 'Close';
      mark.dataset.u2Part = 'close';
    }
    if (tab.icon !== undefined)
      header.append(TabStrip._icon(tab.icon));
    header.append(label, mark);

    const panel = document.createElement('div');
    panel.className = 'u2-tabs-panel';
    panel.id = `${this._idPrefix}-panel-${tab.id}`;
    panel.style.display = 'none';
    panel.setAttribute('role', 'tabpanel');
    panel.setAttribute('aria-labelledby', header.id);

    this._scroller.append(header);
    this._content.append(panel);
    const eager = this._eager && typeof tab.content !== 'function';
    const entry: Tab = {options: tab, header, panel, mark, built: eager, scope: this.scope.child('tab-validity')};
    if (eager)
      panel.append(tab.content instanceof Control ? tab.content.root : tab.content as HTMLElement);
    this._tabs.value = [...this._tabs.peek(), entry];
    if (tab.validity !== undefined)
      this._watchValidity(entry, tab.validity);
    else if (typeof tab.content !== 'function')
      this._watchContent(entry, tab.content);
    if (this.activeTab.peek() === null)
      this.activeTab.value = tab.id;
    this._updateOverflow();
  }

  private _watchValidity(tab: Tab, validity: ReadonlySignal<string | null>): void {
    tab.scope.effect(() => {
      const message = validity.value;
      untracked(() => this.setInvalid(tab.options.id, message));
    });
  }

  private _watchContent(tab: Tab, content: HTMLElement | Control): void {
    // a Form says which field it shows in error ("Name: Value can't be empty")
    const shown = content as {shownProblem?: unknown, error?: unknown};
    const error = shown.shownProblem instanceof Signal ? shown.shownProblem : shown.error;
    if (content instanceof Control && error instanceof Signal)
      this._watchValidity(tab, error as ReadonlySignal<string | null>);
  }

  closeTab(id: string): void {
    const tabs = this._tabs.peek();
    const index = tabs.findIndex((t) => t.options.id === id);
    if (index < 0)
      return;
    tabs[index].header.remove();
    tabs[index].panel.remove();
    tabs[index].scope.dispose();
    const rest = tabs.slice();
    rest.splice(index, 1);
    this._tabs.value = rest;
    if (this.activeTab.peek() === id)
      this.activeTab.value = rest.length > 0 ? rest[Math.min(index, rest.length - 1)].options.id : null;
    this._updateOverflow();
    this.fireEvent('tabClosed', id);
  }

  /** The base state plus the active tab and the tabs as items. */
  getState(): ReturnType<Control['getState']> {
    const active = this.activeTab.peek();
    return {...super.getState(), ...active === null ? {} : {current: active}};
  }

  getItems(): IItemsState {
    const active = this.activeTab.peek();
    const list = this._tabs.peek().map((tab, index) => {
      const {id, title, tooltip} = tab.options;
      return {key: id, label: title || tooltip || id, index, selected: id === active, current: id === active,
        ...tab.invalid ? {invalid: true} : {}, rendered: true};
    });
    return {count: list.length, list, truncated: false};
  }

  setDirty(id: string, dirty: boolean): void {
    const tab = this._find(id);
    if (tab)
      tab.mark.classList.toggle('u2-tabs-dirty', dirty);
  }

  /** Marks a tab whose content has an error the user should go and fix (`message`, shown as the
   * header's tooltip), or clears it with null: red caption, `aria-invalid`, `invalid` in its item. */
  setInvalid(id: string, message: string | null): void {
    const tab = this._find(id);
    if (tab === undefined)
      return;
    tab.invalid = message;
    tab.header.classList.toggle('u2-tabs-invalid', message !== null);
    if (message === null)
      tab.header.removeAttribute('aria-invalid');
    else
      tab.header.setAttribute('aria-invalid', 'true');
    tab.header.title = message ?? tab.options.tooltip ?? tab.options.title;
  }

  private static _icon(source: string | HTMLElement): HTMLElement {
    const wrap = document.createElement('span');
    wrap.className = 'u2-tabs-icon';
    wrap.append(iconOf(source));
    return wrap;
  }

  /** Renames a tab in place — a bound tab title following its source. */
  setTabTitle(id: string, title: string): void {
    const tab = this._find(id);
    if (tab === undefined)
      return;
    tab.options = {...tab.options, title};
    tab.header.querySelector('.u2-tabs-label')!.textContent = title;
    if (tab.options.tooltip === undefined && !tab.invalid)
      tab.header.title = title;
    this._tabs.value = [...this._tabs.peek()];
  }

  private _find(id: string): Tab | undefined {
    return this._tabs.peek().find((t) => t.options.id === id);
  }

  private _applyActive(): void {
    const active = this.activeTab.value;
    for (const tab of this._tabs.value) {
      const on = tab.options.id === active;
      tab.header.classList.toggle('u2-tabs-tab-active', on);
      tab.header.setAttribute('aria-selected', String(on));
      tab.header.tabIndex = on ? 0 : -1;
      tab.panel.style.display = on ? '' : 'none';
      if (!on)
        continue;
      if (!tab.built) {
        tab.built = true;
        const content = tab.options.content;
        // built once and kept: owned by the strip, not by this effect run (disposed on the next switch)
        this.scope.child('tab').build(() => untracked(() => {
          const built = typeof content === 'function' ? content() : content;
          tab.panel.append(built instanceof Control ? built.root : built);
          if (typeof content === 'function' && tab.options.validity === undefined)
            this._watchContent(tab, built);
        }));
      }
      this._reveal(tab);
    }
  }

  /** Scrolls the strip to the tab without `scrollIntoView`, which also scrolls the page around it. */
  private _reveal(tab: Tab): void {
    const strip = this._scroller.getBoundingClientRect();
    const at = tab.header.getBoundingClientRect();
    this._scroller.scrollLeft += Math.min(0, at.left - strip.left) + Math.max(0, at.right - strip.right);
    this._scroller.scrollTop += Math.min(0, at.top - strip.top) + Math.max(0, at.bottom - strip.bottom);
  }

  private _focus(tab: Tab): void {
    for (const t of this._tabs.peek())
      t.header.tabIndex = t === tab ? 0 : -1;
    tab.header.focus();
    tab.header.scrollIntoView({block: 'nearest', inline: 'nearest'});
  }

  private _onPointerDown(e: PointerEvent): void {
    const target = e.target as HTMLElement;
    const header = target.closest('.u2-tabs-tab') as HTMLElement | null;
    if (!header)
      return;
    const tab = this._find(header.dataset.id!);
    if (!tab)
      return;
    if (e.button === 1) {
      e.preventDefault();
      if (tab.options.closable)
        this.closeTab(tab.options.id);
      return;
    }
    if (e.button !== 0)
      return;
    if (target.closest('.u2-tabs-close')) {
      e.preventDefault();
      this.closeTab(tab.options.id);
      return;
    }
    this.activeTab.value = tab.options.id;
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const tabs = this._tabs.peek();
    const current = tabs.findIndex((t) => t.header === document.activeElement);
    if (current < 0)
      return;
    const last = tabs.length - 1;
    const back = current === 0 ? last : current - 1;
    const forward = current === last ? 0 : current + 1;
    const vertical = this.orientation === 'vertical';
    let next: number;
    switch (e.key) {
      case 'ArrowLeft': next = back; break;
      case 'ArrowRight': next = forward; break;
      case 'ArrowUp':
        if (!vertical)
          return;
        next = back;
        break;
      case 'ArrowDown':
        if (!vertical)
          return;
        next = forward;
        break;
      case 'Home': next = 0; break;
      case 'End': next = last; break;
      case 'Enter':
      case ' ':
        e.preventDefault();
        this.activeTab.value = tabs[current].options.id;
        return;
      case 'Delete': {
        if (!tabs[current].options.closable)
          return;
        e.preventDefault();
        this.closeTab(tabs[current].options.id);
        const rest = this._tabs.peek();
        if (rest.length > 0)
          this._focus(rest[Math.min(current, rest.length - 1)]);
        return;
      }
      default: return;
    }
    e.preventDefault();
    this._focus(tabs[next]);
  }

  private _onMenuClick(e: Event): void {
    const item = (e.target as HTMLElement).closest('.u2-tabs-overflow-item') as HTMLElement | null;
    if (!item)
      return;
    const tab = this._find(item.dataset.id!);
    this._menuOpen.value = false;
    if (!tab)
      return;
    this.activeTab.value = tab.options.id;
    tab.header.scrollIntoView({block: 'nearest', inline: 'nearest'});
  }

  private _onDocumentPointerDown(e: Event): void {
    const target = e.target as Node;
    if (!this._menu.contains(target) && !this._overflow.contains(target))
      this._menuOpen.value = false;
  }

  private _onDocumentKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Escape') {
      this._menuOpen.value = false;
      this._overflow.focus();
    }
  }

  private _renderMenu(): void {
    const open = this._menuOpen.value;
    this._menu.classList.toggle('u2-tabs-overflow-menu-open', open);
    if (!open) {
      this._menu.textContent = '';
      return;
    }
    // the page is watched only while the menu is open: this run's scope drops the listeners
    const run = Scope.ambient!;
    run.listen(document, 'pointerdown', (e) => this._onDocumentPointerDown(e));
    run.listen(document, 'keydown', (e) => this._onDocumentKeyDown(e as KeyboardEvent));
    const active = this.activeTab.value;
    this._menu.replaceChildren(...this._tabs.value.map((tab) => {
      const {id, title, icon, tooltip} = tab.options;
      const item = document.createElement('div');
      item.className = 'u2-tabs-overflow-item';
      item.dataset.id = id;
      item.setAttribute('role', 'menuitem');
      item.classList.toggle('u2-tabs-overflow-item-active', id === active);
      // the header owns its icon element, so the menu gets a copy
      if (icon !== undefined)
        item.append(TabStrip._icon(typeof icon === 'string' ? icon : icon.cloneNode(true) as HTMLElement));
      item.append(title || tooltip || id);
      return item;
    }));
  }

  private _updateOverflow(): void {
    const overflowing = this.orientation === 'vertical' ?
      this._scroller.scrollHeight > this._scroller.clientHeight + 1 :
      this._scroller.scrollWidth > this._scroller.clientWidth + 1;
    this._overflow.classList.toggle('u2-tabs-overflow-visible', overflowing);
    this._overflow.tabIndex = overflowing ? 0 : -1;
    if (!overflowing)
      this._menuOpen.value = false;
  }
}
