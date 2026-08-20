/* IDE document tabs: hand-written (no VS Code port — its tab bar is welded to the editor-group
   model and the widget layer we skipped). The overflow dropdown is a plain absolutely-positioned
   div rather than the Overlay layer: it is anchored inside the header and must not outlive it. */
import {signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';

export interface TabOptions {
  id: string;
  label: string;
  closable?: boolean;
  content: HTMLElement | (() => HTMLElement);
}

interface Tab {
  options: TabOptions;
  header: HTMLElement;
  panel: HTMLElement;
  mark: HTMLElement;
  built: boolean;
}

let stripCount = 0;

export class TabStrip extends Control {
  readonly activeTab: Signal<string | null> = signal<string | null>(null);
  readonly onTabClosed: ReadonlySignal<string | null>;

  private readonly _idPrefix = `u2-tabs-${++stripCount}`;
  private readonly _header = document.createElement('div');
  private readonly _scroller = document.createElement('div');
  private readonly _overflow = document.createElement('button');
  private readonly _menu = document.createElement('div');
  private readonly _content = document.createElement('div');
  private readonly _tabs = signal<Tab[]>([]);
  private readonly _menuOpen = signal(false);
  private readonly _closed = signal<string | null>(null);

  constructor(options?: {tabs?: TabOptions[]}) {
    super();
    this.onTabClosed = this._closed;

    this.root.classList.add('u2-tabs');
    this.root.dataset.u2 = 'tabs';

    this._header.className = 'u2-tabs-header';
    this._header.setAttribute('role', 'tablist');
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

    this._listen(this._scroller, 'pointerdown', (e) => this._onPointerDown(e as PointerEvent));
    this._listen(this._scroller, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(this._overflow, 'click', () => this._menuOpen.value = !this._menuOpen.peek());
    this._listen(this._menu, 'click', (e) => this._onMenuClick(e));
    this._listen(document, 'pointerdown', (e) => this._onDocumentPointerDown(e));
    this._listen(document, 'keydown', (e) => this._onDocumentKeyDown(e as KeyboardEvent));

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
    header.tabIndex = -1;
    header.title = tab.label;
    header.setAttribute('role', 'tab');
    header.setAttribute('aria-selected', 'false');
    header.setAttribute('aria-controls', `${this._idPrefix}-panel-${tab.id}`);

    const label = document.createElement('span');
    label.className = 'u2-tabs-label';
    label.textContent = tab.label;
    const mark = document.createElement('span');
    mark.className = tab.closable ? 'u2-tabs-close' : 'u2-tabs-mark';
    if (tab.closable)
      mark.title = 'Close';
    header.append(label, mark);

    const panel = document.createElement('div');
    panel.className = 'u2-tabs-panel';
    panel.id = `${this._idPrefix}-panel-${tab.id}`;
    panel.style.display = 'none';
    panel.setAttribute('role', 'tabpanel');
    panel.setAttribute('aria-labelledby', header.id);

    this._scroller.append(header);
    this._content.append(panel);
    this._tabs.value = [...this._tabs.peek(), {options: tab, header, panel, mark, built: false}];
    if (this.activeTab.peek() === null)
      this.activeTab.value = tab.id;
    this._updateOverflow();
  }

  closeTab(id: string): void {
    const tabs = this._tabs.peek();
    const index = tabs.findIndex((t) => t.options.id === id);
    if (index < 0)
      return;
    tabs[index].header.remove();
    tabs[index].panel.remove();
    const rest = tabs.slice();
    rest.splice(index, 1);
    this._tabs.value = rest;
    if (this.activeTab.peek() === id)
      this.activeTab.value = rest.length > 0 ? rest[Math.min(index, rest.length - 1)].options.id : null;
    this._updateOverflow();
    this._closed.value = id;
  }

  setDirty(id: string, dirty: boolean): void {
    const tab = this._find(id);
    if (tab)
      tab.mark.classList.toggle('u2-tabs-dirty', dirty);
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
        untracked(() => tab.panel.append(typeof content === 'function' ? content() : content));
      }
      tab.header.scrollIntoView({block: 'nearest', inline: 'nearest'});
    }
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
    let next: number;
    switch (e.key) {
      case 'ArrowLeft': next = current === 0 ? last : current - 1; break;
      case 'ArrowRight': next = current === last ? 0 : current + 1; break;
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
    if (!this._menuOpen.peek())
      return;
    const target = e.target as Node;
    if (!this._menu.contains(target) && !this._overflow.contains(target))
      this._menuOpen.value = false;
  }

  private _onDocumentKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Escape' && this._menuOpen.peek()) {
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
    const active = this.activeTab.value;
    this._menu.replaceChildren(...this._tabs.value.map((tab) => {
      const item = document.createElement('div');
      item.className = 'u2-tabs-overflow-item';
      item.dataset.id = tab.options.id;
      item.textContent = tab.options.label;
      item.setAttribute('role', 'menuitem');
      item.classList.toggle('u2-tabs-overflow-item-active', tab.options.id === active);
      return item;
    }));
  }

  private _updateOverflow(): void {
    const overflowing = this._scroller.scrollWidth > this._scroller.clientWidth + 1;
    this._overflow.classList.toggle('u2-tabs-overflow-visible', overflowing);
    this._overflow.tabIndex = overflowing ? 0 : -1;
    if (!overflowing)
      this._menuOpen.value = false;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
