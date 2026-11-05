/* Classic top menu bar: a horizontal row of text items, each opening a fresh Menu anchored to it.
   WAI-ARIA menubar — one tab stop with roving ←/→ focus; while a menu is open, hovering another
   item switches to it instantly (menubar convention, unlike the intent delay of submenus). */
import {signal, ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {Reason} from '../../core/reason.js';
import {Tooltip} from '../../core/tooltip.js';
import {Menu} from './menu.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface MenuBarItem {
  text: string;
  build: (menu: Menu) => void;
  /** Why the item opens nothing now — greyed, the reason as its tooltip; `null` or absent while it opens. */
  disabled?: string | null;
}

interface Entry {
  item: MenuBarItem;
  el: HTMLButtonElement;
  off: string | null;
}

export class MenuBar extends Control {
  /** Label of the open item, null while the bar is closed. */
  readonly openMenu: ReadonlySignal<string | null>;

  private readonly _entries: Entry[] = [];
  private readonly _open = signal<string | null>(null);
  private _menu: Menu | undefined;
  private _watch: Scope | undefined;
  private _openIndex = -1;
  private _active = 0;

  constructor(options?: {items?: MenuBarItem[]}) {
    super();
    this.openMenu = this._open;
    this.root.classList.add('u2-menu-bar');
    this.root.dataset.u2 = 'menu-bar';
    this.root.setAttribute('role', 'menubar');
    this.scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this.scope.listen(this.root, 'focusin', (e) => this._onFocusIn(e));
    this.own(() => this._close());
    for (const item of options?.items ?? [])
      this._add(item);
    this._roving();
  }

  /** `build` runs on every open, so a fresh Menu always reflects current checked state. */
  menu(text: string, build: (menu: Menu) => void): this {
    this._add({text, build});
    this._roving();
    return this;
  }

  /** Turns the item `text` off with `reason`, or on with `null`. */
  setDisabled(text: string, reason: string | null): void {
    const index = this._entries.findIndex((entry) => entry.item.text === text);
    if (index < 0)
      return;
    const entry = this._entries[index];
    MenuBar._apply(entry, reason || null);
    if (entry.off !== null && this._openIndex === index)
      this._close();
    this._roving();
  }

  /** The items, keyed by text; the open one is expanded. */
  getItems(): IItemsState {
    const list = this._entries.map((entry, index) => ({key: entry.item.text, label: entry.item.text, index,
      disabled: entry.off !== null, expanded: index === this._openIndex, rendered: true}));
    return {count: list.length, list, truncated: false};
  }

  private static _apply(entry: Entry, reason: string | null): void {
    entry.off = reason;
    Reason.apply(entry.el, reason);
    if (reason === null)
      entry.el.setAttribute('aria-disabled', 'false');
  }

  private _add(item: MenuBarItem): void {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-menu-bar-item';
    el.tabIndex = -1;
    el.textContent = item.text;
    el.dataset.u2Item = item.text;
    el.setAttribute('role', 'menuitem');
    el.setAttribute('aria-haspopup', 'menu');
    el.setAttribute('aria-expanded', 'false');
    const entry: Entry = {item, el, off: null};
    MenuBar._apply(entry, item.disabled || null);
    Tooltip.bind(el, () => entry.off ?? '', {owner: this.scope});
    this.scope.listen(el, 'click', () => this._toggle(this._entries.indexOf(entry)));
    this.scope.listen(el, 'pointerenter', () => this._hover(this._entries.indexOf(entry)));
    this._entries.push(entry);
    this.root.append(el);
  }

  private _toggle(index: number): void {
    if (this._openIndex === index)
      this._close();
    else
      this._show(index);
  }

  private _show(index: number): void {
    const entry = this._entries[index];
    if (!entry || entry.off !== null)
      return;
    this._close();
    this._active = index;
    this._roving();
    // focused before the menu opens, so Menu returns focus here when it closes
    entry.el.focus();
    const menu = new Menu();
    entry.item.build(menu);
    this._menu = menu;
    this._openIndex = index;
    this._open.value = entry.item.text;
    menu.show({anchor: entry.el});

    const watch = new Scope(undefined, {label: 'menu-bar-open'});
    watch.listen(document, 'keydown', (e) => this._onOpenKeyDown(e as KeyboardEvent));
    watch.effect(() => {
      const open = menu.isOpen.value;
      entry.el.setAttribute('aria-expanded', String(open));
      entry.el.classList.toggle('u2-menu-bar-item-open', open);
      if (open)
        return;
      if (this._menu === menu) {
        this._menu = undefined;
        this._watch = undefined;
        this._openIndex = -1;
        this._open.value = null;
      }
      watch.dispose();
    });
    this._watch = watch;
  }

  private _close(): void {
    if (this._menu)
      this._menu.close();
    if (this._watch)
      this._watch.dispose();
    this._menu = undefined;
    this._watch = undefined;
    this._openIndex = -1;
    this._open.value = null;
  }

  private _hover(index: number): void {
    if (this._openIndex >= 0 && this._openIndex !== index)
      this._show(index);
  }

  private _move(index: number): void {
    if (index < 0)
      return;
    this._active = index;
    this._roving();
    this._entries[index].el.focus();
    if (this._openIndex >= 0 && this._openIndex !== index)
      this._show(index);
  }

  private _roving(): void {
    if (this._active >= this._entries.length || this._disabled(this._active))
      this._active = this._edge(1);
    for (let i = 0; i < this._entries.length; i++)
      this._entries[i].el.tabIndex = i === this._active ? 0 : -1;
  }

  private _disabled(index: number): boolean {
    const entry = this._entries[index];
    return !entry || entry.off !== null;
  }

  private _step(from: number, delta: number): number {
    const count = this._entries.length;
    for (let k = 1; k <= count; k++) {
      const i = ((from + delta * k) % count + count) % count;
      if (!this._disabled(i))
        return i;
    }
    return -1;
  }

  private _edge(delta: number): number {
    const count = this._entries.length;
    for (let i = delta > 0 ? 0 : count - 1; i >= 0 && i < count; i += delta) {
      if (!this._disabled(i))
        return i;
    }
    return -1;
  }

  private _onFocusIn(e: Event): void {
    const index = this._entries.findIndex((entry) => entry.el === e.target);
    if (index < 0)
      return;
    this._active = index;
    this._roving();
  }

  private _onKeyDown(e: KeyboardEvent): void {
    switch (e.key) {
      case 'ArrowLeft':
      case 'ArrowRight':
        e.preventDefault();
        this._move(this._step(this._active, e.key === 'ArrowRight' ? 1 : -1));
        break;
      case 'Home':
        e.preventDefault();
        this._move(this._edge(1));
        break;
      case 'End':
        e.preventDefault();
        this._move(this._edge(-1));
        break;
      case 'ArrowDown':
      case 'Enter':
      case ' ':
        e.preventDefault();
        if (this._openIndex !== this._active)
          this._show(this._active);
        break;
      case 'Escape':
        this._close();
        break;
    }
  }

  /** Focus sits in the open Menu, so ←/→ reach the bar as document keys; anything the menu itself
   * consumed (a submenu open or close) arrives already defaulted-prevented. */
  private _onOpenKeyDown(e: KeyboardEvent): void {
    if (e.defaultPrevented || this._openIndex < 0)
      return;
    if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight')
      return;
    e.preventDefault();
    this._move(this._step(this._active, e.key === 'ArrowRight' ? 1 : -1));
  }

}
