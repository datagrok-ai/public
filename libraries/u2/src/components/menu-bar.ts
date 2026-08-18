/* Classic top menu bar: a horizontal row of text items, each opening a fresh Menu anchored to it.
   WAI-ARIA menubar — one tab stop with roving ←/→ focus; while a menu is open, hovering another
   item switches to it instantly (menubar convention, unlike the intent delay of submenus). */
import {signal, ReadonlySignal} from '../core/signals.js';
import {Scope} from '../core/scope.js';
import {Component} from '../core/component.js';
import {Menu} from './menu.js';

export interface MenuBarItem {
  label: string;
  build: (menu: Menu) => void;
  enabled?: boolean;
}

interface Entry {
  item: MenuBarItem;
  el: HTMLButtonElement;
}

export class MenuBar extends Component {
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
    this._listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(this.root, 'focusin', (e) => this._onFocusIn(e));
    this.own(() => this._close());
    for (const item of options?.items ?? [])
      this._add(item);
    this._roving();
  }

  /** `build` runs on every open, so a fresh Menu always reflects current checked state. */
  menu(label: string, build: (menu: Menu) => void): this {
    this._add({label, build});
    this._roving();
    return this;
  }

  setEnabled(label: string, enabled: boolean): void {
    const index = this._entries.findIndex((entry) => entry.item.label === label);
    if (index < 0)
      return;
    const entry = this._entries[index];
    entry.item.enabled = enabled;
    entry.el.setAttribute('aria-disabled', String(!enabled));
    if (!enabled && this._openIndex === index)
      this._close();
    this._roving();
  }

  private _add(item: MenuBarItem): void {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-menu-bar-item';
    el.tabIndex = -1;
    el.textContent = item.label;
    el.setAttribute('role', 'menuitem');
    el.setAttribute('aria-haspopup', 'menu');
    el.setAttribute('aria-expanded', 'false');
    el.setAttribute('aria-disabled', String(item.enabled === false));
    const entry: Entry = {item, el};
    this._listen(el, 'click', () => this._toggle(this._entries.indexOf(entry)));
    this._listen(el, 'pointerenter', () => this._hover(this._entries.indexOf(entry)));
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
    if (!entry || entry.item.enabled === false)
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
    this._open.value = entry.item.label;
    menu.show({anchor: entry.el});

    const watch = new Scope();
    const onKeyDown = (e: Event) => this._onOpenKeyDown(e as KeyboardEvent);
    document.addEventListener('keydown', onKeyDown);
    watch.own(() => document.removeEventListener('keydown', onKeyDown));
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
    return !entry || entry.item.enabled === false;
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

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
