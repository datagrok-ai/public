/* IDE toolbar: one compact row that collapses whole items into a » menu when it runs out of
   width. WAI-ARIA toolbar pattern — a single tab stop with roving ←/→ focus. */
import {Signal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {Tooltip} from '../../core/tooltip.js';
import {Menu} from './menu.js';

export interface ToolbarButtonOptions {
  tooltip?: string;
  icon?: string;
}

type ItemKind = 'button' | 'toggle' | 'menu' | 'separator' | 'slot';

interface Item {
  kind: ItemKind;
  el: HTMLElement;
  label: string;
  hidden: boolean;
  icon?: string;
  onClick?: () => void;
  bind?: Signal<boolean>;
  build?: (menu: Menu) => void;
}

const HIDDEN = 'u2-toolbar-hidden';
const FOCUSABLE = 'input, button, select, textarea, a[href], [tabindex]';
const CARET = 'input:not([type=button]):not([type=checkbox]):not([type=radio]), textarea';

export class Toolbar extends Control {
  private readonly _items: Item[] = [];
  private readonly _chevron = document.createElement('button');
  private _menu: Menu | undefined;
  private _menuWatch: Scope | undefined;
  private _menuOwner: HTMLElement | undefined;
  private _active = 0;
  private _width = -1;

  constructor(options?: {ariaLabel?: string}) {
    super();
    this.root.classList.add('u2-toolbar');
    this.root.dataset.u2 = 'toolbar';
    this.root.setAttribute('role', 'toolbar');
    if (options && options.ariaLabel)
      this.root.setAttribute('aria-label', options.ariaLabel);

    this._chevron.type = 'button';
    this._chevron.className = `u2-toolbar-btn u2-toolbar-overflow ${HIDDEN}`;
    this._chevron.textContent = '»';
    this._chevron.tabIndex = -1;
    this._chevron.setAttribute('aria-label', 'More');
    this._chevron.setAttribute('aria-haspopup', 'menu');
    this._chevron.setAttribute('aria-expanded', 'false');
    this.root.append(this._chevron);

    this._listen(this._chevron, 'click', () => this._toggleMenu(this._chevron, (m) => this._buildOverflow(m)));
    this._listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(this.root, 'focusin', (e) => this._onFocusIn(e));

    const resize = new ResizeObserver(() => this._layout());
    resize.observe(this.root);
    this.own(() => resize.disconnect());
    this.own(() => this._closeMenu());
  }

  addButton(label: string, onClick: () => void, options?: ToolbarButtonOptions): Toolbar {
    const el = this._button(label, options);
    this._listen(el, 'click', () => onClick());
    return this._add({kind: 'button', el, label, hidden: false, onClick, icon: options ? options.icon : undefined});
  }

  addToggle(label: string, bind: Signal<boolean>, options?: {tooltip?: string}): Toolbar {
    const el = this._button(label, options);
    el.setAttribute('aria-pressed', 'false');
    this._listen(el, 'click', () => bind.value = !bind.peek());
    this.effect(() => {
      const pressed = bind.value;
      el.classList.toggle('u2-toolbar-btn-pressed', pressed);
      el.setAttribute('aria-pressed', String(pressed));
    });
    return this._add({kind: 'toggle', el, label, hidden: false, bind});
  }

  /** `build` runs on every open, so a fresh Menu always reflects current checked state. */
  addMenu(label: string, build: (menu: Menu) => void): Toolbar {
    const el = this._button(label);
    el.classList.add('u2-toolbar-btn-menu');
    el.setAttribute('aria-haspopup', 'menu');
    el.setAttribute('aria-expanded', 'false');
    const arrow = document.createElement('span');
    arrow.className = 'u2-toolbar-arrow';
    arrow.textContent = '▾';
    el.append(arrow);
    this._listen(el, 'click', () => this._toggleMenu(el, build));
    return this._add({kind: 'menu', el, label, hidden: false, build});
  }

  addSeparator(): Toolbar {
    const el = document.createElement('div');
    el.className = 'u2-toolbar-sep';
    el.setAttribute('role', 'separator');
    el.setAttribute('aria-orientation', 'vertical');
    return this._add({kind: 'separator', el, label: '', hidden: false});
  }

  add(el: HTMLElement): Toolbar {
    const slot = document.createElement('div');
    slot.className = 'u2-toolbar-slot';
    slot.append(el);
    return this._add({kind: 'slot', el: slot, label: '', hidden: false});
  }

  private _button(label: string, options?: ToolbarButtonOptions): HTMLButtonElement {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-toolbar-btn';
    el.tabIndex = -1;
    if (options && options.icon) {
      const icon = document.createElement('span');
      icon.className = 'u2-toolbar-icon';
      icon.textContent = options.icon;
      el.append(icon);
    }
    const text = document.createElement('span');
    text.className = 'u2-toolbar-label';
    text.textContent = label;
    el.append(text);
    if (options && options.tooltip)
      Tooltip.bind(el, options.tooltip, this.scope);
    return el;
  }

  private _add(item: Item): Toolbar {
    this._items.push(item);
    this.root.insertBefore(item.el, this._chevron);
    this._width = -1;
    this._layout();
    return this;
  }

  private _layout(): void {
    const width = this.root.clientWidth;
    if (width === this._width)
      return;
    this._width = width;
    for (const item of this._items)
      Toolbar._setHidden(item, false);
    this._chevron.classList.add(HIDDEN);
    if (width === 0 || this._items.length === 0) {
      this._roving();
      return;
    }
    const style = getComputedStyle(this.root);
    const gap = parseFloat(style.columnGap) || 0;
    const available = width - (parseFloat(style.paddingLeft) || 0) - (parseFloat(style.paddingRight) || 0);
    const widths = this._items.map((item) => item.el.offsetWidth);
    let total = widths.reduce((a, b) => a + b, 0) + gap * (widths.length - 1);
    if (total <= available) {
      this._roving();
      return;
    }
    this._chevron.classList.remove(HIDDEN);
    const budget = available - this._chevron.offsetWidth - gap;
    let visible = widths.length;
    while (visible > 0 && total > budget) {
      total -= widths[visible - 1] + gap;
      visible--;
    }
    while (visible > 0 && this._items[visible - 1].kind === 'separator')
      visible--;
    for (let i = visible; i < this._items.length; i++)
      Toolbar._setHidden(this._items[i], true);
    this._roving();
  }

  private _roving(): void {
    const stops = this._stops();
    if (this._active >= stops.length)
      this._active = stops.length - 1;
    if (this._active < 0)
      this._active = 0;
    for (const item of this._items) {
      const target = Toolbar._target(item);
      if (target)
        target.tabIndex = -1;
    }
    this._chevron.tabIndex = -1;
    if (stops.length > 0)
      stops[this._active].tabIndex = 0;
  }

  private _stops(): HTMLElement[] {
    const stops: HTMLElement[] = [];
    for (const item of this._items) {
      const target = item.hidden ? null : Toolbar._target(item);
      if (target)
        stops.push(target);
    }
    if (!this._chevron.classList.contains(HIDDEN))
      stops.push(this._chevron);
    return stops;
  }

  private _onFocusIn(e: Event): void {
    const stops = this._stops();
    const index = stops.findIndex((s) => s === e.target || s.contains(e.target as Node));
    if (index < 0)
      return;
    this._active = index;
    this._roving();
  }

  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight')
      return;
    // arrows belong to the caret when the focus sits in a text editor dropped into a slot
    if ((e.target as HTMLElement).matches(CARET))
      return;
    const stops = this._stops();
    if (stops.length === 0)
      return;
    e.preventDefault();
    this._active = (this._active + (e.key === 'ArrowRight' ? 1 : -1) + stops.length) % stops.length;
    this._roving();
    stops[this._active].focus();
  }

  private _toggleMenu(el: HTMLElement, build: (menu: Menu) => void): void {
    const reopen = this._menuOwner !== el;
    this._closeMenu();
    if (!reopen)
      return;
    const menu = new Menu();
    build(menu);
    this._menu = menu;
    this._menuOwner = el;
    menu.show({anchor: el});
    const watch = new Scope();
    watch.effect(() => {
      const open = menu.isOpen.value;
      el.setAttribute('aria-expanded', String(open));
      if (!open && this._menu === menu)
        this._menuOwner = undefined;
    });
    this._menuWatch = watch;
  }

  private _closeMenu(): void {
    if (this._menu)
      this._menu.close();
    if (this._menuWatch)
      this._menuWatch.dispose();
    this._menu = undefined;
    this._menuWatch = undefined;
    this._menuOwner = undefined;
  }

  private _buildOverflow(menu: Menu): void {
    const hidden: Item[] = [];
    for (const item of this._items) {
      if (!item.hidden || item.kind === 'slot')
        continue;
      if (item.kind === 'separator' &&
          (hidden.length === 0 || hidden[hidden.length - 1].kind === 'separator'))
        continue;
      hidden.push(item);
    }
    while (hidden.length > 0 && hidden[hidden.length - 1].kind === 'separator')
      hidden.pop();
    for (const item of hidden) {
      if (item.kind === 'separator')
        menu.separator();
      else if (item.kind === 'button')
        menu.item(item.label, item.onClick!, {icon: item.icon});
      else if (item.kind === 'toggle')
        menu.item(item.label, () => item.bind!.value = !item.bind!.peek(), {check: item.bind!.peek()});
      else {
        menu.group(item.label);
        item.build!(menu);
        menu.endGroup();
      }
    }
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }

  private static _setHidden(item: Item, hidden: boolean): void {
    item.hidden = hidden;
    item.el.classList.toggle(HIDDEN, hidden);
  }

  private static _target(item: Item): HTMLElement | null {
    if (item.kind === 'separator')
      return null;
    return item.kind === 'slot' ? item.el.querySelector<HTMLElement>(FOCUSABLE) : item.el;
  }
}
