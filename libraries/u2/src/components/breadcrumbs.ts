/* Breadcrumbs. Overflow is resolved by measurement, not by a wrap listener: the row is laid out
   in full, and if it overflows the middle segments are swapped for one `…` segment that expands
   the path inline when clicked (no menu, so the control stays overlay-free). */
import {signal, Signal, ReadonlySignal} from '../core/signals.js';
import {Control} from '../core/component.js';

export class Breadcrumbs extends Control {
  readonly items: ReadonlySignal<string[]>;

  private readonly _items: Signal<string[]>;
  private readonly _onClick?: (index: number, item: string) => void;
  private _expanded = false;
  private _middle: HTMLElement[] = [];
  private _ellipsis: HTMLElement[] = [];

  constructor(options: {items: string[], onClick?: (index: number, item: string) => void}) {
    super(document.createElement('nav'));
    this._items = signal(options.items.slice());
    this.items = this._items;
    this._onClick = options.onClick;

    this.root.classList.add('u2-breadcrumbs');
    this.root.dataset.u2 = 'breadcrumbs';
    this.root.setAttribute('aria-label', 'breadcrumb');

    this._listen(this.root, 'click', (e) => this._onRootClick(e));
    const resize = new ResizeObserver(() => this._layout());
    resize.observe(this.root);
    this.own(() => resize.disconnect());
    this.effect(() => this._render());
  }

  setItems(items: string[]): void {
    this._expanded = false;
    this._items.value = items.slice();
  }

  private _render(): void {
    const items = this._items.value;
    const last = items.length - 1;
    const nodes: HTMLElement[] = [];
    this._middle = [];
    this._ellipsis = [];
    for (let i = 0; i <= last; i++) {
      const separator = i > 0 ? this._separator() : undefined;
      if (separator)
        nodes.push(separator);
      const item = i === last ? this._current(items[i]) : this._segment(items[i], i);
      nodes.push(item);
      if (separator && i < last)
        this._middle.push(separator, item);
      if (i === 0 && last > 1) {
        const ellipsis = this._segment('…', -1);
        ellipsis.classList.add('u2-breadcrumbs-ellipsis');
        ellipsis.title = 'Show the full path';
        ellipsis.setAttribute('aria-label', 'Show the full path');
        this._ellipsis.push(this._separator(), ellipsis);
        nodes.push(...this._ellipsis);
      }
    }
    this.root.replaceChildren(...nodes);
    this._layout();
  }

  private _layout(): void {
    if (this._ellipsis.length === 0)
      return;
    this._collapse(false);
    if (!this._expanded && this.root.scrollWidth > this.root.clientWidth + 1)
      this._collapse(true);
  }

  private _collapse(on: boolean): void {
    for (const el of this._middle)
      el.style.display = on ? 'none' : '';
    for (const el of this._ellipsis)
      el.style.display = on ? '' : 'none';
  }

  private _segment(text: string, index: number): HTMLButtonElement {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-breadcrumbs-item';
    el.dataset.index = String(index);
    el.textContent = text;
    return el;
  }

  private _current(text: string): HTMLElement {
    const el = document.createElement('span');
    el.className = 'u2-breadcrumbs-current';
    el.textContent = text;
    el.setAttribute('aria-current', 'page');
    return el;
  }

  private _separator(): HTMLElement {
    const el = document.createElement('span');
    el.className = 'u2-breadcrumbs-sep';
    el.textContent = '›';
    el.setAttribute('aria-hidden', 'true');
    return el;
  }

  private _onRootClick(e: Event): void {
    const el = (e.target as HTMLElement).closest('.u2-breadcrumbs-item') as HTMLElement | null;
    if (!el)
      return;
    const index = Number(el.dataset.index);
    if (index < 0) {
      this._expanded = true;
      this._collapse(false);
    } else if (this._onClick)
      this._onClick(index, this._items.peek()[index]);
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
