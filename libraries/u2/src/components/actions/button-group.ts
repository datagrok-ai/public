/* Segmented button group: attached buttons that act as an action group, a radio group or a set of
   independent toggles, in the three densities that cover forms, toolbars and ribbons. Icon groups
   are this control in `iconOnly` mode, not a separate one. WAI-ARIA: a single tab stop with roving
   ←/→ focus; `single` is radiogroup/radio, everything else stays plain buttons. */
import {signal, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Tooltip} from '../../core/tooltip.js';
import {icon} from '../display/icon.js';
import type {IItemsState} from '../../core/widget-like.js';

export type ButtonGroupToggle = 'none' | 'single' | 'multi';
export type ButtonGroupDensity = 'form' | 'toolbar' | 'ribbon';

export interface ButtonGroupItem {
  id: string;
  label?: string;
  icon?: string;
  tooltip?: string;
  /** Gets the click and the group. */
  onClick?: (e: MouseEvent, group: ButtonGroup) => void;
  enabled?: boolean;
}

export interface ButtonGroupOptions {
  items: ButtonGroupItem[];
  toggle?: ButtonGroupToggle;
  /** Labels move into the tooltip and the accessible name. */
  iconOnly?: boolean;
  density?: ButtonGroupDensity;
}

export class ButtonGroup extends Control {
  /** Selected ids. `single` keeps 0..1 and never deselects on a repeat click (radio semantics);
   * `none` stays empty. Writes always produce a fresh array. */
  readonly selected: Signal<string[]> = signal<string[]>([]);

  private readonly _items: ButtonGroupItem[];
  private readonly _buttons = new Map<string, HTMLButtonElement>();
  private readonly _toggle: ButtonGroupToggle;
  private _active: string | undefined;

  constructor(options: ButtonGroupOptions) {
    super();
    this._items = options.items;
    this._toggle = options.toggle ?? 'none';
    const iconOnly = options.iconOnly === true;

    this.root.classList.add('u2-btn-group', `u2-btn-group-${options.density ?? 'form'}`);
    if (iconOnly)
      this.root.classList.add('u2-btn-group-icons');
    this.root.dataset.u2 = 'button-group';
    this.root.setAttribute('role', this._toggle === 'single' ? 'radiogroup' : 'group');
    for (const item of this._items)
      this.root.append(this._build(item, iconOnly));

    this.scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this.scope.listen(this.root, 'focusin', (e) => this._onFocusIn(e));
    this.effect(() => this._applySelection());
    this._roving();
  }

  setEnabled(id: string, enabled: boolean): void {
    const el = this._buttons.get(id);
    if (!el)
      return;
    el.disabled = !enabled;
    this._roving();
  }

  /** The base state plus the selection as the value (the id under `single`) and the current id. */
  getState(): ReturnType<Control['getState']> {
    const selected = this.selected.peek();
    return {...super.getState(), value: this._toggle === 'single' ? selected[0] ?? null : selected,
      ...selected.length === 0 ? {} : {current: selected[0]}};
  }

  /** The buttons, keyed by id; under a toggle mode each says whether it is selected. */
  getItems(): IItemsState {
    const selected = this.selected.peek();
    const list = this._items.map((item, index) => ({key: item.id, label: item.label ?? item.tooltip ?? item.id, index,
      ...this._toggle === 'none' ? {} : {selected: selected.includes(item.id)},
      disabled: this._buttons.get(item.id)!.disabled, rendered: true}));
    return {count: list.length, list, truncated: false};
  }

  private _build(item: ButtonGroupItem, iconOnly: boolean): HTMLButtonElement {
    const el = document.createElement('button');
    el.type = 'button';
    el.className = 'u2-btn-group-item';
    el.dataset.id = item.id;
    el.dataset.u2Item = item.id;
    el.tabIndex = -1;
    if (item.icon)
      el.append(icon(item.icon));
    if (item.label && !iconOnly) {
      const label = document.createElement('span');
      label.className = 'u2-btn-group-label';
      label.dataset.u2Part = 'item-label';
      label.textContent = item.label;
      el.append(label);
    }
    const tooltip = item.tooltip ?? (iconOnly ? item.label : undefined);
    if (tooltip) {
      Tooltip.bind(el, tooltip, this.scope);
      if (iconOnly)
        el.setAttribute('aria-label', tooltip);
    }
    if (this._toggle === 'single') {
      el.setAttribute('role', 'radio');
      el.setAttribute('aria-checked', 'false');
    } else if (this._toggle === 'multi')
      el.setAttribute('aria-pressed', 'false');
    if (item.enabled === false)
      el.disabled = true;
    this.scope.listen(el, 'click', (e) => this._activate(item, e as MouseEvent));
    this._buttons.set(item.id, el);
    return el;
  }

  private _activate(item: ButtonGroupItem, e: MouseEvent): void {
    if (this._toggle === 'single')
      this.selected.value = [item.id];
    else if (this._toggle === 'multi') {
      const current = this.selected.peek();
      this.selected.value = current.includes(item.id) ?
        current.filter((id) => id !== item.id) : current.concat([item.id]);
    }
    if (item.onClick)
      item.onClick(e, this);
  }

  private _applySelection(): void {
    const selected = new Set(this.selected.value);
    for (const item of this._items) {
      const el = this._buttons.get(item.id)!;
      const on = selected.has(item.id);
      el.classList.toggle('u2-btn-group-selected', on);
      if (this._toggle === 'single')
        el.setAttribute('aria-checked', String(on));
      else if (this._toggle === 'multi')
        el.setAttribute('aria-pressed', String(on));
    }
  }

  private _stops(): HTMLButtonElement[] {
    const stops: HTMLButtonElement[] = [];
    for (const item of this._items) {
      const el = this._buttons.get(item.id)!;
      if (!el.disabled)
        stops.push(el);
    }
    return stops;
  }

  // the tab stop is tracked by id, not by position: disabling or enabling a cell must not slide
  // the stop onto a different button
  private _roving(): void {
    const stops = this._stops();
    for (const el of this._buttons.values())
      el.tabIndex = -1;
    if (stops.length === 0)
      return;
    const active = stops.find((el) => el.dataset.id === this._active) ?? stops[0];
    this._active = active.dataset.id;
    active.tabIndex = 0;
  }

  private _onFocusIn(e: Event): void {
    const el = e.target as HTMLButtonElement;
    if (!el.dataset || !this._buttons.has(el.dataset.id as string))
      return;
    this._active = el.dataset.id;
    this._roving();
  }

  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight')
      return;
    const stops = this._stops();
    if (stops.length === 0)
      return;
    e.preventDefault();
    const index = Math.max(stops.findIndex((el) => el.dataset.id === this._active), 0);
    const next = stops[(index + (e.key === 'ArrowRight' ? 1 : -1) + stops.length) % stops.length];
    this._active = next.dataset.id;
    this._roving();
    next.focus();
  }
}
