/* Accordion: panes expand independently (the platform's `ui.accordion` has no exclusive mode).
   Lazy content is built once under the pane's own scope, so `removePane` disposes whatever the
   builder constructed; collapsed content is hidden, never detached, so its state survives. */
import {signal, untracked, Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import type {Text} from '../../core/elements.js';
import {iconOf} from '../display/icon.js';
import {Fold} from './fold.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface AccordionOptions {
  /** Titles and content start at the accordion's left edge and the chevron goes to the header's
   * end: panes that carry on a form, their titles and captions in line with the captions above. */
  flush?: boolean;
}

export interface AccordionPaneOptions {
  expanded?: boolean;
  /** A Font Awesome name (rendered through `icon`) or a ready element, shown before the title. */
  icon?: string | HTMLElement;
  /** After the title while the pane is folded: what it holds, or a live summary of it — "ID 0 ·
   * analog C · 5 properties". A signal is followed; an element (a summary with swatches) is placed
   * as it is. */
  hint?: Text | HTMLElement;
  /** Why the content needs the user, followed live: the header turns the error color and, while the
   * pane is folded, the message shows in place of the hint — a fold never hides an error. */
  validity?: ReadonlySignal<string | null>;
  /** Remembers the user's own open/close per browser under this key; the pane's {@link
   * AccordionPane.userSet} says whether there is such a choice, so a host can pick a default (open
   * when there is room) until then. */
  persistKey?: string;
}

export interface AccordionPane {
  readonly title: string;
  readonly icon?: HTMLElement;
  readonly expanded: Signal<boolean>;
  readonly root: HTMLElement;
  /** Whether the user opened or closed it — now, or in an earlier session under its `persistKey`. */
  readonly userSet: boolean;
}

interface Pane extends AccordionPane {
  readonly header: HTMLElement;
  readonly content: HTMLElement;
  readonly scope: Scope;
  readonly persistKey?: string;
  readonly validity?: ReadonlySignal<string | null>;
  userSet: boolean;
  builder?: () => HTMLElement;
}

let accordionCount = 0;

export class Accordion extends Control {
  private readonly _idPrefix = `u2-accordion-${++accordionCount}`;
  private readonly _panes: Pane[] = [];
  private _nextPaneId = 0;

  constructor(options?: AccordionOptions) {
    super();
    this.root.classList.add('u2-accordion');
    this.root.classList.toggle('u2-accordion-flush', options?.flush === true);
    this.root.dataset.u2 = 'accordion';
    this.scope.listen(this.root, 'click', (e) => this._onClick(e));
    this.scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
  }

  /** `options` — or, the older way, `expanded` and then `icon` — as {@link AccordionPaneOptions} says. */
  addPane(title: string, content: HTMLElement | (() => HTMLElement), options?: boolean | AccordionPaneOptions,
    icon?: string | HTMLElement): AccordionPane {
    const opts: AccordionPaneOptions = typeof options === 'object' ? options : {expanded: options, icon};
    icon = opts.icon;
    const stored = Fold.stored('accordion', opts.persistKey);
    const expanded = stored ?? opts.expanded === true;
    const id = `${this._idPrefix}-${this._nextPaneId++}`;
    const root = document.createElement('div');
    root.className = 'u2-accordion-pane';

    const header = document.createElement('div');
    header.className = 'u2-accordion-header';
    header.id = `${id}-header`;
    header.tabIndex = this._panes.length === 0 ? 0 : -1;
    header.title = title;
    header.dataset.u2Item = title;
    header.setAttribute('role', 'button');
    header.setAttribute('aria-controls', `${id}-content`);

    const chevron = document.createElement('span');
    chevron.className = 'u2-accordion-chevron';
    chevron.append(iconOf('chevron-down'));
    header.append(chevron);
    const iconEl = icon === undefined ? undefined : iconOf(icon);
    if (iconEl) {
      const wrap = document.createElement('span');
      wrap.className = 'u2-accordion-icon';
      wrap.append(iconEl);
      header.append(wrap);
    }
    const label = document.createElement('span');
    label.className = 'u2-accordion-title';
    label.dataset.u2Part = 'item-label';
    label.textContent = title;
    header.append(label);
    const hint = opts.hint === undefined && opts.validity === undefined ? undefined : document.createElement('span');
    if (hint !== undefined) {
      hint.className = 'u2-accordion-hint';
      hint.dataset.u2Part = 'pane-hint';
      header.append(hint);
    }

    const host = document.createElement('div');
    host.className = 'u2-accordion-content';
    host.id = `${id}-content`;
    host.setAttribute('role', 'region');
    host.setAttribute('aria-labelledby', header.id);

    root.append(header, host);
    this.root.append(root);

    const scope = this.scope.child('accordion-pane');
    const pane: Pane = {
      title,
      icon: iconEl,
      root,
      header,
      content: host,
      scope,
      expanded: signal(expanded),
      persistKey: opts.persistKey,
      validity: opts.validity,
      userSet: stored !== null,
      builder: typeof content === 'function' ? content : undefined,
    };
    if (!pane.builder)
      host.append(content as HTMLElement);
    this._panes.push(pane);
    scope.effect(() => this._apply(pane));
    if (hint !== undefined) {
      Fold.bindHint(scope, {root, header, hint, cls: 'u2-accordion-pane', text: opts.hint, validity: opts.validity,
        folded: () => !pane.expanded.value, tooltip: () => pane.title});
    }
    return pane;
  }

  getPane(title: string): AccordionPane | undefined {
    return this._panes.find((p) => p.title === title);
  }

  /** Renames the pane at `index` in place — a bound pane title following its source; the title is
   * also the pane's item key. */
  setPaneTitle(index: number, title: string): void {
    const pane = this._panes[index];
    if (pane === undefined)
      return;
    (pane as {title: string}).title = title;
    pane.header.title = title;
    pane.header.dataset.u2Item = title;
    pane.header.querySelector('.u2-accordion-title')!.textContent = title;
  }

  /** By position — unambiguous where {@link getPane} is not (duplicate titles). */
  paneAt(index: number): AccordionPane | undefined {
    return this._panes[index];
  }

  removePane(title: string): void {
    const index = this._panes.findIndex((p) => p.title === title);
    if (index < 0)
      return;
    const pane = this._panes[index];
    pane.scope.dispose();
    pane.root.remove();
    this._panes.splice(index, 1);
    if (this._panes.length > 0 && !this._panes.some((p) => p.header.tabIndex === 0))
      this._panes[0].header.tabIndex = 0;
  }

  /** The panes as items: keyed by title, each saying whether it is expanded. */
  getItems(): IItemsState {
    const list = this._panes.map((pane, index) => ({key: pane.title, label: pane.title, index,
      expanded: pane.expanded.peek(), ...pane.validity?.peek() != null ? {invalid: true} : {}, rendered: true}));
    return {count: list.length, list, truncated: false};
  }

  private _apply(pane: Pane): void {
    const on = pane.expanded.value;
    pane.header.setAttribute('aria-expanded', String(on));
    pane.content.style.display = on ? '' : 'none';
    if (!on || !pane.builder)
      return;
    const build = pane.builder;
    pane.builder = undefined;
    untracked(() => pane.content.append(pane.scope.build(build)));
  }

  private _toggle(index: number): void {
    const pane = this._panes[index];
    pane.expanded.value = !pane.expanded.peek();
    pane.userSet = true;
    Fold.store('accordion', pane.persistKey, pane.expanded.peek());
  }

  private _focus(index: number): void {
    for (let i = 0; i < this._panes.length; i++)
      this._panes[i].header.tabIndex = i === index ? 0 : -1;
    this._panes[index].header.focus();
  }

  private _onClick(e: Event): void {
    const header = (e.target as HTMLElement).closest('.u2-accordion-header') as HTMLElement | null;
    if (!header)
      return;
    const index = this._panes.findIndex((p) => p.header === header);
    if (index < 0)
      return;
    this._focus(index);
    this._toggle(index);
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const current = this._panes.findIndex((p) => p.header === document.activeElement);
    if (current < 0)
      return;
    const last = this._panes.length - 1;
    let next: number;
    switch (e.key) {
      case 'ArrowUp': next = current === 0 ? last : current - 1; break;
      case 'ArrowDown': next = current === last ? 0 : current + 1; break;
      case 'Home': next = 0; break;
      case 'End': next = last; break;
      case 'Enter':
      case ' ':
        e.preventDefault();
        this._toggle(current);
        return;
      default: return;
    }
    e.preventDefault();
    this._focus(next);
  }
}
