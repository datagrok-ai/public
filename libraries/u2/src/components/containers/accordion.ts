/* Accordion: panes expand independently (the platform's `ui.accordion` has no exclusive mode).
   Lazy content is built once under the pane's own scope, so `removePane` disposes whatever the
   builder constructed; collapsed content is hidden, never detached, so its state survives. */
import {signal, untracked, Signal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {iconOf} from '../display/icon.js';

export interface AccordionPane {
  readonly title: string;
  readonly icon?: HTMLElement;
  readonly expanded: Signal<boolean>;
  readonly root: HTMLElement;
}

interface Pane extends AccordionPane {
  readonly header: HTMLElement;
  readonly content: HTMLElement;
  readonly scope: Scope;
  builder?: () => HTMLElement;
}

let accordionCount = 0;

export class Accordion extends Control {
  private readonly _idPrefix = `u2-accordion-${++accordionCount}`;
  private readonly _panes: Pane[] = [];
  private _nextPaneId = 0;

  constructor() {
    super();
    this.root.classList.add('u2-accordion');
    this.root.dataset.u2 = 'accordion';
    this._listen(this.root, 'click', (e) => this._onClick(e));
    this._listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
  }

  /** `icon` is a Font Awesome name (rendered through `icon`) or a ready element, shown before the title. */
  addPane(title: string, content: HTMLElement | (() => HTMLElement), expanded = false,
    icon?: string | HTMLElement): AccordionPane {
    const id = `${this._idPrefix}-${this._nextPaneId++}`;
    const root = document.createElement('div');
    root.className = 'u2-accordion-pane';

    const header = document.createElement('div');
    header.className = 'u2-accordion-header';
    header.id = `${id}-header`;
    header.tabIndex = this._panes.length === 0 ? 0 : -1;
    header.title = title;
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
    label.textContent = title;
    header.append(label);

    const host = document.createElement('div');
    host.className = 'u2-accordion-content';
    host.id = `${id}-content`;
    host.setAttribute('role', 'region');
    host.setAttribute('aria-labelledby', header.id);

    root.append(header, host);
    this.root.append(root);

    const scope = new Scope();
    this.own(() => scope.dispose());
    const pane: Pane = {
      title,
      icon: iconEl,
      root,
      header,
      content: host,
      scope,
      expanded: signal(expanded),
      builder: typeof content === 'function' ? content : undefined,
    };
    if (!pane.builder)
      host.append(content as HTMLElement);
    this._panes.push(pane);
    scope.effect(() => this._apply(pane));
    return pane;
  }

  getPane(title: string): AccordionPane | undefined {
    return this._panes.find((p) => p.title === title);
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

  private _apply(pane: Pane): void {
    const on = pane.expanded.value;
    pane.header.setAttribute('aria-expanded', String(on));
    pane.content.style.display = on ? '' : 'none';
    if (!on || !pane.builder)
      return;
    const build = pane.builder;
    pane.builder = undefined;
    untracked(() => pane.content.append(Scope.runWith(pane.scope, build)));
  }

  private _toggle(index: number): void {
    const expanded = this._panes[index].expanded;
    expanded.value = !expanded.peek();
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

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
