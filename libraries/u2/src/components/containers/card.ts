/* Card: a bordered surface with optional header, media, body and footer sections. `clickable`
   makes the whole card a button, `selectable` turns a click into a selection toggle — and a
   `Signal` handed as `selected` is adopted, not copied, so an outside owner keeps driving it. */
import {signal, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {div, span, Child, Text} from '../../core/elements.js';
import {iconOf} from '../display/icon.js';

/* a click that lands on a control of its own — a nested clickable card included — belongs to it,
   not to the card around it */
const INTERACTIVE = 'button, a, input, textarea, select, [contenteditable="true"], [role="button"]';

export interface CardOptions {
  title?: Text;
  subtitle?: Text;
  /** A Font Awesome name, or a ready element (an avatar, a structure thumbnail). */
  icon?: string | HTMLElement;
  actions?: (HTMLElement | Control)[];
  /** An image URL, or a ready element rendered in the media band. */
  media?: string | HTMLElement;
  children?: Child[];
  footer?: Child[];
  clickable?: boolean;
  onClick?: () => void;
  selectable?: boolean;
  selected?: boolean | Signal<boolean>;
}

export class Card extends Control {
  readonly selected: Signal<boolean>;
  /** The body element, so content can be appended after construction. */
  readonly body: HTMLElement;

  constructor(options: CardOptions = {}) {
    super();
    this.selected = options.selected instanceof Signal ? options.selected : signal(!!options.selected);
    this.root.classList.add('u2-card');
    this.root.dataset.u2 = 'card';
    this.body = this.run(() => this._build(options));
    this._wire(options);
  }

  private _build(options: CardOptions): HTMLElement {
    const header = this._header(options);
    if (header)
      this.root.append(header);
    if (options.media !== undefined)
      this.root.append(Card._media(options.media));
    const body = div(options.children, 'u2-card-body');
    this.root.append(body);
    if (options.footer)
      this.root.append(div(options.footer, 'u2-card-footer'));
    return body;
  }

  private _header(options: CardOptions): HTMLElement | undefined {
    const {title, subtitle, icon, actions} = options;
    if (title === undefined && subtitle === undefined && icon === undefined && !actions?.length)
      return undefined;
    const parts: Child[] = [];
    if (icon !== undefined) {
      const glyph = iconOf(icon);
      glyph.classList.add('u2-card-icon');
      parts.push(glyph);
    }
    if (title !== undefined || subtitle !== undefined) {
      const titles: Child[] = [];
      if (title !== undefined)
        titles.push(span(title, 'u2-card-title'));
      if (subtitle !== undefined)
        titles.push(span(subtitle, 'u2-card-subtitle'));
      parts.push(div(titles, 'u2-card-titles'));
    }
    if (actions?.length)
      parts.push(div(actions, 'u2-card-actions'));
    return div(parts, 'u2-card-header');
  }

  private static _media(source: string | HTMLElement): HTMLElement {
    if (typeof source !== 'string')
      return div([source], 'u2-card-media');
    const image = document.createElement('img');
    image.className = 'u2-card-media-img';
    image.src = source;
    image.alt = '';
    return div([image], 'u2-card-media');
  }

  private _wire(options: CardOptions): void {
    const onClick = options.onClick;
    const selectable = !!options.selectable;
    const clickable = !!options.clickable || onClick !== undefined;
    if (!clickable && !selectable)
      return;
    if (clickable)
      this.root.classList.add('u2-card-clickable');
    this.root.setAttribute('role', 'button');
    this.root.tabIndex = 0;

    const activate = () => {
      if (selectable)
        this.selected.value = !this.selected.peek();
      if (onClick)
        onClick();
    };
    this._listen(this.root, 'click', (e) => {
      const inner = (e.target as Element).closest(INTERACTIVE);
      if (inner === null || inner === this.root)
        activate();
    });
    this._listen(this.root, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (e.target !== this.root || (key !== 'Enter' && key !== ' '))
        return;
      e.preventDefault();
      activate();
    });

    if (!selectable)
      return;
    this.root.classList.add('u2-card-selectable');
    this.effect(() => {
      const on = this.selected.value;
      this.root.classList.toggle('u2-card-selected', on);
      this.root.setAttribute('aria-pressed', String(on));
    });
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
