/* Card: a bordered surface with optional header, media, body and footer sections. `clickable`
   makes the whole card a button, `selectable` turns a click into a selection — a toggle, or with
   `selectMode: 'select'` a pick that never unselects — and a `Signal` handed as `selected` is
   adopted, not copied, so an outside owner keeps driving it; a card given `selected` shows it
   whether or not a click can change it. */
import {signal, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {debug} from '../../core/scope.js';
import {div, span, Child, Text} from '../../core/elements.js';
import {iconOf} from '../display/icon.js';

/* a click that lands on a control of its own — a nested clickable card included — belongs to it,
   not to the card around it */
const INTERACTIVE = 'button, a, input, textarea, select, [contenteditable="true"], [role="button"]';

export type CardSelectMode = 'toggle' | 'select';

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
  /** Gets the click — or the Enter / Space keydown that activated the card — and the card. */
  onClick?: (e: MouseEvent | KeyboardEvent, card: Card) => void;
  selectable?: boolean;
  /** What a click does to {@link Card.selected}: `toggle` (default) flips it, `select` sets it. */
  selectMode?: CardSelectMode;
  selected?: boolean | Signal<boolean>;
  /** Tighter padding and gaps, a smaller title: for dense grids of cards. */
  compact?: boolean;
}

export class Card extends Control {
  readonly selected: Signal<boolean>;
  /** The body element, so content can be appended after construction. */
  readonly body: HTMLElement;

  private _title: HTMLElement | undefined;

  constructor(options: CardOptions = {}) {
    super();
    this.selected = options.selected instanceof Signal ? options.selected : signal(!!options.selected);
    this.root.classList.add('u2-card');
    this.root.classList.toggle('u2-card-compact', options.compact === true);
    this.root.dataset.u2 = 'card';
    this.body = this.build(() => this._build(options));
    this._wire(options);
  }

  /** The base state plus the title and the selection. */
  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), selected: this.selected.peek()};
  }

  /** The title — read off the element, so a card not in the document still has it. */
  protected accessibleLabel(): string | undefined {
    const title = this._title === undefined ? '' : Control.text(this._title);
    return Control.labelOf(this.root) ?? (title || undefined);
  }

  private _build(options: CardOptions): HTMLElement {
    const header = this._header(options);
    if (header)
      this.root.append(header);
    if (options.media !== undefined)
      this.root.append(Card._media(options.media));
    const body = Card._part(div(options.children, 'u2-card-body'), 'body');
    this.root.append(body);
    if (options.footer)
      this.root.append(Card._part(div(options.footer, 'u2-card-footer'), 'footer'));
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
        titles.push(this._title = Card._part(span(title, 'u2-card-title'), 'title'));
      if (subtitle !== undefined)
        titles.push(Card._part(span(subtitle, 'u2-card-subtitle'), 'subtitle'));
      parts.push(div(titles, 'u2-card-titles'));
    }
    if (actions?.length)
      parts.push(Card._part(div(actions, 'u2-card-actions'), 'actions'));
    return Card._part(div(parts, 'u2-card-header'), 'header');
  }

  private static _media(source: string | HTMLElement): HTMLElement {
    if (typeof source !== 'string')
      return Card._part(div([source], 'u2-card-media'), 'media');
    const image = document.createElement('img');
    image.className = 'u2-card-media-img';
    image.src = source;
    image.alt = '';
    return Card._part(div([image], 'u2-card-media'), 'media');
  }

  private static _part<T extends HTMLElement>(el: T, name: string): T {
    el.dataset.u2Part = name;
    return el;
  }

  private _wire(options: CardOptions): void {
    const onClick = options.onClick;
    const selectable = !!options.selectable;
    const clickable = !!options.clickable || onClick !== undefined;
    const interactive = clickable || selectable;
    if (selectable || options.selected !== undefined)
      this._showSelection(interactive);
    if (!interactive)
      return;
    if (clickable)
      this.root.classList.add('u2-card-clickable');
    if (selectable)
      this.root.classList.add('u2-card-selectable');
    this.root.setAttribute('role', 'button');
    this.root.tabIndex = 0;

    const toggle = options.selectMode !== 'select';
    const activate = (e: MouseEvent | KeyboardEvent) => {
      if (selectable)
        this.selected.value = toggle ? !this.selected.peek() : true;
      onClick?.(e, this);
    };
    const scope = this.scope;
    scope.listen(this.root, 'click', (e) => {
      const inner = (e.target as Element).closest(INTERACTIVE);
      if (inner === null || inner === this.root)
        activate(e as MouseEvent);
    });
    scope.listen(this.root, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (e.target !== this.root || (key !== 'Enter' && key !== ' '))
        return;
      e.preventDefault();
      activate(e as KeyboardEvent);
    });
    // a card disposed with its owner (a dialog's scope) can stay on the page and stop answering
    if (debug.enabled) {
      this.root.addEventListener('click', () => {
        if (scope.isDisposed && this.root.isConnected)
          console.warn(`u2: a click on card "${this.accessibleLabel() ?? ''}", which was disposed with its owner`);
      });
    }
  }

  /** A card a click can change says so with `aria-pressed`; one only driven from outside, with
   * `aria-current`. */
  private _showSelection(interactive: boolean): void {
    this.effect(() => {
      const on = this.selected.value;
      this.root.classList.toggle('u2-card-selected', on);
      if (interactive)
        this.root.setAttribute('aria-pressed', String(on));
      else if (on)
        this.root.setAttribute('aria-current', 'true');
      else
        this.root.removeAttribute('aria-current');
    });
  }
}
