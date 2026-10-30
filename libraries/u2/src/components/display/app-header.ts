/* The app home's hero: a big icon beside a bullet list of what the app does, and a "Learn more"
   link under it. With a title, the name and a one-line subtitle stand next to a smaller icon. Each
   item is text, an element, or a link; `description` takes the markdown bullets the platform
   headers were written in. Width response follows the header's own width, so it holds in any host: under
   480px the icon is at most 96px, under 400px the bullets go under the icon. */
import {Control} from '../../core/component.js';
import {div, span, link, Child, Text} from '../../core/elements.js';
import {Input} from '../../core/input-base.js';
import type {IItemsState} from '../../core/widget-like.js';
import {icon} from './icon.js';

export interface AppHeaderLink {
  text: string;
  /** Opens in a new tab. */
  href?: string;
  onClick?: (e: MouseEvent | KeyboardEvent) => void;
}

export type AppHeaderItem = Child | AppHeaderLink;

export interface AppHeaderOptions {
  /** An image URL, a Font Awesome name ('flask'), or a ready element (an svg). */
  icon: string | HTMLElement;
  /** In pixels; 150, or 90 beside a title. */
  iconSize?: number;
  title?: Text;
  subtitle?: Text;
  items?: AppHeaderItem[];
  /** Markdown bullets, one item per line ("- a\n- b"), `[text](url)` a link; after `items`. */
  description?: string;
  /** A URL opens in a new tab; `{onClick}` runs instead (the app's own help). "Learn more" by default. */
  learnMore?: string | Partial<AppHeaderLink>;
  /** A rule under the header. */
  bottomLine?: boolean;
  name?: string;
}

const MARKDOWN_LINK = /\[([^\]]+)\]\(([^)\s]+)\)/g;
const FA_NAME = /^[a-z0-9-]+$/;

export class AppHeader extends Control {
  private readonly _items: HTMLLIElement[] = [];
  private readonly _learnMore: (Partial<AppHeaderLink> & {text: string}) | undefined;
  private _title: HTMLElement | undefined;

  constructor(options: AppHeaderOptions) {
    super();
    this.root.classList.add('u2-app-header');
    this.root.classList.toggle('u2-app-header-named', options.title !== undefined);
    this.root.classList.toggle('u2-app-header-rule', options.bottomLine === true);
    this.root.dataset.u2 = 'app-header';
    this.root.style.setProperty('--u2-app-header-icon-size',
      `${options.iconSize ?? (options.title !== undefined ? 90 : 150)}px`);
    const more = typeof options.learnMore === 'string' ? {href: options.learnMore} : options.learnMore;
    this._learnMore = more === undefined ? undefined : {...more, text: more.text ?? 'Learn more'};
    this.build(() => this.root.append(this._lead(options), this._body(options)));
    if (options.name !== undefined)
      this.name = options.name;
    // a container query would need a container: `container-type` on the root collapses it in a fit-content parent.
    // The layout changes in the next frame: written in the callback it resizes the header, and the observers of
    // what is around it are notified again in the same frame ("ResizeObserver loop completed …")
    let frame = 0;
    const fit = new ResizeObserver(() => {
      cancelAnimationFrame(frame);
      frame = requestAnimationFrame(() => this._fit());
    });
    fit.observe(this.root);
    this.own(() => {
      fit.disconnect();
      cancelAnimationFrame(frame);
    });
  }

  /** Narrow hosts only — a docked panel, a phone: every platform landing at 1366px is wider (Chem MPO 569px). */
  private _fit(): void {
    const width = this.root.clientWidth;
    if (width === 0)
      return;
    this.root.classList.toggle('u2-app-header-compact', width < 480);
    this.root.classList.toggle('u2-app-header-narrow', width < 400);
  }

  /** The base state plus the "Learn more" link: its text, and the page it opens when it opens one. */
  getState(): ReturnType<Control['getState']> {
    const more = this._learnMore;
    return {...super.getState(), ...more === undefined ? {} :
      {extra: {learnMore: more.text, ...more.href === undefined ? {} : {learnMoreHref: more.href}}}};
  }

  getItems(opts?: {all?: boolean}): IItemsState {
    return Input.itemsState(this._items.map((li, index) =>
      ({key: String(index), label: Control.text(li), index, rendered: true})), opts?.all);
  }

  /** The title — read off the element, so a header not in the document still has it. */
  protected accessibleLabel(): string | undefined {
    const title = this._title === undefined ? '' : Control.text(this._title);
    return Control.labelOf(this.root) ?? (title || undefined);
  }

  private _lead(options: AppHeaderOptions): HTMLElement {
    const lead = div([AppHeader._part(AppHeader._icon(options.icon), 'icon')], 'u2-app-header-lead');
    if (options.title === undefined)
      return lead;
    this._title = AppHeader._part(span(options.title), 'title');
    const names = div([this._title], 'u2-app-header-names');
    if (options.subtitle !== undefined)
      names.append(AppHeader._part(span(options.subtitle), 'subtitle'));
    lead.append(names);
    return lead;
  }

  private _body(options: AppHeaderOptions): HTMLElement {
    const list = document.createElement('ul');
    list.className = 'u2-app-header-items';
    for (const item of [...options.items ?? [], ...AppHeader._bullets(options.description ?? '')]) {
      const li = AppHeader._part(document.createElement('li'), 'item');
      li.dataset.u2Item = String(this._items.length);
      li.append(AppHeader._content(item));
      this._items.push(li);
      list.append(li);
    }
    const body = div([list], 'u2-app-header-body');
    const more = this._learnMore;
    if (more === undefined)
      return body;
    const a = AppHeader._part(AppHeader._link(more), 'learn-more');
    if (more.href !== undefined) {
      const mark = icon('external-link-square', {cls: 'u2-app-header-external'});
      mark.setAttribute('aria-hidden', 'true');
      a.append(mark);
    }
    body.append(a);
    return body;
  }

  private static _content(item: AppHeaderItem): Node | string {
    if (typeof item === 'string' || item instanceof HTMLElement)
      return item;
    if (Control.is(item))
      return item.root;
    if (AppHeader._isLink(item))
      return AppHeader._link(item);
    return span(item);
  }

  private static _isLink(x: unknown): x is AppHeaderLink {
    return typeof (x as AppHeaderLink).text === 'string';
  }

  private static _link(spec: Partial<AppHeaderLink> & {text: string}): HTMLAnchorElement {
    const onClick = spec.onClick;
    return link(spec.text, spec.href ?? ((e) => onClick?.(e)));
  }

  /** One item per non-empty line, its bullet marker dropped and `[text](url)` made a link. */
  private static _bullets(markdown: string): HTMLElement[] {
    return markdown.split('\n').map((line) => line.trim().replace(/^[-*+]\s+/, '')).filter((line) => line !== '')
      .map((line) => {
        const host = span('');
        let at = 0;
        for (const m of line.matchAll(MARKDOWN_LINK)) {
          host.append(line.slice(at, m.index), link(m[1], m[2]));
          at = m.index! + m[0].length;
        }
        host.append(line.slice(at));
        return host;
      });
  }

  private static _icon(source: string | HTMLElement): HTMLElement {
    if (typeof source !== 'string')
      return source;
    if (FA_NAME.test(source))
      return icon(source);
    // a background, as the platform header drew it: a shrunk icon resamples the same way, and an image that does not
    // load leaves its square empty instead of the broken-image glyph
    const image = div([]);
    image.style.backgroundImage = `url(${JSON.stringify(source)})`;
    return image;
  }

  private static _part<T extends HTMLElement>(el: T, part: string): T {
    el.classList.add(`u2-app-header-${part}`);
    el.dataset.u2Part = part;
    return el;
  }
}
