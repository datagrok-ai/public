import {Control} from '../../core/component.js';
import {div, h1, Text} from '../../core/elements.js';
import {LabelColumn} from '../forms/form.js';

type Child = Control | HTMLElement;

export interface PageSection {
  /** The heading over the content. */
  title?: Text;
  /** Beside the title: the icons and links that act on the content (a table's columns, group and refresh). */
  actions?: Child[];
  content: Child | Child[];
  /** The content's widest, in px: a form that should not stretch across a wide view. */
  maxWidth?: number;
  /** The section's automation id (`data-u2-name`). */
  name?: string;
}

export interface PageOptions {
  /** Above the body and not scrolled with it: an app's header, a view's own toolbar. */
  header?: Child;
  /** The body: sections one under the other, each a heading with its actions over its content. */
  sections?: PageSection[];
  /** The body's content when it has no sections — a form; after the sections when it has both. */
  content?: Child | Child[];
  /** Under the body, pinned while the body scrolls: the action row, its buttons on the right, a `statusLine` taking
   * the room on their left. `form.addButtons(build, {host: page.footer})` puts a form's button row there. */
  footer?: Child[];
  /** The body's widest, in px (`--u2-page-max-width`). */
  maxWidth?: number;
  /** The forms on the page line up on one label column, this one when given; `false` leaves each its own. */
  labelColumn?: LabelColumn | false;
  name?: string;
}

const list = (x: Child | Child[] | undefined): Child[] => x === undefined ? [] : Array.isArray(x) ? x : [x];

/** A view's content: a header, a body that alone scrolls, and an action row pinned under it. Spacing and widths
 * come from the token scale and the `--u2-page-*` properties, so the page needs no CSS of its own. */
export class Page extends Control {
  /** What scrolls: the sections, then the content. */
  readonly body: HTMLElement;
  /** The label column the page's forms share; none with `labelColumn: false`. */
  readonly labelColumn: LabelColumn | undefined;

  private _header: HTMLElement | undefined;
  private _footer: HTMLElement | undefined;

  constructor(options: PageOptions = {}) {
    super();
    const root = this.root;
    root.classList.add('u2-page');
    root.dataset.u2 = 'page';
    if (options.name !== undefined)
      this.name = options.name;
    if (options.maxWidth !== undefined)
      root.style.setProperty('--u2-page-max-width', `${options.maxWidth}px`);
    if (options.header !== undefined) {
      this._header = this._part('header', list(options.header));
      root.append(this._header);
    }
    this.body = this._part('body', [...(options.sections ?? []).map((s) => this._section(s)),
      ...list(options.content)]);
    root.append(this.body);
    if (options.footer !== undefined)
      this.footer.append(...options.footer.map((c) => this._adopt(c)));
    if (options.labelColumn !== false) {
      this.labelColumn = options.labelColumn ?? new LabelColumn();
      LabelColumn.share(root, this.labelColumn);
    }
  }

  get header(): HTMLElement | undefined {
    return this._header;
  }

  /** The pinned action row: the `footer` given, else made on first ask. */
  get footer(): HTMLElement {
    if (this._footer === undefined) {
      this._footer = this._part('footer', []);
      this.root.append(this._footer);
    }
    return this._footer;
  }

  private _part(part: string, children: Child[]): HTMLElement {
    const el = div(children.map((c) => this._adopt(c)), `u2-page-${part}`);
    el.dataset.u2Part = part;
    return el;
  }

  private _section(section: PageSection): HTMLElement {
    const parts: HTMLElement[] = [];
    if (section.title !== undefined || section.actions !== undefined) {
      const title = this.build(() => h1(section.title ?? ''));
      title.className = 'u2-page-section-title';
      title.dataset.u2Part = 'section-title';
      const header = div([title, ...(section.actions ?? []).map((c) => this._adopt(c))], 'u2-page-section-header');
      header.dataset.u2Part = 'section-header';
      parts.push(header);
    }
    const content = div(list(section.content).map((c) => this._adopt(c)), 'u2-page-section-content');
    content.dataset.u2Part = 'section-content';
    if (section.maxWidth !== undefined)
      content.style.maxWidth = `${section.maxWidth}px`;
    parts.push(content);
    const el = div(parts, 'u2-page-section');
    el.dataset.u2Part = 'section';
    if (section.name !== undefined)
      el.dataset.u2Name = section.name;
    return el;
  }

  private _adopt(child: Child): HTMLElement {
    if (!Control.is(child))
      return child;
    this.own(() => child.dispose());
    return child.root;
  }
}
