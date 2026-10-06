/* Section: one standalone heading-styled header over a collapsible body — the form-category
   look, unlike Accordion's managed pane set. The chevron sits in the host's left padding gutter
   and appears on hover/focus only, so the title stays flush with form captions; collapsed content
   is hidden, never detached, so its state survives. */
import {signal, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {div, span, Text} from '../../core/elements.js';
import {iconOf} from '../display/icon.js';

export interface SectionOptions {
  title: Text;
  /** A `Signal` handed here is adopted, not copied, so an outside owner keeps driving it. */
  expanded?: boolean | Signal<boolean>;
  /** `false` renders a plain heading + body: no chevron, no click, {@link expanded} inert. */
  collapsible?: boolean;
}

let sectionCount = 0;

export class Section extends Control {
  readonly expanded: Signal<boolean>;
  readonly header: HTMLElement;
  readonly body: HTMLElement;

  constructor(options: SectionOptions) {
    super();
    const id = `u2-section-${++sectionCount}`;
    this.root.classList.add('u2-section');
    this.root.dataset.u2 = 'section';
    this.expanded = options.expanded instanceof Signal ? options.expanded :
      signal(options.expanded !== false);
    const collapsible = options.collapsible !== false;
    const title = this.runInScope(() => span(options.title, 'u2-section-title'));
    this.header = div([title], 'u2-section-header');
    this.header.id = `${id}-header`;
    this.body = div([], 'u2-section-body');
    this.body.id = `${id}-body`;
    this.body.setAttribute('role', 'region');
    this.body.setAttribute('aria-labelledby', this.header.id);
    this.root.append(this.header, this.body);
    if (!collapsible)
      return;

    const chevron = document.createElement('span');
    chevron.className = 'u2-section-chevron';
    chevron.setAttribute('aria-hidden', 'true');
    chevron.append(iconOf('chevron-down'));
    this.header.prepend(chevron);
    this.header.setAttribute('role', 'button');
    this.header.setAttribute('aria-controls', this.body.id);
    this.header.tabIndex = 0;
    const toggle = () => this.expanded.value = !this.expanded.peek();
    this._listen(this.header, 'click', toggle);
    this._listen(this.header, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key !== 'Enter' && key !== ' ')
        return;
      e.preventDefault();
      toggle();
    });
    this.effect(() => {
      const on = this.expanded.value;
      this.header.setAttribute('aria-expanded', String(on));
      this.body.style.display = on ? '' : 'none';
    });
  }

  add(...children: (Control | HTMLElement)[]): Section {
    for (const child of children)
      this.body.append(Control.is(child) ? child.root : child);
    return this;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
