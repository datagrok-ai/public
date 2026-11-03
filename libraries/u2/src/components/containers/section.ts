/* Section: one standalone heading-styled header over a collapsible body — the form-category
   look, unlike Accordion's managed pane set. The chevron sits in the host's left padding gutter
   and appears on hover/focus only, so the title stays flush with form captions; collapsed content
   is hidden, never detached, so its state survives. */
import {signal, Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {div, span, Text} from '../../core/elements.js';
import {Reason} from '../../core/reason.js';
import {iconOf} from '../display/icon.js';
import {Fold} from './fold.js';

export interface SectionOptions {
  title: Text;
  /** A `Signal` handed here is adopted, not copied, so an outside owner keeps driving it. */
  expanded?: boolean | Signal<boolean>;
  /** `false` renders a plain heading + body: no chevron, no click, {@link expanded} inert. */
  collapsible?: boolean;
  /** A line after the title saying what the body holds — "ID, colors, properties", or a live summary
   * of it — so a folded section still tells what is inside. An element (a summary with swatches) is
   * placed as it is. */
  hint?: Text | HTMLElement;
  /** Why the body needs the user, followed live: the header turns the error color and, while the
   * section is folded, the hint shows the message — an error is never hidden by a fold. */
  validity?: ReadonlySignal<string | null>;
  /** Remembers the user's own open/close per browser under this key; {@link Section.userSet} says
   * whether there is such a choice, so a host can pick a default (open when there is room) until then. */
  persistKey?: string;
  /** `'end'`: the chevron after the title and hint, always shown, pointing down while folded — the look
   * of a toggle placed in a row of its own. The {@link Section.header} stays the section's toggle wherever
   * it is put: `row.append(section.header)` sets Details beside a Save, its body opening below. */
  chevron?: 'start' | 'end';
  /** After the title: the icons and links that act on the body (a table's columns and refresh); a press on one
   * does not fold the section. */
  actions?: (Control | HTMLElement)[];
  /** Why the body is off: dimmed and out of reach, the reason said after the title; `null` while on. */
  disabled?: Reason;
}

let sectionCount = 0;

export class Section extends Control {
  readonly expanded: Signal<boolean>;
  readonly header: HTMLElement;
  readonly body: HTMLElement;

  private readonly _title: HTMLElement;
  private readonly _collapsible: boolean;
  private readonly _validity: ReadonlySignal<string | null> | undefined;
  private _userSet = false;

  constructor(options: SectionOptions) {
    super();
    const id = `u2-section-${++sectionCount}`;
    this.root.classList.add('u2-section');
    this.root.dataset.u2 = 'section';
    this.root.setAttribute('role', 'group');
    const stored = Fold.stored('section', options.persistKey);
    this._userSet = stored !== null;
    this.expanded = options.expanded instanceof Signal ? options.expanded :
      signal(stored ?? options.expanded !== false);
    if (stored !== null && options.expanded instanceof Signal)
      this.expanded.value = stored;
    this._validity = options.validity;
    const collapsible = this._collapsible = options.collapsible !== false;
    const title = this._title = this.build(() => span(options.title, 'u2-section-title'));
    title.id = `${id}-title`;
    title.dataset.u2Part = 'title';
    this.root.setAttribute('aria-labelledby', title.id);
    // only asked for: a section's content may have a hint part of its own
    const hint = options.hint === undefined && options.validity === undefined ? undefined :
      span('', 'u2-section-hint');
    if (hint !== undefined)
      hint.dataset.u2Part = 'section-hint';
    this.header = div(hint === undefined ? [title] : [title, hint], 'u2-section-header');
    this.header.id = `${id}-header`;
    this.header.dataset.u2Part = 'header';
    this.body = div([], 'u2-section-body');
    this.body.id = `${id}-body`;
    this.body.dataset.u2Part = 'body';
    this.body.setAttribute('role', 'region');
    this.body.setAttribute('aria-labelledby', this.header.id);
    this.root.append(this.header, this.body);
    if (hint !== undefined) {
      Fold.bindHint(this.scope, {root: this.root, header: this.header, hint, cls: 'u2-section', text: options.hint,
        validity: options.validity, folded: () => collapsible && !this.expanded.value});
    }
    if (options.disabled !== undefined) {
      const reason = span('', 'u2-section-reason');
      reason.dataset.u2Part = 'reason';
      this.header.append(reason);
      Reason.follow(options.disabled, (why) => {
        reason.textContent = why ?? '';
        Reason.region(this.body, why);
      }, this.scope);
    }
    if (options.actions !== undefined) {
      const actions = span('', 'u2-section-actions');
      actions.dataset.u2Part = 'actions';
      for (const action of options.actions)
        actions.append(Control.is(action) ? action.root : action);
      // a press or a key on an action is the action's, not the header's toggle
      for (const type of ['click', 'keydown'])
        this.scope.listen(actions, type, (e) => e.stopPropagation());
      this.header.append(actions);
    }
    if (!collapsible)
      return;

    const chevron = document.createElement('span');
    chevron.className = 'u2-section-chevron';
    chevron.setAttribute('aria-hidden', 'true');
    chevron.append(iconOf('chevron-down'));
    if (options.chevron === 'end') {
      this.header.classList.add('u2-section-chevron-end');
      this.header.append(chevron);
    }
    else
      this.header.prepend(chevron);
    this.header.setAttribute('role', 'button');
    this.header.setAttribute('aria-controls', this.body.id);
    this.header.tabIndex = 0;
    const toggle = () => {
      this.expanded.value = !this.expanded.peek();
      this._userSet = true;
      Fold.store('section', options.persistKey, this.expanded.peek());
    };
    this.scope.listen(this.header, 'click', toggle);
    this.scope.listen(this.header, 'keydown', (e) => {
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

  /** Whether the user opened or closed it — now, or in an earlier session under its `persistKey`. */
  get userSet(): boolean {
    return this._userSet;
  }

  /** The base state plus whether the body shows (a plain heading's always does) and, with a
   * `validity`, whether it says something. */
  getState(): ReturnType<Control['getState']> {
    const invalid = this._validity?.peek() ?? null;
    return {...super.getState(), expanded: !this._collapsible || this.expanded.peek(),
      ...invalid !== null ? {valid: false, error: invalid} : {}};
  }

  /** The title — read off the element, so a section not in the document still has it. */
  protected accessibleLabel(): string | undefined {
    return this.root.getAttribute('aria-label') ?? (Control.text(this._title) || undefined);
  }
}
