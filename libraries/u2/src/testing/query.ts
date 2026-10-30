import {state as controlState} from 'datagrok-api/u2core';
import type {IControlState} from 'datagrok-api/u2core';
import {U2_KINDS} from '../core/kinds.js';

/** Where a finder looks: an element, the document, or a control (its root). */
export type Within = Element | Document | {root: Element};

const OPERABLE = 'button, a[href], input, select, textarea, summary, [contenteditable="true"], [tabindex], ' +
  '[role="button"], [role="link"], [role="checkbox"], [role="radio"], [role="switch"], [role="tab"], ' +
  '[role="menuitem"], [role="menuitemcheckbox"], [role="menuitemradio"], [role="option"], [role="combobox"], ' +
  '[role="textbox"], [role="slider"], [role="spinbutton"], [role="treeitem"]';
const NAMED_BY_TEXT = 'button, a, summary, [role="button"], [role="link"], [role="tab"], [role="menuitem"], ' +
  '[role="menuitemcheckbox"], [role="menuitemradio"], [role="option"], [role="treeitem"]';
const GENERIC_CLASSES = ['u2-div', 'u2-component', 'ui-div'];

/** The element under `within` named `name`: its `data-u2-name` (the automation id a control is given), else the name
 * a user reads on it — `aria-label`, its label, its text for a button-like element, its `title` or placeholder
 * (case aside). Throws when none or several are, listing the names there.
 * @example
 * user.click(byName(dialog, 'Edit Campaign Name'));
 * const status = byName(dialog, 'status'); */
export function byName(within: Within, name: string): HTMLElement {
  const root = Query.root(within);
  const named = Query.all(root, '[data-u2-name]').filter((el) => el.getAttribute('data-u2-name') === name);
  if (named.length === 1)
    return named[0];
  if (named.length > 1) {
    throw new Error(`byName: ${named.length} elements in ${Names.of(root)} are named "${name}": ` +
      named.map((el) => Names.describe(el)).join(', '));
  }
  const labels = Query.labels();
  const read = Query.all(root, OPERABLE).filter((el) => Names.same(Names.accessible(el, labels), name));
  const innermost = read.filter((el) => !read.some((other) => other !== el && el.contains(other)));
  const shown = innermost.filter((el) => Query.hidden(el) === undefined);
  if (innermost.length === 1 || shown.length === 1)
    return innermost.length === 1 ? innermost[0] : shown[0];
  const several = shown.length > 0 ? shown : innermost;
  if (several.length > 1) {
    throw new Error(`byName: ${several.length} elements in ${Names.of(root)} read "${name}": ` +
      several.map((el) => Names.describe(el)).join(', ') + '; give the one you mean a name');
  }
  throw new Error(`byName: nothing in ${Names.of(root)} is named "${name}"; named there: ` +
    `${Query.names(root, labels).join(', ') || 'nothing'}`);
}

/** The control roots of kind `kind` (their `data-u2`) under `within`, in document order; throws when there are none.
 * @example
 * const [grid] = byKind(dialog, 'access-grid'); */
export function byKind(within: Within, kind: string): HTMLElement[] {
  const root = Query.root(within);
  const found = Query.all(root, '[data-u2]').filter((el) => el.getAttribute('data-u2') === kind);
  if (found.length > 0)
    return found;
  const known = U2_KINDS.some((k) => k.kind === kind);
  const there = [...new Set(Query.all(root, '[data-u2]').map((el) => el.getAttribute('data-u2')))];
  throw new Error(`byKind: no ${kind} in ${Names.of(root)}` + (known ? '' : ` (u2 has no kind "${kind}")`) +
    `; kinds there: ${there.join(', ') || 'none'}`);
}

/** The state of the u2 control owning `of` — what automation reads (`DG.U2.state`): label, value, enabled with the
 * reason, readOnly, valid, error, items… Throws outside a u2 control.
 * @example
 * expect(state(byName(dialog, 'Name')).readOnly, false); */
export function state(of: Element | {root: Element}): IControlState {
  const el = Query.element(of);
  const found = controlState(el);
  if (found === null)
    throw new Error(`state: ${Names.describe(el)} is not inside a u2 control`);
  return found;
}

export class Query {
  static root(within: Within): Element | Document {
    return 'root' in within && within.root instanceof Element ? within.root : within as Element | Document;
  }

  static element(of: Element | {root: Element}): HTMLElement {
    return Query.root(of) as HTMLElement;
  }

  static all(root: Element | Document, selector: string): HTMLElement[] {
    const self = root instanceof Element && root.matches(selector) ? [root as HTMLElement] : [];
    return [...self, ...Array.from(root.querySelectorAll<HTMLElement>(selector))];
  }

  static hidden(el: Element): string | undefined {
    const who = () => Names.called(el);
    if (!el.isConnected)
      return `${who()} is not on the page`;
    for (let at: Element | null = el; at !== null; at = at.parentElement) {
      const where = at === el ? '' : ` (in ${Names.describe(at)})`;
      if (at.hasAttribute('hidden'))
        return `${who()} is hidden${where}`;
      if (getComputedStyle(at).display === 'none')
        return `${who()} is not shown: display: none${where}`;
    }
    const visibility = getComputedStyle(el).visibility;
    return visibility === 'hidden' || visibility === 'collapse' ? `${who()} is not shown: visibility: ${visibility}` :
      undefined;
  }

  static labels(): Map<string, string> {
    const labels = new Map<string, string>();
    for (const label of Array.from(document.querySelectorAll('label[for]')))
      labels.set(label.getAttribute('for')!, Names.squash(label.textContent ?? ''));
    return labels;
  }

  static names(root: Element | Document, labels: Map<string, string>): string[] {
    const ids = Query.all(root, '[data-u2-name]').map((el) => el.getAttribute('data-u2-name')!);
    const read = Query.all(root, OPERABLE).map((el) => Names.accessible(el, labels)).filter((n) => n !== '');
    return [...new Set([...ids, ...read.map((n) => `"${n}"`)])].slice(0, 40);
  }
}

export class Names {
  static squash(text: string): string {
    return text.replace(/\s+/g, ' ').trim();
  }

  static same(read: string, wanted: string): boolean {
    return Names.squash(read).toLowerCase() === Names.squash(wanted).toLowerCase();
  }

  static accessible(el: Element, labels?: Map<string, string>): string {
    const aria = Names.squash(el.getAttribute('aria-label') ?? '');
    if (aria !== '')
      return aria;
    const by = (el.getAttribute('aria-labelledby') ?? '').split(/\s+/).filter((id) => id !== '')
      .map((id) => document.getElementById(id)?.textContent ?? '').join(' ');
    const label = Names.squash(by) || (el.id ? labels?.get(el.id) ?? '' : '');
    if (label !== '')
      return label;
    const text = el.matches(NAMED_BY_TEXT) ? Names.squash(el.textContent ?? '') : '';
    return text || Names.squash(el.getAttribute('title') ?? '') || Names.squash(el.getAttribute('placeholder') ?? '');
  }

  static called(el: Element): string {
    const name = Names.accessible(el, Query.labels()).slice(0, 40);
    return name === '' ? Names.describe(el) : `${el.tagName.toLowerCase()} "${name}"`;
  }

  static describe(el: Element): string {
    const classes = typeof el.className === 'string' ? el.className.trim().split(/\s+/).filter((c) => c !== '') : [];
    const own = classes.filter((c) => !GENERIC_CLASSES.includes(c));
    const cls = own[own.length - 1] ?? classes[classes.length - 1] ?? '';
    const named = el.getAttribute('data-u2-name');
    return `${el.tagName.toLowerCase()}${cls ? `.${cls}` : ''}${named ? `[data-u2-name="${named}"]` : ''}`;
  }

  static of(root: Element | Document): string {
    return root instanceof Element ? Names.describe(root) : 'the document';
  }
}
