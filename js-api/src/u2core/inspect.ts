/* The in-page accessors automation reads u2 controls through (`DG.U2.state(el)` …): resolve the
   control owning an element, answer its status, state block, items and declared properties. They
   resolve u2 controls only — a platform `DG.Widget` answers `DG.Widget.find(el)?.getWidgetStatus()`. */
import {Control} from './component.js';
import {untracked} from './signals.js';
import type {IControlState, IItemsState, IWidgetStatus, JsonValue} from './protocol.js';

/** One declared property of a control, as {@link props} reports it. */
export interface PropState {
  /** The property name `setProp` takes. */
  name: string;
  /** The property type (`string`, `bool`, `int` …); `object` when undeclared. */
  type: string;
  /** The current value, JSON-projected. */
  value: JsonValue | undefined;
  /** Whether `setProp` can write it. */
  writable: boolean;
}

/** The u2 control owning `el`: the nearest control root at or above it, or the owner of the popup it
 * sits in. Platform widget wrappers are passed over — whether one exists depends on who asked the
 * platform for it. Undefined for an element outside any u2 control.
 * @example
 * const input = DG.U2.control(document.querySelector('[data-u2="text-input"]'));
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function control(el: Element | null): Control | undefined {
  for (let found = Control.forElement(el); found !== undefined; found = Control.forElement(found.root.parentElement)) {
    if (typeof (found as {toDart?: unknown}).toDart !== 'function')
      return found;
  }
  return undefined;
}

/** Full status of the control owning `el`; its parts are elements, so use it inside the page only.
 * @example
 * const parts = DG.U2.status(el)?.parts;
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function status(el: Element | null): IWidgetStatus | null {
  return control(el)?.getWidgetStatus() ?? null;
}

/** JSON-safe state of the control owning `el` — what `page.evaluate` hands to a test runner; null
 * outside a u2 control. Synchronous and side-effect free.
 * @example
 * const {label, value, valid, enabled} = DG.U2.state(el)!;
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function state(el: Element | null): IControlState | null {
  return control(el)?.getState() ?? null;
}

/** Items of the control owning `el` (the state's `items`), uncapped with `all: true`; null when it
 * lists none.
 * @example
 * const keys = DG.U2.items(el, {all: true})?.list.map((item) => item.key);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function items(el: Element | null, opts?: {all?: boolean}): IItemsState | null {
  return control(el)?.getItems(opts) ?? null;
}

/** Brings an item of the control owning `el` into view by key or label (see `Control.revealItem`);
 * rejects naming the candidates when a label matches several.
 * @example
 * const row = await DG.U2.revealItem(el, 'Discovery');
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export async function revealItem(el: Element | null, keyOrLabel: string): Promise<HTMLElement | null> {
  return await control(el)?.revealItem(keyOrLabel) ?? null;
}

/** The declared properties of the control owning `el`, values JSON-projected.
 * @example
 * const writable = DG.U2.props(el).filter((p) => p.writable).map((p) => p.name);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function props(el: Element | null): PropState[] {
  const c = control(el);
  if (c === undefined)
    return [];
  return c.getProperties().map((p) => ({name: p.name, type: p.type ?? p.propertyType ?? 'object',
    value: Control.jsonValue(p.get?.(c.propertyTarget)), writable: p.set != null}));
}

/** A rule {@link operability} checks: `covered` — at its centre the page answers another element; `no-pointer` —
 * `pointer-events: none` on it or above it in its control that no `:hover` of an ancestor taking the pointer lifts;
 * `inert` — inside an `inert` region that says no reason;
 * `disabled-without-reason` — off with no reason on it, its control or a control around it that is off too;
 * `name-mismatch` — an icon-only button whose
 * accessible name is not its tooltip; `oversized-option` — an option taller than the one-row rail it sits on. */
export type OperabilityRule = 'covered' | 'no-pointer' | 'inert' | 'disabled-without-reason' | 'name-mismatch' |
  'oversized-option';

/** A control the page offers that a user cannot operate as offered, as {@link operability} finds it. */
export interface OperabilityFinding {
  /** The rule it breaks. */
  rule: OperabilityRule;
  /** The element (a control root for `disabled-without-reason`, the region for `inert`); in the page only. */
  element: Element;
  /** What a person reads: `button "Edit" is covered by div.u2-input-box`. */
  detail: string;
}

/** Options of {@link operability}. */
export interface OperabilityOptions {
  /** Also checks `covered`, with `document.elementFromPoint` at each element's centre: it needs layout (a browser). */
  hitTest?: boolean;
  /** More selectors of what a user operates (a kind's `interactive` parts). */
  interactive?: string[];
}

const OPERABLE = 'button, a[href], input, select, textarea, summary, [contenteditable], [tabindex], [role="button"], ' +
  '[role="link"], [role="checkbox"], [role="radio"], [role="switch"], [role="tab"], [role="menuitem"], ' +
  '[role="menuitemcheckbox"], [role="menuitemradio"], [role="option"], [role="combobox"], [role="textbox"], ' +
  '[role="slider"], [role="spinbutton"], [role="treeitem"], [data-u2-part="options"] > *';
const TOOLTIP = Symbol.for('u2.tooltip');
const CONTROL = Symbol.for('u2.control');
// the containers' generic classes: a finding names an element by a class of its own
const GENERIC_CLASSES = ['u2-div', 'u2-component', 'ui-div'];

const HOVER = /:hover(?![\w-])/g;

/** One selector of a page rule that sets `pointer-events`, as it reads with the pointer over an ancestor. */
interface PointerRule {
  /** The selector with `:hover` taken off. */
  selector: string;
  /** What an element at or inside the hovered ancestor matches: the selector through its last `:hover` compound,
   * `:hover` off; null when nothing in it is hovered. */
  hovered: string | null;
  none: boolean;
}

/** What under `root` a user is offered and cannot operate as offered: every visible interactive element (native
 * controls, `[tabindex]`, interactive roles, `interactive`, the children of an options rail) is checked against the
 * {@link OperabilityRule}s. An `aria-hidden` element is not offered; what a platform widget below `root` holds is the
 * platform's to judge; what a control holds inside another that is off and says why (a disabled filter builder, a
 * designer canvas in Design mode) is off with that reason; an element under a sticky part of its own scroller is
 * scrolled out of view, not covered; one
 * the page's style rules reveal on hover (a row's actions: `.row:hover .actions {pointer-events: auto}`) takes the
 * pointer once it is over the row, so it is no finding and its hit test at rest is skipped.
 * Platform-free; without layout (Node) `covered` is skipped and rectangles are whatever the test set.
 * @example
 * for (const f of DG.U2.operability(dialog.root, {hitTest: true}))
 *   console.log(f.rule, f.detail);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-operability} */
export function operability(root: Element, options?: OperabilityOptions): OperabilityFinding[] {
  const selector = [OPERABLE, ...options?.interactive ?? []].join(', ');
  const candidates = [...root.matches(selector) ? [root] : [], ...root.querySelectorAll(selector)];
  const findings: OperabilityFinding[] = [];
  let pointer: PointerRule[] | undefined;
  const reported = new Map<Element, Set<OperabilityRule>>();
  const report = (rule: OperabilityRule, element: Element, detail: string) => {
    const rules = reported.get(element) ?? new Set<OperabilityRule>();
    if (!rules.has(rule)) {
      reported.set(element, rules.add(rule));
      findings.push({rule, element, detail});
    }
  };
  for (const el of new Set(candidates)) {
    if (el.getAttribute('contenteditable') === 'false' || !Probe.visible(el) || Probe.hosted(el, root))
      continue;
    const region = el.closest('[inert]');
    if (region !== null) {
      if (Probe.reason(region) === '')
        report('inert', region, `${Probe.describe(region)} is inert and says no reason: ${Probe.called(el)} does nothing`);
      continue;
    }
    const control = Control.forElement(el)?.root;
    const boundary = control !== undefined && control.contains(el) ? control : root;
    const off = Probe.offAt(el, boundary);
    const said = Probe.offSaying(boundary, root);
    if (off !== undefined) {
      if (Probe.reason(off) === '' && Probe.reason(boundary) === '' && !said)
        report('disabled-without-reason', off, `${Probe.called(off)} is disabled and says no reason`);
      continue;
    }
    if (said)
      continue;
    const source = Probe.noPointer(el, boundary);
    if (source !== undefined && !Probe.revealed(el, source, pointer ??= Probe.pointerRules())) {
      report('no-pointer', el, `${Probe.called(el)} takes no pointer: pointer-events: none` +
        (source === el ? '' : ` on ${Probe.describe(source)}`));
      continue;
    }
    if (options?.hitTest && source === undefined)
      Probe.covered(el, report);
    Probe.mismatch(el, report);
    Probe.oversized(el, report);
  }
  return findings;
}

class Probe {
  static visible(el: Element): boolean {
    if (el.closest('[hidden], [aria-hidden="true"]') !== null || el.getClientRects().length === 0)
      return false;
    const style = getComputedStyle(el);
    return style.display !== 'none' && style.visibility !== 'hidden' && style.visibility !== 'collapse';
  }

  // the nearest owner mark below root decides: a u2 control, or a platform widget (`data-widget`, or its wrapper)
  static hosted(el: Element, root: Element): boolean {
    for (let at: Element | null = el; at !== null && at !== root; at = at.parentElement) {
      const owner = (at as unknown as Record<symbol, {toDart?: unknown} | undefined>)[CONTROL];
      if (owner !== undefined)
        return typeof owner.toDart === 'function';
      if (at.hasAttribute('data-widget'))
        return true;
    }
    return false;
  }

  static offAt(el: Element, boundary: Element): Element | undefined {
    let off: Element | undefined;
    for (let at: Element | null = el; at !== null; at = at === boundary ? null : at.parentElement) {
      if (at.getAttribute('aria-disabled') === 'true' || at.hasAttribute('disabled'))
        off = at;
    }
    return off;
  }

  // above el's own control, up to the probed root: a control nested in one that is off and says why is off with it
  static offSaying(boundary: Element, root: Element): boolean {
    if (boundary === root || !root.contains(boundary))
      return false;
    for (let at = boundary.parentElement; at !== null; at = at === root ? null : at.parentElement) {
      if ((at.getAttribute('aria-disabled') === 'true' || at.hasAttribute('disabled')) && Probe.reason(at) !== '')
        return true;
    }
    return false;
  }

  static reason(el: Element): string {
    return Probe.squash(el.getAttribute('aria-description') ?? '') || Probe.squash(el.getAttribute('title') ?? '');
  }

  static tooltip(el: Element): string {
    let content = (el as unknown as Record<symbol, unknown>)[TOOLTIP];
    try {
      if (typeof content === 'function')
        content = untracked(content as () => unknown);
    }
    catch (_) {
      content = '';
    }
    const text = content instanceof Element ? content.textContent : typeof content === 'string' ? content : '';
    return Probe.squash(text || el.getAttribute('title') || '');
  }

  static noPointer(el: Element, boundary: Element): Element | undefined {
    let source: Element | undefined;
    for (let at: Element | null = el; at !== null; at = at === boundary ? null : at.parentElement) {
      const value = getComputedStyle(at).pointerEvents;
      if (value === 'none')
        source = at;
      else if (value)
        break;
    }
    return source;
  }

  // the none is lifted on the way in: a rule applies to el, or to a part above it holding no none of its own, while
  // an ancestor of `source` is hovered — those take the pointer, `source` and what is inside it do not
  static revealed(el: Element, source: Element, rules: PointerRule[]): boolean {
    const above = source.parentElement;
    for (let at: Element | null = el; at !== null; at = at.parentElement) {
      const part = at;
      if ((part as HTMLElement).style?.pointerEvents === 'none')
        return false;
      if (above !== null && rules.some((r) => !r.none && r.hovered !== null && Probe.matches(part, r.selector) &&
          Probe.matches(above, r.hovered)))
        return true;
      if (part === source || rules.some((r) => r.none && Probe.matches(part, r.selector)))
        return false;
    }
    return false;
  }

  static pointerRules(): PointerRule[] {
    const rules: PointerRule[] = [];
    const visit = (list: CSSRuleList) => {
      for (const rule of Array.from(list)) {
        const {selectorText, style} = rule as CSSStyleRule;
        const value = typeof selectorText === 'string' ? style?.pointerEvents : '';
        for (const selector of value ? Probe.split(selectorText, ',') : []) {
          const parts = Probe.split(selector, ' ');
          const last = parts.map((part) => Probe.mask(part).search(HOVER) >= 0).lastIndexOf(true);
          const head = Probe.unhover(parts.slice(0, last + 1));
          rules.push({selector: Probe.unhover(parts), hovered: last < 0 ? null : `${head}, ${head} *`,
            none: value === 'none'});
        }
        if ((rule as CSSGroupingRule).cssRules !== undefined)
          visit((rule as CSSGroupingRule).cssRules);
      }
    };
    for (const sheet of [...Array.from(document.styleSheets ?? []), ...document.adoptedStyleSheets ?? []]) {
      try {
        visit(sheet.cssRules);
      }
      catch (_) {
        // a cross-origin sheet keeps its rules to itself
      }
    }
    return rules;
  }

  // a selector the page cannot match (a nested rule's `&`, a pseudo-class it does not know) matches nothing
  static matches(el: Element, selector: string): boolean {
    try {
      return el.matches(selector);
    }
    catch (_) {
      return false;
    }
  }

  // `text` split where `separator` stands outside brackets and parentheses (` ` stands for any whitespace)
  static split(text: string, separator: ',' | ' '): string[] {
    const masked = Probe.mask(text);
    const parts: string[] = [];
    let from = 0;
    for (let i = 0; i <= masked.length; i++) {
      if (i === masked.length || (separator === ' ' ? /\s/.test(masked[i]) : masked[i] === separator)) {
        const part = text.slice(from, i).trim();
        if (part !== '')
          parts.push(part);
        from = i + 1;
      }
    }
    return parts;
  }

  // what brackets and parentheses hold blanked out, at the same offsets: a `:hover` or a comma inside `:not(…)` stays
  static mask(text: string): string {
    let depth = 0;
    let masked = '';
    for (let i = 0; i < text.length; i++) {
      if (text[i] === ')' || text[i] === ']')
        depth--;
      masked += depth > 0 ? '_' : text[i];
      if (text[i] === '(' || text[i] === '[')
        depth++;
    }
    return masked;
  }

  static unhover(parts: string[]): string {
    return parts.map((part) => {
      let bare = '';
      let from = 0;
      for (const m of Probe.mask(part).matchAll(HOVER)) {
        bare += part.slice(from, m.index);
        from = m.index! + m[0].length;
      }
      return bare + part.slice(from) || '*';
    }).join(' ');
  }

  static covered(el: Element, report: (rule: OperabilityRule, el: Element, detail: string) => void): void {
    const r = el.getBoundingClientRect();
    const x = r.left + r.width / 2;
    const y = r.top + r.height / 2;
    if (r.width === 0 || r.height === 0 || x < 0 || y < 0 || x >= innerWidth || y >= innerHeight)
      return;
    const clips: Element[] = [];
    for (let at = el.parentElement; at !== null; at = at.parentElement) {
      const style = getComputedStyle(at);
      if ((style.overflowX || 'visible') === 'visible' && (style.overflowY || 'visible') === 'visible')
        continue;
      const box = at.getBoundingClientRect();
      if (x < box.left || x >= box.right || y < box.top || y >= box.bottom)
        return;
      clips.push(at);
    }
    const hit = document.elementFromPoint(x, y);
    if (hit === null || hit === el || el.contains(hit) || hit.closest('label')?.control === el)
      return;
    // a sticky part of el's own scroller is its edge: el is scrolled under it, as past the clip edge
    const stuck = Probe.stuck(hit, el);
    if (stuck !== undefined && clips.some((clip) => clip.contains(stuck)))
      return;
    report('covered', el, `${Probe.called(el)} is covered by ${Probe.describe(hit)}`);
  }

  static stuck(hit: Element, el: Element): Element | undefined {
    for (let at: Element | null = hit; at !== null && !at.contains(el); at = at.parentElement) {
      if (getComputedStyle(at).position === 'sticky')
        return at;
    }
    return undefined;
  }

  static mismatch(el: Element, report: (rule: OperabilityRule, el: Element, detail: string) => void): void {
    if (!el.matches('button, [role="button"], a') || Probe.squash(el.textContent ?? '') !== '')
      return;
    const tooltip = Probe.tooltip(el);
    const name = Probe.squash(el.getAttribute('aria-label') ?? '');
    if (tooltip === '' || name === tooltip)
      return;
    report('name-mismatch', el, name === '' ? `${Probe.describe(el)} has no accessible name; its tooltip says ` +
      `"${tooltip}"` : `${Probe.describe(el)} is "${name}" to a screen reader and "${tooltip}" on hover`);
  }

  static oversized(el: Element, report: (rule: OperabilityRule, el: Element, detail: string) => void): void {
    const rail = el.parentElement;
    if (rail?.getAttribute('data-u2-part') !== 'options')
      return;
    const row = rail.getBoundingClientRect().height;
    const height = el.getBoundingClientRect().height;
    if (row > 0 && height > row + 1) {
      report('oversized-option', el, `${Probe.describe(el)} is ${Math.round(height)}px tall on a ${Math.round(row)}px ` +
        'options rail: an element shown as the value belongs in a DynamicInput row');
    }
  }

  static called(el: Element): string {
    const label = Probe.squash(el.getAttribute('aria-label') ?? '') || Probe.squash(el.getAttribute('title') ?? '') ||
      Probe.squash(el.textContent ?? '').slice(0, 40);
    return label === '' ? Probe.describe(el) : `${el.tagName.toLowerCase()} "${label}"`;
  }

  static describe(el: Element): string {
    const classes = typeof el.className === 'string' ? el.className.trim().split(/\s+/) : [];
    const own = classes.filter((c) => !GENERIC_CLASSES.includes(c));
    const cls = own[own.length - 1] ?? classes[classes.length - 1] ?? '';
    const named = el.getAttribute('data-u2-name');
    return `${el.tagName.toLowerCase()}${cls ? `.${cls}` : ''}${named ? `[data-u2-name="${named}"]` : ''}`;
  }

  static squash(text: string): string {
    return text.replace(/\s+/g, ' ').trim();
  }
}

/** Writes a declared property of the control owning `el` through its setter — for putting a control in
 * a state from a test. Throws when the control has no such property or it is read-only.
 * @example
 * DG.U2.setProp(el, 'enabled', false);
 * @see {@link https://public.datagrok.ai/js/samples/ui/widgets/u2-status} */
export function setProp(el: Element | null, name: string, value: unknown): void {
  const c = control(el);
  if (c === undefined)
    throw new Error('u2: not inside a u2 control');
  const kind = c.root.dataset.u2 ?? c.componentMeta?.tag ?? 'control';
  const prop = c.getProperties().find((p) => p.name === name);
  if (prop === undefined)
    throw new Error(`u2: ${kind} has no property ${name}`);
  if (prop.set == null)
    throw new Error(`u2: ${kind}.${name} is read-only`);
  prop.set(c.propertyTarget, value);
}
