/* Noun plan → Playwright Locator. Registry entries are plain selectors; kinds try their qualifier
   strategies in order and keep the first that finds something (the union of all of them when
   nothing does, so negative assertions still have a locator to fail against); scopes nest
   locators, with the u2 owner edge as the fallback for popups portaled out of their owner. Inside
   a context, the context's own names and generic kinds are looked up in its root first, then on
   the whole page (dialogs and notifications are portaled out of it); global names are global. */
import type {Locator, Page} from '@playwright/test';
import {contextFirst, describeNoun, NounRef, parseNoun} from '../nouns.js';
import type {KindEntry} from '../registry.js';
import {contextOf, ElementRef} from './args.js';

export {describeNoun, parseNoun};

type Base = Page | Locator;

const lastRefs = new WeakMap<Page, NounRef>();

/** The phrase parsed with the page's active context — the one parse every runtime path shares.
 * Remembered per page: the phrase a step failed on is the last one it parsed. */
export function refOf(page: Page, target: ElementRef | string): NounRef {
  const ref = parseNoun(typeof target === 'string' ? target : target.phrase, contextOf(page));
  lastRefs.set(page, ref);
  return ref;
}

const LABEL_MAX = 40;
const SHOWN_MAX = 20;

/** What the page shows where the last phrase looked, for a failure report: the scope that is not
 * open, or the visible elements of the phrase's kind (labelled as the platform labels them) so the
 * author sees the name to write instead. Empty when nothing was looked up. */
export async function explain(page: Page): Promise<string> {
  const ref = lastRefs.get(page);
  if (!ref)
    return '';
  let base: Base = page;
  if (ref.scope) {
    const scope = await locateRef(page, ref.scope);
    if (await scope.filter({visible: true}).count() === 0)
      return `${ref.scope.raw}: not open`;
    base = scope.first();
  }
  const plan = ref.plan;
  const selector = plan.type === 'kind' ? [plan, ...plan.alternatives].map((s) => s.kind.selector).join(', ') :
    plan.type === 'entry' ? plan.entry.selector : plan.selector;
  const labels = await base.locator(selector).evaluateAll((els, max) => {
    const isShown = (e: Element) => e.getClientRects().length > 0 && getComputedStyle(e).visibility !== 'hidden';
    const shown = els.filter(isShown);
    const label = (e: Element) => e.getAttribute('data-u2-name') || e.getAttribute('aria-label') || e.getAttribute('title') ||
      ((e as HTMLElement).innerText ?? e.textContent ?? '').trim().split('\n')[0].trim() || e.tagName.toLowerCase();
    return {total: els.length, visible: shown.length, labels: [...new Set(shown.map((e) => label(e).slice(0, max)))],
      hidden: [...new Set(els.filter((e) => !isShown(e)).map((e) => (e.getAttribute('name') || label(e)).slice(0, max)))]};
  }, LABEL_MAX);
  // an element that is in the DOM but not shown is a different finding from one that is not there
  const hidden = labels.hidden.length > 0 ? ` (in the DOM but not shown: ${labels.hidden.slice(0, SHOWN_MAX).join(' | ')})` : '';
  const where = ref.scope ? ` in ${ref.scope.raw}` : '';
  const what = plan.type === 'kind' ? `${plan.kind.name}s` : plan.type === 'entry' ? `"${plan.entry.name}"` : `"${plan.part}"`;
  // the context panel renders the current object: what that is says whether the click landed
  const object = ref.scope && /context panel|property panel/i.test(ref.scope.raw) ? `; ${await currentObject(page)}` : '';
  if (labels.visible === 0)
    return `${what}${where}: ${labels.total === 0 ? 'none on the page' : `${labels.total} present, none visible`}${hidden}${object}`;
  const list = labels.labels.slice(0, SHOWN_MAX).join(' | ') + (labels.labels.length > SHOWN_MAX ? ` | … ${labels.labels.length - SHOWN_MAX} more` : '');
  return `visible ${what}${where}: ${list}${hidden}${object}`;
}

/** `grok.shell.o` as one phrase: `current object: Project "BDD-CP-Root"`. */
export function currentObject(page: Page): Promise<string> {
  return page.evaluate(() => {
    const o = (window as any).grok?.shell?.o;
    if (o == null)
      return 'current object: none';
    const name = o.friendlyName ?? o.name ?? o.caption ?? '';
    return `current object: ${o.constructor?.name ?? typeof o}${name ? ` "${name}"` : ''}`;
  }).catch(() => 'current object: unreadable');
}

export async function locate(page: Page, target: ElementRef | string, within?: Locator): Promise<Locator> {
  const ref = refOf(page, target);
  const loc = await locateRef(page, ref, within);
  return loc.describe(describeNoun(ref));
}

/** The element a gesture or a state check acts on: the visible matches of the phrase (a Dart menu
 * keeps a zero-size mirror of every item under "Properties..."; a closed view leaves its viewers
 * behind). Several visible matches stay ambiguous — Playwright's strict mode reports them; an
 * ordinal names its element as counted, visible or not. No roundtrip of its own. */
export async function locateActionable(page: Page, target: ElementRef | string, within?: Locator): Promise<Locator> {
  const ref = refOf(page, target);
  const loc = await locateRef(page, ref, within);
  return (ref.ordinal === undefined ? loc.filter({visible: true}) : loc).describe(describeNoun(ref));
}

export async function locateRef(page: Page, ref: NounRef, within?: Locator): Promise<Locator> {
  let base: Base = within ?? page;
  let scope: Locator | undefined;
  if (ref.scope) {
    scope = await locateRef(page, ref.scope, within);
    base = scope;
  }
  else if (ref.plan.type === 'entry' && ref.plan.entry.in && !within)
    base = await locate(page, ref.plan.entry.in);
  else if (!within && contextFirst(ref)) {
    const inRoot = await inBase(page, page.locator(ref.context!.selector), ref);
    if (await inRoot.count() > 0)
      return pick(inRoot, ref);
  }
  let loc = await inBase(page, base, ref);
  // a scope that is not on the page has no owner edge to try — and `getAttribute` on it would
  // wait the whole action timeout for it to appear
  if (scope && await loc.count() === 0 && await scope.count() > 0) {
    const owner = await scope.first().getAttribute('data-u2-name').catch(() => null);
    if (owner) {
      const alt = await inBase(page, page.locator(`[data-u2-owner="${cssString(owner)}"]`), ref);
      if (await alt.count() > 0)
        loc = alt;
    }
  }
  return pick(loc, ref);
}

function pick(loc: Locator, ref: NounRef): Locator {
  if (ref.ordinal === 'last')
    return loc.last();
  return ref.ordinal === undefined ? loc : loc.nth(ref.ordinal);
}

async function inBase(page: Page, base: Base, ref: NounRef): Promise<Locator> {
  const plan = ref.plan;
  if (plan.type === 'entry') {
    const loc = base.locator(plan.entry.selector);
    return plan.entry.text === undefined ? loc : loc.filter({hasText: plan.entry.text});
  }
  if (plan.type === 'part')
    return base.locator(plan.selector);
  const candidates: Locator[] = [];
  for (const split of [plan, ...plan.alternatives])
    candidates.push(...(split.qualifier ? strategies(page, base, split.kind, split.qualifier) : [base.locator(split.kind.selector)]));
  for (const c of candidates) {
    if (await c.count() > 0)
      return c;
  }
  return candidates.reduce((a, b) => a.or(b));
}

/** The element's whole text is the qualifier, allowing decoration around it (an icon glyph, a
 * trailing colon). */
export function exactText(q: string): RegExp {
  return new RegExp(`^\\W*${escapeRegExp(q)}\\W*$`, 'i');
}

function strategies(page: Page, base: Base, kind: KindEntry, q: string): Locator[] {
  const compact = q.replace(/\s+/g, '');
  // Dart names join a menu path with "---": "Markers > Size" is div-Markers---Size
  const dashed = q.replace(/\s*>\s*/g, '---').replace(/\s+/g, '-');
  const out: Locator[] = [];
  for (const strategy of kind.match) {
    switch (strategy) {
      case 'name':
        out.push(base.locator(withAttr(kind.selector, `[data-u2-name="${cssString(compact)}" i]`)));
        if (dashed !== compact)
          out.push(base.locator(withAttr(kind.selector, `[data-u2-name="${cssString(dashed)}" i]`)));
        break;
      case 'label':
      case 'title':
        if (kind.labelSelector)
          out.push(byLabel(page, base, kind, q));
        break;
      case 'text':
        out.push(base.locator(kind.selector).filter({hasText: exactText(q)}));
        break;
      case 'exact-text':
        out.push(base.getByText(exactText(q)));
        break;
      case 'aria':
        out.push(base.locator(withAttr(kind.selector, `[aria-label="${cssString(q)}" i]`)));
        out.push(base.locator(withAttr(kind.selector, `[title="${cssString(q)}" i]`)));
        break;
      case 'placeholder':
        out.push(base.locator(kind.selector).filter({has: page.locator(`[placeholder="${cssString(q)}" i]`)}));
        break;
      case 'dart':
        for (const template of kind.dartNames ?? []) {
          out.push(base.locator(withAttr(kind.selector, `[name="${cssString(template.replace('{q}', dashed))}" i]`)));
          if (dashed !== q)
            out.push(base.locator(withAttr(kind.selector, `[name="${cssString(template.replace('{q}', q))}" i]`)));
        }
        break;
    }
  }
  return out;
}

/** The elements of the kind whose label reads `q`. A label that is a direct child (`:scope > …`)
 * is found first and its parent taken: on a popup of a few hundred menu items a `has:` filter
 * over every item costs ~35 ms per query, the label's parent ~2 ms. */
function byLabel(page: Page, base: Base, kind: KindEntry, q: string): Locator {
  const label = kind.labelSelector!;
  const direct = label.split(',').map((s) => s.trim());
  if (direct.every((s) => s.startsWith(':scope >')))
    return base.locator(direct.map((s) => s.slice(':scope >'.length).trim()).join(', '), {hasText: exactText(q)})
      .locator('xpath=parent::*').and(page.locator(kind.selector));
  return base.locator(kind.selector).filter({has: page.locator(label, {hasText: exactText(q)})});
}

/** Appends an attribute selector to every alternative of a comma-separated selector list. */
export function withAttr(selector: string, attr: string): string {
  return selector.split(',').map((s) => s.trim() + attr).join(', ');
}

export function cssString(s: string): string {
  return s.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
}

export function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
