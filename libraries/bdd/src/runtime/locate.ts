/* Noun plan → Playwright Locator. Registry entries are plain selectors; kinds try their qualifier
   strategies in the R9 order (`effectiveOrder`: what the user sees first, `data-u2-name` after it,
   Dart names where they were) and keep the first that finds something (the union of all of them when
   nothing does, so negative assertions still have a locator to fail against); scopes nest
   locators, with the u2 owner edge as the fallback for popups portaled out of their owner. Inside
   a context, the context's own names and generic kinds are looked up in its root first, then on
   the whole page (dialogs and notifications are portaled out of it); global names are global. */
import {appendFileSync} from 'node:fs';
import {type Locator, type Page, test} from '@playwright/test';
import {contextFirst, describeNoun, KindSplit, NounRef, parseNoun} from '../nouns.js';
import type {KindEntry, MatchStrategy} from '../registry.js';
import {contextOf, ElementRef} from './args.js';
import * as guide from './guide.js';

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
 * open, or the visible elements of the phrase's kind, named as the user sees them
 * (`describeElements`), so the author sees the name to write instead — and the two matches of an
 * ambiguous phrase are told apart by their sections. Empty when nothing was looked up. */
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
  const found = await base.locator(selector).evaluateAll(describeElements, LABEL_MAX);
  // an element that is in the DOM but not shown is a different finding from one that is not there
  const hidden = found.hidden.length > 0 ? ` (in the DOM but not shown: ${found.hidden.slice(0, SHOWN_MAX).join(' | ')})` : '';
  const where = ref.scope ? ` in ${ref.scope.raw}` : '';
  const what = plan.type === 'kind' ? `${plan.kind.name}s` : plan.type === 'entry' ? `"${plan.entry.name}"` : `"${plan.part}"`;
  // the context panel renders the current object: what that is says whether the click landed
  const object = ref.scope && /context panel|property panel/i.test(ref.scope.raw) ? `; ${await currentObject(page)}` : '';
  if (found.visible === 0)
    return `${what}${where}: ${found.total === 0 ? 'none on the page' : `${found.total} present, none visible`}${hidden}${object}`;
  return `visible ${what}${where}: ${listed(found.labels)}${hidden}${object}`;
}

/** The visible elements of `loc` named as the user sees them, one line, `|`-separated:
 * `"Name" in section "Patient" | "Name" in section "Doctor" [name=doctorName]` — what an ambiguity
 * report lists. One roundtrip. */
export async function describeCandidates(loc: Locator): Promise<string> {
  return listed((await loc.evaluateAll(describeElements, LABEL_MAX)).labels);
}

function listed(labels: string[]): string {
  return labels.slice(0, SHOWN_MAX).join(' | ') + (labels.length > SHOWN_MAX ? ` | … ${labels.length - SHOWN_MAX} more` : '');
}

/** In the page. A shown element is `"<label>"` (the u2 status label of a control root, else its label
 * part, aria-label, data-u2-name, first line of text, title), ` in <container> "<title>"` for the nearest
 * u2 section, dialog, form, tabs, card or accordion around it, and ` [name=<data-u2-name>]` when the
 * label is another word; identical descriptions are counted. One that is in the DOM but not shown goes
 * by its `name` or label. */
const describeElements = (els: Element[], max: number) => {
  const U2 = (window as any).DG?.U2;
  // bound by destructuring, not by name: tsx wraps a named function in a `__name` call the page lacks
  const [status, label, container, isShown]:
    [(e: Element) => string, (e: Element) => string, (e: Element) => string, (e: Element) => boolean] = [
    (e) => {
      try {
        return U2?.control?.(e)?.root === e ? U2.state(e)?.label ?? '' : '';
      }
      catch {
        return '';
      }
    },
    (e) => (status(e) || e.querySelector(':scope > [data-u2-part="label"]')?.textContent?.trim() ||
      e.getAttribute('aria-label') || e.getAttribute('data-u2-name') ||
      ((e as HTMLElement).innerText ?? e.textContent ?? '').trim().split('\n')[0].trim() || e.getAttribute('title') ||
      e.tagName.toLowerCase()).slice(0, max),
    (e) => {
      const c = e.parentElement?.closest(['section', 'dialog', 'form', 'object-form', 'func-form', 'tabs', 'card', 'accordion']
        .map((k) => `[data-u2="${k}"]`).join(', '));
      const title = c ? status(c) || c.querySelector('[data-u2-part="title"], .u2-section-title, .u2-dialog-title-text, .u2-card-title')
        ?.textContent?.trim() || c.getAttribute('data-u2-name') : '';
      return c && title ? ` in ${c.getAttribute('data-u2')!.replace(/-/g, ' ')} "${title.slice(0, max)}"` : '';
    },
    (e) => e.getClientRects().length > 0 && getComputedStyle(e).visibility !== 'hidden',
  ];
  const shown = els.filter(isShown);
  const counts = new Map<string, number>();
  for (const e of shown) {
    const l = label(e);
    const name = e.getAttribute('data-u2-name');
    const d = `"${l}"${container(e)}${name && name !== l ? ` [name=${name}]` : ''}`;
    counts.set(d, (counts.get(d) ?? 0) + 1);
  }
  return {total: els.length, visible: shown.length, labels: [...counts].map(([d, n]) => n > 1 ? `${d} ×${n}` : d),
    hidden: [...new Set(els.filter((e) => !isShown(e)).map((e) => (e.getAttribute('name') || label(e)).slice(0, max)))]};
};

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
  const loc = (await locateRef(page, ref, within)).describe(describeNoun(ref));
  await guide.located(page, loc);
  return loc;
}

/** The element a gesture or a state check acts on: the visible matches of the phrase (a Dart menu
 * keeps a zero-size mirror of every item under "Properties..."; a closed view leaves its viewers
 * behind). Several visible matches stay ambiguous — Playwright's strict mode reports them, and the
 * failure lists them (`explain`); an ordinal counts the visible ones too. An item of a u2 collection
 * that is not rendered is revealed first (`revealed`); `locate` never reveals. No roundtrip of its
 * own otherwise. */
export async function locateActionable(page: Page, target: ElementRef | string, within?: Locator): Promise<Locator> {
  const ref = refOf(page, target);
  let loc = await locateRef(page, ref, within);
  if (await revealed(page, ref, loc, within))
    loc = await locateRef(page, ref, within);
  const actionable = (ref.ordinal === undefined ? loc.filter({visible: true}) : loc).describe(describeNoun(ref));
  await guide.located(page, actionable);
  return actionable;
}

/** An item of a u2 collection that is not rendered (a virtual list's row scrolled away, a collapsed
 * tree's node, a closed popup's option) is brought into view by u2 (`DG.U2.revealItem`: key first,
 * then label) when the phrase names its collection and finds no visible match. Whether u2 was asked;
 * u2's rejection ("2 items match …") is the step's failure. No-op where u2 has no `revealItem` and on
 * a collection that is not u2. */
async function revealed(page: Page, ref: NounRef, loc: Locator, within?: Locator): Promise<boolean> {
  const plan = ref.plan;
  if (plan.type !== 'kind' || !ref.scope)
    return false;
  const item = [plan, ...plan.alternatives].find((s) => s.qualifier && s.kind.items);
  if (!item || await (ref.ordinal === undefined ? loc.filter({visible: true}) : loc).count() > 0)
    return false;
  const scope = await locateRef(page, ref.scope, within);
  return scope.evaluateAll(async (els, q) => {
    const U2 = (window as any).DG?.U2;
    if (!U2?.revealItem || !els[0]?.closest('[data-u2]'))
      return false;
    await U2.revealItem(els[0], q);
    return true;
  }, asWritten(ref.raw, item.qualifier));
}

/** The qualifier as the feature spells it — parsing lowercases it, and u2 matches an item key exactly. */
function asWritten(raw: string, q: string): string {
  const i = raw.toLowerCase().indexOf(q);
  const s = i < 0 ? q : raw.slice(i, i + q.length);
  return s.toLowerCase() === q ? s : q;
}

/** Of the visible matches, the one a pointer can reach: a visible match can lie under another element (a tree
 * row repeated in a group the platform lays out under its parent's row, as Favorites under My stuff), and a
 * gesture aimed at it lands on what covers it. One reachable match is the element; none or several leave the
 * matches as they are, still ambiguous. */
export async function reachable(page: Page, target: ElementRef | string): Promise<Locator> {
  const loc = await locateActionable(page, target);
  const hits: number[] = await loc.evaluateAll((els) => els.map((el, i) => {
    const r = el.getBoundingClientRect();
    return [0.5, 0.25, 0.75, 0.1].some((f) => {
      const hit = document.elementFromPoint(r.left + r.width * f, r.top + r.height / 2);
      return hit !== null && (hit === el || el.contains(hit));
    }) ? i : -1;
  }).filter((i) => i >= 0));
  return hits.length === 1 ? loc.nth(hits[0]) : loc;
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
  // a scope read while its container rebuilds can resolve to its label (a section's header before
  // its pane is back): the target is then looked for in every candidate the scope has
  if (scope && ref.scope!.ordinal === undefined && await loc.count() === 0 && await isKindLabel(scope, ref.scope!)) {
    scope = (await candidates(page, await scopeBase(page, ref.scope!, within), ref.scope!)).reduce((a, b) => a.or(b));
    loc = await inBase(page, scope, ref);
  }
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

/** An ordinal counts the visible matches: the Home page keeps viewers of its own in the page,
 * hidden, under a table view. */
function pick(loc: Locator, ref: NounRef): Locator {
  if (ref.ordinal === undefined)
    return loc;
  const shown = loc.filter({visible: true});
  return ref.ordinal === 'last' ? shown.last() : shown.nth(ref.ordinal);
}

async function inBase(page: Page, base: Base, ref: NounRef): Promise<Locator> {
  const plan = ref.plan;
  if (plan.type === 'entry') {
    const loc = base.locator(plan.entry.selector);
    return plan.entry.text === undefined ? loc : loc.filter({hasText: plan.entry.text});
  }
  if (plan.type === 'part')
    return base.locator(plan.selector);
  const all = kindCandidates(page, base, plan);
  const hit = await firstHit(all);
  if (process.env.BDD_RESOLVE_TRACE)
    await trace(page, base, ref, plan, hit);
  return hit?.loc ?? all.map((c) => c.loc).reduce((a, b) => a.or(b));
}

async function candidates(page: Page, base: Base, ref: NounRef): Promise<Locator[]> {
  const plan = ref.plan;
  return plan.type === 'kind' ? kindCandidates(page, base, plan).map((c) => c.loc) : [await inBase(page, base, ref)];
}

interface Candidate {
  loc: Locator;
  kind: KindEntry;
  qualifier: string;
  /** `kind`: an unqualified split, every element of the kind. */
  strategy: MatchStrategy | 'kind';
}

type KindPlan = KindSplit & {alternatives: KindSplit[]};

/** Every split, longest kind first, each in its strategy order: the R9 order, or (`legacy`) the order
 * the kind declares without the strategies R9 added — what the trace compares with. */
function kindCandidates(page: Page, base: Base, plan: KindPlan, legacy = false): Candidate[] {
  const out: Candidate[] = [];
  for (const {kind, qualifier} of [plan, ...plan.alternatives]) {
    if (!qualifier) {
      out.push({loc: base.locator(kind.selector), kind, qualifier, strategy: 'kind'});
      continue;
    }
    for (const strategy of legacy ? kind.match.filter((s) => s !== 'accessible' && s !== 'item') : orderOf(kind)) {
      for (const loc of byStrategy(page, base, kind, strategy, qualifier))
        out.push({loc, kind, qualifier, strategy});
    }
  }
  return out;
}

async function firstHit(all: Candidate[]): Promise<Candidate | undefined> {
  for (const c of all) {
    if (await c.loc.count() > 0)
      return c;
  }
  return undefined;
}

/** `BDD_RESOLVE_TRACE=<file>`: one JSON line for every resolution the R9 order decides by another
 * strategy than the kind's declared order did, saying whether both land on the same elements — the
 * proof that no phrase of a suite changed its target. */
async function trace(page: Page, base: Base, ref: NounRef, plan: KindPlan, chosen: Candidate | undefined): Promise<void> {
  const legacy = await firstHit(kindCandidates(page, base, plan, true));
  if (!chosen || !legacy || (legacy.strategy === chosen.strategy && legacy.kind === chosen.kind))
    return;
  const ids = (loc: Locator): Promise<number[]> => loc.evaluateAll((els) => els.map((e) =>
    (e as any).__bddRid ??= (window as any).__bddRidSeq = ((window as any).__bddRidSeq ?? 0) + 1));
  const before = await ids(legacy.loc);
  const after = await ids(chosen.loc);
  appendFileSync(process.env.BDD_RESOLVE_TRACE!, JSON.stringify({test: testTitle(), phrase: ref.raw, kind: chosen.kind.name,
    qualifier: chosen.qualifier,
    legacy: {strategy: legacy.strategy, count: before.length, ...(legacy.kind === chosen.kind ? {} : {kind: legacy.kind.name})},
    chosen: {strategy: chosen.strategy, count: after.length},
    same: before.length === after.length && before.every((id) => after.includes(id))}) + '\n');
}

function testTitle(): string {
  try {
    return test.info().titlePath.join(' › ');
  }
  catch {
    return '';
  }
}

const VISIBLE: MatchStrategy[] = ['label', 'accessible', 'title', 'text', 'exact-text', 'aria', 'placeholder', 'item'];

/** R9: what the user sees first, data-u2-name as the fallback; dart keeps its place. `accessible`
 * joins right after the first of `label`/`aria`. */
export function effectiveOrder(match: MatchStrategy[]): MatchStrategy[] {
  const out: MatchStrategy[] = match.filter((s) => s !== 'name');
  if (!out.includes('accessible')) {
    const i = out.findIndex((s) => s === 'label' || s === 'aria');
    if (i >= 0)
      out.splice(i + 1, 0, 'accessible');
  }
  if (!match.includes('name'))
    return out;
  let last = -1;
  for (const [i, s] of out.entries()) {
    if (VISIBLE.includes(s))
      last = i;
  }
  out.splice(last >= 0 ? last + 1 : match.indexOf('name'), 0, 'name');
  return out;
}

const orders = new WeakMap<KindEntry, MatchStrategy[]>();

function orderOf(kind: KindEntry): MatchStrategy[] {
  if (!orders.has(kind))
    orders.set(kind, effectiveOrder(kind.match));
  return orders.get(kind)!;
}

async function scopeBase(page: Page, scopeRef: NounRef, within?: Locator): Promise<Base> {
  if (scopeRef.scope)
    return locateRef(page, scopeRef.scope, within);
  if (!within && contextFirst(scopeRef)) {
    const root = page.locator(scopeRef.context!.selector);
    if (await root.count() > 0)
      return root;
  }
  return within ?? page;
}

async function isKindLabel(scope: Locator, scopeRef: NounRef): Promise<boolean> {
  const label = scopeRef.plan.type === 'kind' ? scopeRef.plan.kind.labelSelector : undefined;
  if (!label || await scope.count() === 0)
    return false;
  const selector = label.split(',').map((s) => s.trim().replace(/^:scope\s*>\s*/, '')).join(', ');
  return scope.first().evaluate((e, s) => e.matches(s), selector).catch(() => false);
}

/** The element's whole text is the qualifier, allowing decoration around it (an icon glyph, a
 * trailing colon). */
export function exactText(q: string): RegExp {
  return new RegExp(`^\\W*${escapeRegExp(q)}\\W*$`, 'i');
}

function byStrategy(page: Page, base: Base, kind: KindEntry, strategy: MatchStrategy, q: string): Locator[] {
  const compact = q.replace(/\s+/g, '');
  // Dart names join a menu path with "---": "Markers > Size" is div-Markers---Size
  const dashed = q.replace(/\s*>\s*/g, '---').replace(/\s+/g, '-');
  switch (strategy) {
    case 'name':
      return [compact, ...(dashed !== compact ? [dashed] : [])]
        .map((n) => base.locator(withAttr(kind.selector, `[data-u2-name="${cssString(n)}" i]`)));
    case 'label':
    case 'title':
      return kind.labelSelector ? [byLabel(page, base, kind, q)] : [];
    case 'accessible':
      // the label the browser associates (`<label for>`, aria-labelledby, aria-label): with the
      // element inside the kind's root (an input's editor), or on the root itself
      return [base.locator(kind.selector).filter({has: page.getByLabel(exactText(q))}),
        base.getByLabel(exactText(q)).and(page.locator(kind.selector))];
    case 'text':
      return [base.locator(kind.selector).filter({hasText: exactText(q)})];
    case 'exact-text':
      return [base.getByText(exactText(q))];
    case 'aria':
      return [base.locator(withAttr(kind.selector, `[aria-label="${cssString(q)}" i]`)),
        base.locator(withAttr(kind.selector, `[title="${cssString(q)}" i]`))];
    case 'placeholder':
      return [base.locator(kind.selector).filter({has: page.locator(`[placeholder="${cssString(q)}" i]`)})];
    case 'item':
      return [base.locator(withAttr(kind.selector, `[data-u2-item="${cssString(q)}" i]`))];
    case 'dart': {
      // Dart's `annotate` turns these into dashes as well (html_utils.dart): molecule_dictionary
      // is div-table-molecule-dictionary
      const annotated = dashed.replace(/[:_;*\\[\]{}|]/g, '-');
      return (kind.dartNames ?? []).flatMap((template) => [...new Set([annotated, dashed, q])]
        .map((name) => base.locator(withAttr(kind.selector, `[name="${cssString(template.replace('{q}', name))}" i]`))));
    }
  }
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
