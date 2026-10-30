import type {Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {NOT_SHOWN, until} from './layout.js';
import {describeCandidates, locateActionable} from './locate.js';

/** The rules of the engine's probe (`OperabilityRule` in u2core inspect.ts). */
export type OperabilityRule = 'covered' | 'no-pointer' | 'inert' | 'disabled-without-reason' | 'name-mismatch' |
  'oversized-option';

/** A finding of `rule` is dropped when its element is inside `selector` (`closest`). */
export interface OperabilityAllow {rule: OperabilityRule; selector: string; reason: string}

/** `detail` is the probe's sentence and, outside a claim, where it was found (`… — in dialog "Submit"`). */
export interface OperabilityFinding {rule: OperabilityRule; detail: string}

/** `probed` — the roots checked; `controls` — whether they hold anything to operate; `missing` — why nothing was
 * probed (a js-api without the probe). */
export interface OperabilityRead {findings: OperabilityFinding[]; probed: number; controls: boolean; missing?: string}

/** `new` — the u2 overlays no probe of this document has seen (after a step); `all` — the current view's u2 controls
 * and every u2 overlay (a scenario's end); `scoped` — the element itself. */
export type OperabilityScope = 'new' | 'all' | 'scoped';

interface ProbeArg {scope: OperabilityScope; allow: OperabilityAllow[]}

/** In the page; sees only its own source. The engine's probe (`DG.U2.operability`) with hit-testing on: in a browser
 * `elementFromPoint` answers, so `covered` is real. In the view, only what a u2 control owns (`DG.U2.control`) is
 * judged. Overlays are numbered in the order a probe first saw them: a tooltip, a balloon or a layer opened later
 * lying over a control is the user's doing, an older layer over a newer one is a finding. */
const probe = async (host: Element, a: ProbeArg): Promise<OperabilityRead> => {
  const w = window as any;
  const U2 = w.DG?.U2;
  if (typeof U2?.operability !== 'function')
    return {findings: [], probed: 0, controls: false, missing: 'the platform has no DG.U2.operability (a js-api before the u2 engine 1.2)'};
  const OVERLAYS = '[data-u2="dialog"], [data-u2-part="popup"], [data-u2="menu"]';
  const LAYERS = `${OVERLAYS}, .d4-dialog, .d4-menu-popup, .u2-dialog-backdrop, .d4-modal-background`;
  const TRANSIENT = '[role="tooltip"], .d4-tooltip, .d4-balloon, .d4-balloon-container, [data-u2="notify"], .ui-hint-popup';
  const CONTROLS = 'button, a[href], input, select, textarea, summary, [contenteditable], [tabindex], [role]';
  const state: {seen: WeakMap<Element, number>; next: number} = w.__bddOperability ??= {seen: new WeakMap(), next: 0};
  const shown = (e: Element) => e.getClientRects().length > 0 && e.closest('[hidden]') === null &&
    getComputedStyle(e).visibility !== 'hidden';
  let current: unknown = null;
  try {
    current = w.grok?.shell?.v?.root;
  }
  catch {
    current = null;
  }
  const view = current instanceof Element && current.isConnected ? current : null;
  const page = view ?? document.body;
  const overlays = [...document.querySelectorAll(OVERLAYS)].filter((o) => shown(o) && !view?.contains(o));
  const tops = overlays.filter((o) => !overlays.some((p) => p !== o && p.contains(o)));
  const fresh = tops.filter((o) => !state.seen.has(o));
  for (const o of fresh)
    state.seen.set(o, state.next++);
  const roots = a.scope === 'scoped' ? [host] : a.scope === 'new' ? fresh : [page, ...tops];
  if (roots.length === 0)
    return {findings: [], probed: 0, controls: false};
  const frame = () => new Promise<void>((resolve) => {
    requestAnimationFrame(() => resolve());
    setTimeout(resolve, 100);
  });
  await frame();
  await frame();
  const busyIn = () => roots.flatMap((r) => [...r.matches('[aria-busy="true"]') ? [r] : [], ...r.querySelectorAll('[aria-busy="true"]')])
    .filter((e) => shown(e) && U2.control(e)?.root === e);
  if (busyIn().length > 0) {
    await new Promise<void>((resolve) => {
      const done = () => {
        observer.disconnect();
        clearTimeout(timer);
        resolve();
      };
      const observer = new MutationObserver(() => {
        if (busyIn().length === 0)
          done();
      });
      for (const r of roots)
        observer.observe(r, {subtree: true, childList: true, attributes: true});
      const timer = setTimeout(done, 2000);
    });
  }
  const busy = busyIn();

  const layerOf = (e: Element) => tops.find((o) => o.contains(e));
  const named = (e: Element): string => {
    let label = '';
    try {
      label = U2.control(e)?.root === e ? U2.state(e)?.label ?? '' : '';
    }
    catch {
      label = '';
    }
    return label || e.getAttribute('aria-label') || e.getAttribute('data-u2-name') || '';
  };
  const where = (root: Element): string => {
    const layer = layerOf(root);
    if (layer === undefined) {
      let name = '';
      try {
        name = String(w.grok?.shell?.v?.name ?? '');
      }
      catch {
        name = '';
      }
      return view ? `view "${name}"` : 'the page';
    }
    const label = named(layer);
    const owner = layer.getAttribute('data-u2-owner');
    return `${layer.getAttribute('data-u2') ?? 'popup'}${label ? ` "${label}"` : owner ? ` of "${owner}"` : ''}`;
  };
  const covers = (el: Element, root: Element): boolean => {
    const r = el.getBoundingClientRect();
    const hit = document.elementFromPoint(r.left + r.width / 2, r.top + r.height / 2);
    if (hit === null || hit.closest(TRANSIENT) !== null)
      return false;
    let order: number | undefined;
    for (let at: Element | null = hit; at !== null && order === undefined; at = at.parentElement)
      order = state.seen.get(at);
    const own = layerOf(root);
    if (own === undefined)
      return order === undefined && hit.closest(LAYERS) === null && (view === null || view.contains(hit));
    if (root.contains(hit))
      return true;
    return order !== undefined ? order < state.seen.get(own)! : hit.closest(LAYERS) === null;
  };
  const judged = (el: Element, root: Element) => root !== page ||
    (U2.control(el) !== undefined && el.closest(`${LAYERS}, ${TRANSIENT}`) === null);

  const findings: OperabilityFinding[] = [];
  const reported = new Set<string>();
  let controls = false;
  for (const root of roots) {
    if (!root.isConnected)
      continue;
    controls ||= root.matches(CONTROLS) || root.querySelector(CONTROLS) !== null;
    for (const f of U2.operability(root, {hitTest: true}) as {rule: OperabilityRule; element: Element; detail: string}[]) {
      if (!judged(f.element, root) || busy.some((b) => b.contains(f.element)) ||
        (f.rule === 'covered' && !covers(f.element, root)) ||
        a.allow.some((e) => e.rule === f.rule && f.element.closest(e.selector) !== null))
        continue;
      const detail = a.scope === 'scoped' ? f.detail : `${f.detail} — in ${where(root)}`;
      if (!reported.has(`${f.rule}: ${detail}`)) {
        reported.add(`${f.rule}: ${detail}`);
        findings.push({rule: f.rule, detail});
      }
    }
  }
  return {findings, probed: roots.length, controls};
};

/** The probe over the page in one roundtrip: the waits (two frames, busy controls up to 2 s) run in the page. */
export function readOperability(page: Page, scope: Exclude<OperabilityScope, 'scoped'>,
  allow: OperabilityAllow[] = []): Promise<OperabilityRead> {
  return page.locator(':root').evaluate(probe, {scope, allow}, {timeout: 10000});
}

/** Every control in the element operable as offered, polled; an element holding no control fails, since nothing
 * would be checked, and so does a platform without the probe. */
export async function expectOperable(page: Page, target: ElementRef): Promise<void> {
  const loc = await locateActionable(page, target);
  await until(async () => {
    const n = await loc.count();
    if (n !== 1)
      return n === 0 ? target.phrase + NOT_SHOWN : `${target.phrase} names ${n} shown elements: ${await describeCandidates(loc)}`;
    const read = await loc.evaluate(probe, {scope: 'scoped' as const, allow: []});
    if (read.missing !== undefined)
      return read.missing;
    if (!read.controls)
      return `${target.phrase} holds no control to operate`;
    return read.findings.length === 0 ? 'none' : `${read.findings.length} control(s) in ${target.phrase} cannot be operated ` +
      `as offered: ${read.findings.map((f) => `${f.rule}: ${f.detail}`).join('; ')}`;
  }, 'none', `the controls in ${target.phrase}`);
}
