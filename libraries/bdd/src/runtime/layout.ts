import type {Locator, Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {describeCandidates, locateActionable} from './locate.js';
import {expect} from './patience.js';
import type {Box} from './viewer-runtime.js';

export type LayoutRule = 'zero-size' | 'offscreen' | 'overlap' | 'clipped-text' | 'truncated-text' |
  'overflow' | 'label-misaligned' | 'label-column';

/** A finding of `rule` is dropped when one of its elements is inside `selector` (`closest`). */
export interface AllowEntry {rule: LayoutRule; selector: string; reason: string}

/** `what` names the elements as the user sees them (`"Dose, mg" number input`, `label of "Dose, mg"
 * number input`); `detail` is one sentence with the measurements. */
export interface LayoutFinding {rule: LayoutRule; what: string; detail: string}

/** `u2: false` — the scope shows no u2 control and nothing was audited. `ms` is the in-page
 * measurement, waits excluded; `capped` — more than 400 controls, the first 400 audited; `busy` —
 * controls still busy after 2 s, left out with everything inside them. */
export interface LayoutResult {u2: boolean; controls: number; findings: LayoutFinding[]; ms: number; capped: boolean; busy: number}

const CAP = 400;

interface AuditArg {scoped: boolean; allow: AllowEntry[]; cap: number}

/** In the page; sees only its own source. `host` is the scope when `scoped`, else the current view's
 * root and every visible u2 overlay make it, or the whole body when there is no view. */
const audit = async (host: Element, a: AuditArg): Promise<LayoutResult> => {
  const frame = () => new Promise<void>((resolve) => {
    requestAnimationFrame(() => resolve());
    setTimeout(resolve, 100);
  });
  await frame();
  await frame();
  const shown = (e: Element) => (typeof e.checkVisibility === 'function' ? e.checkVisibility({checkVisibilityCSS: true}) :
    e.getClientRects().length > 0) && e.closest('[hidden]') === null;
  const scope = (): Element[] => {
    if (a.scoped)
      return [host];
    let root: unknown = null;
    try {
      root = (window as any).grok?.shell?.v?.root;
    }
    catch {
      root = null;
    }
    if (!(root instanceof Element) || !root.isConnected)
      return [document.body];
    const view: Element = root;
    const overlays = [...document.querySelectorAll('[data-u2="dialog"], [data-u2-part="popup"], [data-u2="menu"]')]
      .filter((o) => shown(o) && !view.contains(o));
    return [view, ...overlays.filter((o) => !overlays.some((p) => p !== o && p.contains(o)))];
  };
  const within = (roots: Element[], selector: string): Element[] => {
    const found = new Set<Element>();
    for (const r of roots) {
      if (r.matches(selector))
        found.add(r);
      for (const e of r.querySelectorAll(selector))
        found.add(e);
    }
    return [...found].filter(shown);
  };
  if (within(scope(), '[data-u2]').length === 0)
    return {u2: false, controls: 0, findings: [], ms: 0, capped: false, busy: 0};
  const deadline = performance.now() + 2000;
  while (within(scope(), '[data-u2][aria-busy="true"]').length > 0 && performance.now() < deadline)
    await new Promise((resolve) => setTimeout(resolve, 100));

  const t0 = performance.now();
  const roots = scope();
  const busy = within(roots, '[data-u2][aria-busy="true"]');
  const all = within(roots, '[data-u2]').filter((e) => !busy.some((b) => b.contains(e)));
  const controls = all.slice(0, a.cap);
  const isControl = new Set(controls);
  const styles = new Map<Element, CSSStyleDeclaration>();
  const css = (e: Element): CSSStyleDeclaration => {
    let s = styles.get(e);
    if (!s)
      styles.set(e, s = getComputedStyle(e));
    return s;
  };
  const rects = new Map<Element, DOMRect>();
  const rect = (e: Element): DOMRect => {
    let r = rects.get(e);
    if (!r)
      rects.set(e, r = e.getBoundingClientRect());
    return r;
  };
  const px = (n: number) => Math.round(n);
  const boxless = (e: Element) => /^(inline|contents)$/.test(css(e).display);
  const escapes = (e: Element) => /^(absolute|fixed)$/.test(css(e).position);
  const layered = (e: Element) => /^(absolute|fixed|sticky)$/.test(css(e).position);

  const U2 = (window as any).DG?.U2;
  const ownPart = (e: Element, selector: string) => [...e.querySelectorAll(selector)]
    .find((p) => p.parentElement?.closest('[data-u2]') === e);
  const firstText = (e: Element) => ((e as HTMLElement).innerText ?? e.textContent ?? '').trim().replace(/\s+/g, ' ').slice(0, 30);
  const nameOf = (e: Element): string => {
    let label = '';
    try {
      label = U2?.control?.(e)?.root === e ? U2.state(e)?.label ?? '' : '';
    }
    catch {
      label = '';
    }
    label = label || ownPart(e, '[data-u2-part="label"], [data-u2-part="title"]')?.textContent?.trim() ||
      e.getAttribute('aria-label') || e.getAttribute('data-u2-name') || firstText(e);
    const noun = (e.getAttribute('data-u2') || e.getAttribute('role') || e.tagName.toLowerCase()).replace(/-/g, ' ');
    return label ? `"${label}" ${noun}` : noun;
  };
  const textName = (t: Element): string => {
    if (t.hasAttribute('data-u2'))
      return nameOf(t);
    const owner = t.parentElement?.closest('[data-u2]');
    const of = owner ? nameOf(owner) : 'the page';
    const part = t.getAttribute('data-u2-part');
    if (part === 'item-label')
      return `"${firstText(t)}" item of ${of}`;
    return part ? `${part.replace(/-/g, ' ')} of ${of}` : `"${firstText(t)}" ${t.getAttribute('role') || t.tagName.toLowerCase()} of ${of}`;
  };

  type Edges = {l: number; t: number; r: number; b: number};
  // pairs compare within one layer (the nearest positioned ancestor), so a dialog never overlaps the page under it
  const layers = new Map<Element, Element | null>();
  const layerOf = (e: Element | null): Element | null => {
    if (!e)
      return null;
    if (!layers.has(e))
      layers.set(e, layered(e) ? e : layerOf(e.parentElement));
    return layers.get(e)!;
  };
  // where the content of `e` can be seen: its own padding box where it clips, cut by its ancestors' up
  // to the nearest absolutely or fixed positioned one
  const clips = new Map<Element, Edges>();
  const clipOf = (e: Element | null): Edges => {
    if (!e || e === document.documentElement)
      return {l: -Infinity, t: -Infinity, r: Infinity, b: Infinity};
    let c = clips.get(e);
    if (c)
      return c;
    c = escapes(e) ? {l: -Infinity, t: -Infinity, r: Infinity, b: Infinity} : {...clipOf(e.parentElement)};
    const s = css(e);
    if (!boxless(e) && (s.overflowX !== 'visible' || s.overflowY !== 'visible')) {
      const r = rect(e);
      const he = e as HTMLElement;
      const l = r.left + he.clientLeft;
      const t = r.top + he.clientTop;
      c = {l: Math.max(c.l, l), t: Math.max(c.t, t), r: Math.min(c.r, l + he.clientWidth), b: Math.min(c.b, t + he.clientHeight)};
    }
    clips.set(e, c);
    return c;
  };
  const scrolls = (e: Element): boolean => {
    for (let p = e.parentElement; p; p = p.parentElement) {
      const s = css(p);
      if (/auto|scroll/.test(s.overflowX + s.overflowY) && (p.scrollWidth > p.clientWidth + 1 || p.scrollHeight > p.clientHeight + 1))
        return true;
    }
    const doc = document.scrollingElement;
    return doc !== null && (doc.scrollWidth > doc.clientWidth + 1 || doc.scrollHeight > doc.clientHeight + 1);
  };

  const findings: (LayoutFinding & {els: Element[]})[] = [];
  const add = (rule: LayoutRule, els: Element[], what: string, detail: string) => {
    if (!a.allow.some((x) => x.rule === rule && els.some((e) => e.closest(x.selector) !== null)))
      findings.push({rule, what, detail, els});
  };

  // what a box shows: text, or a shown element that holds it, or a shown leaf other than an empty div or span
  // (an editor, an icon, an image)
  const holds = (e: Element): boolean => [...e.childNodes].some((n) => n.nodeType === Node.TEXT_NODE ? (n.textContent ?? '').trim() !== '' :
    n instanceof Element && shown(n) && (n.childElementCount > 0 ? holds(n) : (n.textContent ?? '').trim() !== '' || !/^(DIV|SPAN)$/.test(n.tagName)));
  const sized: Element[] = [];
  for (const e of controls) {
    const r = rect(e);
    if (r.width < 1 || r.height < 1) {
      // a container that shows nothing takes no room by right: an empty one, or a form section whose fields are all hidden
      if (holds(e))
        add('zero-size', [e], nameOf(e), `${nameOf(e)} has no size on the page (${px(r.width)}×${px(r.height)} px)`);
      continue;
    }
    sized.push(e);
    if ((r.right <= 0 || r.bottom <= 0 || r.left >= innerWidth || r.top >= innerHeight) && !scrolls(e))
      add('offscreen', [e], nameOf(e), `${nameOf(e)} lies outside the window at (${px(r.left)}, ${px(r.top)}) and nothing scrolls to it`);
  }

  const parties = sized.filter((e) => !layered(e)).map((e) => {
    const r = rect(e);
    const c = clipOf(e.parentElement);
    return {e, layer: layerOf(e), l: Math.max(r.left, c.l), t: Math.max(r.top, c.t), r: Math.min(r.right, c.r), b: Math.min(r.bottom, c.b)};
  }).filter((p) => p.r > p.l && p.b > p.t);
  const overlapping: Element[][] = [];
  for (let i = 0; i < parties.length; i++) {
    for (let j = i + 1; j < parties.length; j++) {
      const p = parties[i];
      const q = parties[j];
      if (p.layer !== q.layer || p.e.contains(q.e) || q.e.contains(p.e))
        continue;
      const w = Math.min(p.r, q.r) - Math.max(p.l, q.l);
      const h = Math.min(p.b, q.b) - Math.max(p.t, q.t);
      // the controls inside two overlapping ones overlap too: the outermost pair is the finding
      if (w <= 2 || h <= 2 || overlapping.some(([x, y]) => (x.contains(p.e) && y.contains(q.e)) || (x.contains(q.e) && y.contains(p.e))))
        continue;
      overlapping.push([p.e, q.e]);
      add('overlap', [p.e, q.e], `${nameOf(p.e)} and ${nameOf(q.e)}`,
        `${nameOf(p.e)} overlaps ${nameOf(q.e)} by ${px(w)}×${px(h)} px at (${px(Math.max(p.l, q.l))}, ${px(Math.max(p.t, q.t))})`);
    }
  }

  let ctx: CanvasRenderingContext2D | null | undefined;
  const fitting = (n: HTMLElement, ellipsis: boolean): string => {
    const text = (n.textContent ?? '').trim().replace(/\s+/g, ' ');
    if (ctx === undefined)
      ctx = document.createElement('canvas').getContext('2d');
    if (!ctx)
      return text;
    const s = css(n);
    ctx.font = `${s.fontStyle} ${s.fontWeight} ${s.fontSize} ${s.fontFamily}`;
    const room = n.clientWidth - parseFloat(s.paddingLeft) - parseFloat(s.paddingRight);
    const tail = ellipsis ? '…' : '';
    let k = text.length;
    while (k > 0 && ctx.measureText(text.slice(0, k) + tail).width > room)
      k--;
    return text.slice(0, k) + tail;
  };
  const texts = within(roots, '[data-u2-part="label"], [data-u2-part="title"], [data-u2-part="item-label"], button, [role="tab"], [role="button"]')
    .filter((t) => {
      const owner = t.closest('[data-u2]');
      const r = rect(t);
      return owner !== null && isControl.has(owner) && r.width >= 1 && r.height >= 1 && (t.textContent ?? '').trim() !== '';
    });
  for (const t of texts) {
    const own = [t, ...[...t.querySelectorAll('*')].slice(0, 20)].find((n) => !boxless(n) && /hidden|clip/.test(css(n).overflowX) &&
      n.scrollWidth > n.clientWidth + 1 && (n.textContent ?? '').trim() !== '') as HTMLElement | undefined;
    if (own) {
      const ellipsis = css(own).textOverflow === 'ellipsis';
      // a list ellipsizes its item labels by design
      if (!(ellipsis && t.getAttribute('data-u2-part') === 'item-label'))
        add(ellipsis ? 'truncated-text' : 'clipped-text', [t], textName(t),
          `${textName(t)} is cut to "${fitting(own, ellipsis)}" (needs ${px(own.scrollWidth)} px, has ${px(own.clientWidth)})`);
      continue;
    }
    // cut by the first ancestor that clips it; one that scrolls is the user's to scroll, and a text
    // wholly outside its clipping ancestor is put away (a collapsed pane), not cut. A scroller that
    // runs past the ancestor cutting the text (`height: 100%` plus padding) hides a strip of its own
    // view for good: that is the defect the finding names.
    const r = rect(t);
    let scroller: Element | null = null;
    for (let p = t.parentElement; p && p !== document.documentElement; p = p.parentElement) {
      const s = css(p);
      if (!boxless(p) && (s.overflowX !== 'visible' || s.overflowY !== 'visible')) {
        const pr = rect(p);
        const l = pr.left + p.clientLeft;
        const top = pr.top + p.clientTop;
        const cutX = r.left < l - 1 || r.right > l + p.clientWidth + 1;
        const cutY = r.top < top - 1 || r.bottom > top + p.clientHeight + 1;
        if (cutX || cutY) {
          const w = Math.min(r.right, l + p.clientWidth) - Math.max(r.left, l);
          const h = Math.min(r.bottom, top + p.clientHeight) - Math.max(r.top, top);
          if (!/auto|scroll/.test(cutX ? s.overflowX : s.overflowY) && w >= 1 && h >= 1) {
            const by = p.closest('[data-u2]');
            const sr = scroller ? rect(scroller) : null;
            const over = !sr ? 0 : cutX ? Math.max(sr.right - l - p.clientWidth, l - sr.left) : Math.max(sr.bottom - top - p.clientHeight, top - sr.top);
            add('clipped-text', [t], textName(t), `${textName(t)} is cut off by ${by ? nameOf(by) : 'its container'}: ` +
              `${px(w)}×${px(h)} of ${px(r.width)}×${px(r.height)} px show` +
              (over > 1 ? `; the scroller inside it runs ${px(over)} px past its edge` : ''));
          }
          break;
        }
        if (/auto|scroll/.test(s.overflowX + s.overflowY))
          scroller = p;
      }
      if (escapes(p))
        break;
    }
  }

  const spilled: {e: Element; sw: number; sh: number; cw: number; ch: number}[] = [];
  const range = document.createRange();
  for (const e of sized) {
    const s = css(e);
    const he = e as HTMLElement;
    if (boxless(e) || !((s.overflowX === 'visible' && he.scrollWidth > he.clientWidth + 1) ||
      (s.overflowY === 'visible' && he.scrollHeight > he.clientHeight + 1)))
      continue;
    // the scroll size counts positioned decorations too: confirm with the content in the flow
    const r = rect(e);
    const right = r.left + he.clientLeft + he.clientWidth + 1;
    const bottom = r.top + he.clientTop + he.clientHeight + 1;
    const stack: Node[] = [...e.childNodes];
    let spills = false;
    for (let seen = 0; stack.length > 0 && seen < 300 && !spills; seen++) {
      const n = stack.pop()!;
      let nr: DOMRect;
      if (n.nodeType === Node.TEXT_NODE) {
        if ((n.textContent ?? '').trim() === '')
          continue;
        range.selectNodeContents(n);
        nr = range.getBoundingClientRect();
      }
      else if (n.nodeType === Node.ELEMENT_NODE) {
        const c = n as Element;
        const cs = css(c);
        if (escapes(c) || cs.display === 'none')
          continue;
        nr = rect(c);
        if (cs.display === 'contents' || (cs.overflowX === 'visible' && cs.overflowY === 'visible'))
          stack.push(...c.childNodes);
      }
      else
        continue;
      spills = nr.width > 0 && nr.height > 0 && (nr.right > right || nr.bottom > bottom);
    }
    if (spills)
      spilled.push({e, sw: he.scrollWidth, sh: he.scrollHeight, cw: he.clientWidth, ch: he.clientHeight});
  }
  // a control whose child control spills spills too: the innermost is the finding
  for (const o of spilled) {
    if (!spilled.some((d) => d !== o && o.e.contains(d.e)))
      add('overflow', [o.e], nameOf(o.e), `content spills out of ${nameOf(o.e)}: it needs ${px(o.sw)}×${px(o.sh)} px, has ${px(o.cw)}×${px(o.ch)}`);
  }

  const labelLeft: {e: Element; editor: DOMRect}[] = [];
  for (const e of sized) {
    const label = [...e.children].find((c) => c.getAttribute('data-u2-part') === 'label');
    const editor = ownPart(e, '[data-u2-part="editor"]');
    if (!label || !editor || !shown(label) || !shown(editor))
      continue;
    const lr = rect(label);
    const er = rect(editor);
    if (lr.width < 1 || lr.height < 1 || er.width < 1 || er.height < 1)
      continue;
    const beside = lr.top < er.bottom && er.top < lr.bottom && (lr.right <= er.left + 2 || lr.left >= er.right - 2);
    if (beside) {
      if (lr.right <= er.left + 2)
        labelLeft.push({e, editor: er});
      // beside a multi-line editor the label belongs to the first line; beside one taller than a row
      // (a two-row check list) to its first line or its middle, whichever it is nearer
      const tall = er.height > 2 * lr.height;
      const middle = (lr.top + lr.height / 2) - (er.top + er.height / 2);
      const first = lr.top - er.top;
      const onFirst = tall || (er.height > lr.height + 4 && Math.abs(first) < Math.abs(middle));
      const d = onFirst ? first : middle;
      if (Math.abs(d) > (tall ? lr.height / 2 + 2 : 2))
        add('label-misaligned', [e], `label of ${nameOf(e)}`, `label of ${nameOf(e)} sits ${px(Math.abs(d))} px ${d < 0 ? 'above' : 'below'} ` +
          `the ${onFirst ? 'first line' : 'middle'} of its editor`);
    }
    else if (lr.bottom <= er.top + 2 && Math.abs(lr.left - er.left) > 2)
      add('label-misaligned', [e], `label of ${nameOf(e)}`,
        `label of ${nameOf(e)} starts ${px(Math.abs(lr.left - er.left))} px ${lr.left > er.left ? 'right' : 'left'} of its editor`);
  }
  const forms = new Map<Element, {e: Element; editor: DOMRect}[]>();
  for (const i of labelLeft) {
    const form = i.e.parentElement?.closest('[data-u2="form"], [data-u2="object-form"], [data-u2="func-form"]');
    if (!form || !isControl.has(form))
      continue;
    // an input inside another control's editor (a columns map's key rows) is that control's content,
    // not a row of the form
    const host = i.e.parentElement?.closest('[data-u2-part="editor"]');
    if (host && form.contains(host))
      continue;
    const inputs = forms.get(form) ?? [];
    inputs.push(i);
    forms.set(form, inputs);
  }
  // a cell of a form's columns (Form.addColumns) lines up with its own column: measured from the grid's
  // first column, its editor joins the form's one column like any other row's
  const leftOf = (i: {e: Element; editor: DOMRect}) => {
    const cell = i.e.parentElement?.matches('[data-u2-part="columns"]') ? i.e : null;
    return cell === null ? i.editor.left :
      i.editor.left - (cell.getBoundingClientRect().left - cell.parentElement!.getBoundingClientRect().left);
  };
  for (const [form, inputs] of forms) {
    const lefts = inputs.map(leftOf);
    let column = lefts[0];
    let most = 0;
    for (const left of lefts) {
      const n = lefts.filter((l) => Math.abs(l - left) <= 2).length;
      if (n > most) {
        most = n;
        column = left;
      }
    }
    const out = inputs.map((i, k) => ({...i, left: lefts[k]})).filter((i) => Math.abs(i.left - column) > 2);
    if (out.length > 0)
      add('label-column', [form, ...out.map((o) => o.e)], nameOf(form), `the editors of ${nameOf(form)} do not start in one column: ` +
        `${out.map((o) => `${nameOf(o.e)} at ${px(o.left)} px`).join(', ')}; the other ${inputs.length - out.length} at ${px(column)} px`);
  }

  return {u2: true, controls: controls.length, findings: findings.map(({rule, what, detail}) => ({rule, what, detail})),
    ms: Math.round((performance.now() - t0) * 10) / 10, capped: all.length > a.cap, busy: busy.length};
};

/** The layout audit of `opts.root` (one element), else of the current view and the u2 overlays over it,
 * else of the page. One roundtrip: the waits (two frames, busy controls up to 2 s) run in the page. */
export async function auditLayout(page: Page, opts: {root?: Locator; allow?: AllowEntry[]} = {}): Promise<LayoutResult> {
  return (opts.root ?? page.locator('body')).evaluate(audit, {scoped: opts.root !== undefined, allow: opts.allow ?? [], cap: CAP});
}

export const NOT_SHOWN = ' is not shown';

/** Polls a reading until it is `ok`; a timeout fails with the last reading as the sentence, except an
 * element not on the page, which keeps Playwright's wait failure so the report lists what is shown. */
export async function until(read: () => Promise<string>, ok: string, message: string): Promise<void> {
  let last = '';
  await expect.poll(async () => last = await read(), {message}).toBe(ok).catch((e: unknown) => {
    throw last.endsWith(NOT_SHOWN) ? e : new Error(last);
  });
}

/** The box of the one visible match, or why there is none to read. */
async function readBox(loc: Locator, phrase: string): Promise<Box | string> {
  const boxes = await loc.evaluateAll((els) => els.map((e) => {
    const r = e.getBoundingClientRect();
    return {x: r.x, y: r.y, width: r.width, height: r.height};
  }));
  if (boxes.length === 1)
    return boxes[0];
  return boxes.length === 0 ? phrase + NOT_SHOWN : `${phrase} names ${boxes.length} shown elements: ${await describeCandidates(loc)}`;
}

/** The box of the element a phrase names (its one visible match) and that element as the user sees it. */
export async function rectOf(page: Page, target: ElementRef): Promise<Box & {what: string}> {
  const loc = await locateActionable(page, target);
  const box = await readBox(loc, target.phrase);
  if (typeof box === 'string')
    throw new Error(box);
  return {...box, what: await describeCandidates(loc)};
}

export type Placement = 'left of' | 'right of' | 'above' | 'below' | 'inside' | 'clear of' | 'left-aligned with' |
  'top-aligned with' | 'vertically centered with' | 'as tall as' | 'as wide as';

const num = (n: number) => String(Math.round(n * 10) / 10);
const span = (b: Box) => `${num(b.x)}..${num(b.x + b.width)} × ${num(b.y)}..${num(b.y + b.height)}`;
const side = (d: number, plus: string, minus: string) => `${num(Math.abs(d))} px ${d > 0 ? plus : minus}`;

/** Whether `a` stands so against `b` (1 px of rounding allowed), and if not, where it is instead. */
const PLACEMENTS: Record<Placement, (a: Box, b: Box) => string | null> = {
  'left of': (a, b) => a.x + a.width <= b.x + 1 ? null : `its right edge is ${num(a.x + a.width - b.x)} px past the left edge of the other`,
  'right of': (a, b) => a.x >= b.x + b.width - 1 ? null : `its left edge is ${num(b.x + b.width - a.x)} px short of the right edge of the other`,
  'above': (a, b) => a.y + a.height <= b.y + 1 ? null : `its bottom is ${num(a.y + a.height - b.y)} px below the top of the other`,
  'below': (a, b) => a.y >= b.y + b.height - 1 ? null : `its top is ${num(b.y + b.height - a.y)} px above the bottom of the other`,
  'inside': (a, b) => {
    const out = ([[b.x - a.x, 'on the left'], [a.x + a.width - b.x - b.width, 'on the right'], [b.y - a.y, 'at the top'],
      [a.y + a.height - b.y - b.height, 'at the bottom']] as [number, string][]).filter(([d]) => d > 1);
    return out.length === 0 ? null : `it sticks out ${out.map(([d, where]) => `${num(d)} px ${where}`).join(', ')}`;
  },
  'clear of': (a, b) => {
    const w = Math.min(a.x + a.width, b.x + b.width) - Math.max(a.x, b.x);
    const h = Math.min(a.y + a.height, b.y + b.height) - Math.max(a.y, b.y);
    return w <= 1 || h <= 1 ? null : `they overlap by ${num(w)}×${num(h)} px`;
  },
  'left-aligned with': (a, b) => Math.abs(a.x - b.x) <= 1 ? null : `its left edge is ${side(a.x - b.x, 'right', 'left')} of the other's`,
  'top-aligned with': (a, b) => Math.abs(a.y - b.y) <= 1 ? null : `its top is ${side(a.y - b.y, 'below', 'above')} the other's`,
  'vertically centered with': (a, b) => {
    const d = (a.y + a.height / 2) - (b.y + b.height / 2);
    return Math.abs(d) <= 1 ? null : `its middle is ${side(d, 'below', 'above')} the other's`;
  },
  'as tall as': (a, b) => Math.abs(a.height - b.height) <= 1 ? null : `it is ${num(a.height)} px tall, the other ${num(b.height)} px`,
  'as wide as': (a, b) => Math.abs(a.width - b.width) <= 1 ? null : `it is ${num(a.width)} px wide, the other ${num(b.width)} px`,
};

/** `a` stands as `placement` says against `b`, polled; the failure gives both boxes. */
export async function expectPlaced(page: Page, a: ElementRef, placement: Placement, b: ElementRef): Promise<void> {
  const la = await locateActionable(page, a);
  const lb = await locateActionable(page, b);
  await until(async () => {
    const ra = await readBox(la, a.phrase);
    if (typeof ra === 'string')
      return ra;
    const rb = await readBox(lb, b.phrase);
    if (typeof rb === 'string')
      return rb;
    const miss = PLACEMENTS[placement](ra, rb);
    return miss === null ? 'holds' : `${a.phrase} is not ${placement} ${b.phrase} — ${miss}: ${a.phrase} spans ${span(ra)}, ${b.phrase} ${span(rb)}`;
  }, 'holds', `${a.phrase} ${placement} ${b.phrase}`);
}

export async function expectHeight(page: Page, target: ElementRef, height: number): Promise<void> {
  const loc = await locateActionable(page, target);
  await until(async () => {
    const box = await readBox(loc, target.phrase);
    if (typeof box === 'string')
      return box;
    return Math.abs(box.height - height) <= 1 ? 'holds' : `${target.phrase} is ${num(box.height)} px tall, not ${height}: it spans ${span(box)}`;
  }, 'holds', `the height of ${target.phrase}`);
}

/** Every edge within the browser window: a fixed or portaled element (a dialog, a popup) that no
 * ancestor clips is checked against the window itself. */
export async function expectInWindow(page: Page, target: ElementRef): Promise<void> {
  const loc = await locateActionable(page, target);
  await until(async () => {
    const box = await readBox(loc, target.phrase);
    if (typeof box === 'string')
      return box;
    const view = await page.evaluate(() => ({x: 0, y: 0, width: innerWidth, height: innerHeight}));
    const miss = PLACEMENTS.inside(box, view);
    return miss === null ? 'holds' : `${target.phrase} is not inside the window — ${miss}: it spans ${span(box)}, ` +
      `the window ${num(view.width)}×${num(view.height)}`;
  }, 'holds', `${target.phrase} inside the window`);
}

/** The box width within 1 px of `width`, or with `bound` no narrower (`at least`) or no wider (`at most`) than it. */
export async function expectWidth(page: Page, target: ElementRef, width: number, bound?: 'at least' | 'at most'): Promise<void> {
  const loc = await locateActionable(page, target);
  await until(async () => {
    const box = await readBox(loc, target.phrase);
    if (typeof box === 'string')
      return box;
    const holds = bound === 'at least' ? box.width >= width - 1 : bound === 'at most' ? box.width <= width + 1 :
      Math.abs(box.width - width) <= 1;
    return holds ? 'holds' : `${target.phrase} is ${num(box.width)} px wide, not ${bound ? bound + ' ' : ''}${width}: it spans ${span(box)}`;
  }, 'holds', `the width of ${target.phrase}`);
}

/** The element whole inside every ancestor that clips it (an absolutely positioned one only by its
 * containing block and above, a fixed one by none), and its own text not overflowing its box. */
export async function expectNotClipped(page: Page, target: ElementRef): Promise<void> {
  const loc = await locateActionable(page, target);
  await until(async () => {
    const reads = await loc.evaluateAll((els) => els.map((e) => {
      const r = e.getBoundingClientRect();
      let escaping = getComputedStyle(e).position;
      for (let p = e.parentElement; p && p !== document.documentElement && escaping !== 'fixed'; p = p.parentElement) {
        const s = getComputedStyle(p);
        if (escaping === 'absolute' && s.position === 'static')
          continue;
        escaping = s.position === 'absolute' || s.position === 'fixed' ? s.position : '';
        if (/^(inline|contents)$/.test(s.display) || (s.overflowX === 'visible' && s.overflowY === 'visible'))
          continue;
        const pr = p.getBoundingClientRect();
        const l = pr.left + p.clientLeft;
        const t = pr.top + p.clientTop;
        const cut = ([[l - r.left, 'on the left'], [r.right - l - p.clientWidth, 'on the right'], [t - r.top, 'at the top'],
          [r.bottom - t - p.clientHeight, 'at the bottom']] as [number, string][]).filter(([d]) => d > 1)
          .map(([d, where]) => `${Math.round(d)} px ${where}`);
        if (cut.length > 0) {
          const by = p.closest('[data-u2]');
          const name = by?.getAttribute('data-u2-name') || by?.getAttribute('data-u2') || p.tagName.toLowerCase();
          return `it is cut off by ${by ? `the ${name.replace(/-/g, ' ')}` : `a ${name}`} it sits in: ${cut.join(', ')}`;
        }
      }
      const own = [e, ...[...e.querySelectorAll('*')].slice(0, 20)].find((n) => /hidden|clip/.test(getComputedStyle(n).overflowX) &&
        n.scrollWidth > n.clientWidth + 1 && (n.textContent ?? '').trim() !== '');
      return own ? `its text needs ${own.scrollWidth} px and has ${own.clientWidth}${getComputedStyle(own).textOverflow === 'ellipsis' ?
        ', ending in an ellipsis' : ''}` : '';
    }));
    if (reads.length !== 1)
      return reads.length === 0 ? target.phrase + NOT_SHOWN : `${target.phrase} names ${reads.length} shown elements: ${await describeCandidates(loc)}`;
    return reads[0] === '' ? 'whole' : `${target.phrase} is clipped: ${reads[0]}`;
  }, 'whole', `whether ${target.phrase} is clipped`);
}

/** The audit over the element's subtree finds nothing, polled: a busy control is waited for; an
 * element holding no u2 control fails, since nothing would be checked. */
export async function expectNoLayoutProblems(page: Page, target: ElementRef, allow: AllowEntry[] = []): Promise<void> {
  const loc = await locateActionable(page, target);
  await until(async () => {
    const n = await loc.count();
    if (n !== 1)
      return n === 0 ? target.phrase + NOT_SHOWN : `${target.phrase} names ${n} shown elements: ${await describeCandidates(loc)}`;
    const result = await auditLayout(page, {root: loc, allow});
    if (!result.u2)
      return `${target.phrase} holds no u2 control the layout audit reads`;
    if (result.busy > 0)
      return `${target.phrase} is still busy: ${result.busy} control(s) report aria-busy`;
    return result.findings.length === 0 ? 'none' : `${target.phrase} has ${result.findings.length} layout problem(s): ` +
      result.findings.map((f) => `${f.rule}: ${f.detail}`).join('; ');
  }, 'none', `the layout problems in ${target.phrase}`);
}
