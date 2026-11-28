/* A stage's main panels on a twelve-column grid, laid out from a list of (panel, width): the layout is data, so a
   user's or a lab's own is a different list, not different code. The page is customized in place (Customize this
   page, in the toolbox and File > View): each panel's header then offers a handle to drag it and × to hide it, and its
   right edge a grip to widen it (to a quarter, a third, a half, two thirds, the row; arrows too, a double-click
   the full row and back); Done keeps it per user (grok.userSettings), Save as the lab's for every user who has
   none, Reset goes back to the lab's or the stage's own. A panel the stage folds by
   default shows one line of what it holds; what a user unfolds or folds is remembered with their layout. Hidden
   panels stay in reach: the toolbox's index lists them, and a quiet line at the end of the page shows them again. As
   the view narrows, the panels widen to half and then to the full row (container queries in studio.css). */
import '../css/studio.css';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {div, span} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {block} from './ui';
import {_package} from './package';
import {canEditLabRules, saveLabSetting} from './outlier-rules';
import {Box, SPRING, SPRING_MS, autoScrollStep, boxesOf, fadeOut, flip, nearestBox, fitSpan, reducedMotion, scrollerOf,
  snapSpan, stepSpan} from './layout-motion';

/** A panel's width in twelfths: a quarter, a third, a half, two thirds, the row. */
export type Span = 3 | 4 | 6 | 8 | 12;
export const SPANS: Span[] = [3, 4, 6, 8, 12];
const SPAN_WORDS: {[n: number]: [string, string]} = {3: ['¼', 'A quarter of the row'], 4: ['⅓', 'A third of the row'],
  6: ['½', 'Half the row'], 8: ['⅔', 'Two thirds of the row'], 12: ['Full row', 'The whole row']};

/** A panel: what it is (title, what it shows, on its (i)), its content and, for one that folds, the line it shows
 *  while folded. */
export interface Panel {id: string, title: string, hint: string | HTMLElement, content: HTMLElement[],
  summary?: ReadonlySignal<string>}

/** Where a panel goes: its id and its width; `hidden`: left out; `folded` (true or false) makes it a panel that
 *  folds, folded by default when true. */
export interface Slot {id: string, span: Span, hidden?: boolean, folded?: boolean}

/** What follows a panel's summary: the view's form or host. */
export interface Owner {effect: (fn: () => void) => unknown}

const SETTINGS = 'PmxStudio';
const keyOf = (stage: string) => `layout-${stage}`;
const foldKey = (stage: string) => `fold-${stage}`;

/** A stored layout over the stage's panels: unknown panels dropped, panels it does not name (added since) at the
 *  end as the stage places them; null: none stored, or unreadable. */
function readLayout(text: unknown, defaults: Slot[]): Slot[] | null {
  if (typeof text !== 'string' || text.trim() === '')
    return null;
  try {
    const x = JSON.parse(text);
    if (!Array.isArray(x))
      return null;
    const known = new Map(defaults.map((s) => [s.id, s]));
    const kept: Slot[] = x.filter((s) => known.has(s?.id) && SPANS.includes(s.span))
      .map((s) => ({...known.get(s.id), id: s.id, span: s.span, hidden: !!s.hidden}));
    const named = new Set(kept.map((s) => s.id));
    return [...kept, ...defaults.filter((s) => !named.has(s.id))];
  } catch {
    return null;
  }
}

/** The folds a user left (panel id: folded), {} when none. */
function readFolds(stage: string): {[id: string]: boolean} {
  try {
    const x = JSON.parse(grok.userSettings.getValue(SETTINGS, foldKey(stage)) ?? '{}');
    return x && typeof x === 'object' && !Array.isArray(x) ? x : {};
  } catch {
    return {};
  }
}

const labLayout = (stage: string, d: Slot[]) => readLayout(_package.settings?.[keyOf(stage)], d);
const userLayout = (stage: string, d: Slot[]) => readLayout(grok.userSettings.getValue(SETTINGS, keyOf(stage)), d);
const stored = (l: Slot[]) => JSON.stringify(l.map((s) => ({id: s.id, span: s.span, hidden: !!s.hidden})));
const same = (a: Slot[], b: Slot[]) => stored(a) === stored(b);

/** A panel as the page's index lists it: `tab` the tab it is on (panelTabs), `folded` null for one that does not
 *  fold; `go` opens its tab, unfolds it and scrolls to it, `show` shows a hidden one again. */
export interface PageEntry {id: string, title: string, tab?: string, hidden: boolean, folded: boolean | null,
  go: () => void, show: () => void}

/** A page of panels (panelGrid, panelTabs): its index, its changes, and customizing it in place. */
export interface Page {
  entries(): PageEntry[];
  /** `f` runs after every change of the page's layout or folds. */
  onChange(f: () => void): void;
  customizing(): boolean;
  /** Starts customizing (true), or ends it keeping the changes (false). */
  customize(on: boolean): void;
  cancel(): void;
  /** Takes back the last change made while customizing (Ctrl+Z); false: none left. */
  undo(): boolean;
  reset(): void;
  saveForLab(): Promise<void>;
}
const pages = new WeakMap<HTMLElement, Page>();

/** The page a panelGrid or panelTabs root is. */
export const pageOf = (root: HTMLElement): Page | undefined => pages.get(root);

/** Toggles customizing the page `root` is in place (the toolbox's and File's Customize this page). */
export function customizePage(root: HTMLElement): void {
  const p = pages.get(root);
  p?.customize(!p.customizing());
}

/** Whether a key press is the page's Ctrl+Z (Cmd+Z) while it is customized and `shown`: not Ctrl+Shift+Z, not
 *  typing in a field (the field's own undo). */
export function undoKey(e: {ctrlKey: boolean, metaKey: boolean, shiftKey: boolean, altKey?: boolean, key: string,
  target: EventTarget | null}, shown: boolean): boolean {
  return shown && (e.ctrlKey || e.metaKey) && !e.shiftKey && !e.altKey && e.key?.toLowerCase() === 'z' &&
    !(e.target as HTMLElement | null)?.closest?.('input, textarea, select, [contenteditable="true"]');
}

const button = (text: string, title: string, cls: string, run: () => void) => {
  const b = document.createElement('button');
  b.type = 'button';
  b.className = cls;
  b.textContent = text;
  b.title = title;
  b.addEventListener('click', (e) => {
    e.stopPropagation();
    run();
  });
  return b;
};

/** The panels of `layout`, in its order and widths (a panel it leaves out is not shown); with `stage`, the user's
 *  layout or the lab's in its place, the folds the user left, and the page customized in place. `bar`: false when a
 *  page around it (panelTabs) shows the customizing bar. */
export function panelGrid(panels: Panel[], layout: Slot[], cls: string, stage?: string, owner?: Owner,
  bar = true): HTMLElement {
  const defaults = layout.filter((s) => panels.some((p) => p.id === s.id));
  const foldable = new Set(defaults.filter((s) => s.folded !== undefined).map((s) => s.id));
  const listeners: (() => void)[] = [];
  const changed = () => listeners.forEach((f) => f());
  const folds = new Map<string, boolean>();
  const defaultFolds = () => {
    folds.clear();
    defaults.filter((s) => foldable.has(s.id)).forEach((s) => folds.set(s.id, !!s.folded));
  };
  defaultFolds();
  if (stage)
    Object.entries(readFolds(stage)).forEach(([id, f]) => foldable.has(id) && folds.set(id, !!f));
  const titles = new Map(panels.map((p) => [p.id, p.title]));
  const grid = div([], 'pmx-grid');
  let current = defaults;
  const apply = (l: Slot[], quiet = false) => {
    current = l;
    grid.replaceChildren(...l.filter((s) => !s.hidden).map((s) => {
      const c = cells.get(s.id)!;
      // (shown again while it fades out: the fade lets go of it, fully visible)
      if (c.parentElement === host)
        c.getAnimations?.().forEach((a) => a.cancel());
      c.dataset.span = String(s.span);
      const edge = c.querySelector<HTMLElement>('.pmx-resize')!;
      edge.setAttribute('aria-valuenow', String(s.span));
      // (on a narrow page, how it shows too, as the bar's note says: QA 0.1.3 round 4 technician R4-11)
      const w = grid.clientWidth;
      const shown = w > 640 && w <= 980 ? s.span < 6 ? ', shown as half on this narrow page' : s.span === 8 ?
        ', shown as the whole row on this narrow page' : '' : '';
      edge.setAttribute('aria-valuetext', SPAN_WORDS[s.span][1] + shown);
      return c;
    }));
    drawHidden();
    // the charts follow their panel's new width
    if (!quiet)
      window.dispatchEvent(new Event('resize'));
    changed();
  };
  // the page shown changed on the spring: the panels slide from where they were (the lifted one follows the pointer)
  const visible = (): HTMLElement[] => current.filter((s) => !s.hidden).map((s) => cells.get(s.id)!);
  const animate = (l: Slot[], quiet = false, first = boxesOf(visible())) => {
    apply(l, quiet);
    flip(visible().filter((el) => el !== lifted), first);
  };
  // what customizing changed, newest last, for Ctrl+Z
  const history: Slot[][] = [];
  const copy = () => current.map((s) => ({...s}));
  const remember = () => host.classList.contains('is-customizing') && history.push(copy());
  // a panel's box as laid out, whatever it is drawn at meanwhile (a transform: lifted, or sliding home)
  const laidOut = (el: HTMLElement): Box => {
    const g = grid.getBoundingClientRect();
    return {left: g.left + el.offsetLeft, top: g.top + el.offsetTop, width: el.offsetWidth, height: el.offsetHeight};
  };
  // `id` placed `at` among the panels shown (the hidden keep their places among the rest)
  const placed = (id: string, at: number, span?: Span) => {
    const rest = current.filter((s) => s.id !== id);
    const shown = rest.filter((s) => !s.hidden);
    const was = current.find((s) => s.id === id)!;
    const me = span ? {...was, span} : was;
    const i = at >= shown.length ? rest.length : rest.indexOf(shown[at]);
    return [...rest.slice(0, i), me, ...rest.slice(i)];
  };
  let lifted: HTMLElement | null = null;
  // a panel dragged by its header: it lifts and follows the pointer a touch behind; the others slide out of the way as
  // the nearest slot changes; released anywhere, it settles into the slot shown. Escape puts everything back
  const drag = (id: string, cell: HTMLElement, e0: PointerEvent) => {
    const start = copy();
    const instant = reducedMotion();
    const home0 = laidOut(cell);
    const grab = {x: e0.clientX - home0.left, y: e0.clientY - home0.top};
    const pointer = {x: e0.clientX, y: e0.clientY};
    const pos = {x: home0.left, y: home0.top};
    const scroller = scrollerOf(grid);
    lifted = cell;
    cell.classList.add('is-lifted');
    let raf = 0;
    // the slot under the pointer now, and the width the panel takes there: its own, or the gap left in that row (it
    // springs back to its own away from the gap); taken once it holds for a few frames and the last move has settled
    // (a big panel passing a small one would otherwise swing back and forth between two slots)
    const own = start.find((s) => s.id === id)!.span;
    const spanOf = (el: HTMLElement) => current.find((s) => s.id === el.dataset.panel)!.span;
    const slotNow = (): {at: number, span: Span} => {
      const others = visible().filter((el) => el !== cell);
      const boxes = others.map(laidOut);
      const n = nearestBox(boxes, pointer.x, pointer.y);
      if (n.i < 0)
        return {at: 0, span: own};
      const row = others.filter((_, k) => Math.abs(boxes[k].top - boxes[n.i].top) < 1).map(spanOf);
      return {at: n.i + (n.after ? 1 : 0), span: fitSpan(row, own) as Span};
    };
    const differs = (s: {at: number, span: Span}) => visible().indexOf(cell) !== s.at || spanOf(cell) !== s.span;
    let [seen, held, movedAt] = ['', 0, 0];
    const frame = () => {
      const s = scroller.getBoundingClientRect();
      const step = autoScrollStep(pointer.y, Math.max(0, s.top), Math.min(innerHeight, s.bottom));
      if (step)
        scroller.scrollTop += step;
      const k = instant ? 1 : 0.35;
      pos.x += (pointer.x - grab.x - pos.x) * k;
      pos.y += (pointer.y - grab.y - pos.y) * k;
      const slot = slotNow();
      const key = `${slot.at}:${slot.span}`;
      held = key === seen ? held + 1 : 0;
      seen = key;
      if (differs(slot) && (instant || (held >= 2 && performance.now() - movedAt > SPRING_MS / 2))) {
        animate(placed(id, slot.at, slot.span), true);
        movedAt = performance.now();
      }
      const home = laidOut(cell);
      cell.style.transform = `translate(${pos.x - home.left}px, ${pos.y - home.top}px) scale(1.02)`;
      raf = requestAnimationFrame(frame);
    };
    const onMove = (e: PointerEvent) => {
      pointer.x = e.clientX;
      pointer.y = e.clientY;
    };
    const end = (cancel: boolean) => {
      cancelAnimationFrame(raf);
      removeEventListener('pointermove', onMove);
      removeEventListener('pointerup', onUp);
      removeEventListener('pointercancel', onUp);
      removeEventListener('keydown', onKey, true);
      // released: where the pointer is, even before a frame caught up with a fast throw
      const slot = cancel ? null : slotNow();
      if (!slot)
        animate(start, true);
      else if (differs(slot))
        animate(placed(id, slot.at, slot.span), true);
      const home = laidOut(cell);
      const from = `translate(${pos.x - home.left}px, ${pos.y - home.top}px) scale(1.02)`;
      cell.style.transform = '';
      lifted = null;
      const done = () => cell.classList.remove('is-lifted');
      if (instant)
        done();
      else {
        cell.animate([{transform: from}, {transform: 'none'}], {duration: SPRING_MS, easing: SPRING}).finished
          .then(done, done);
      }
      if (!cancel && stored(start) !== stored(current))
        history.push(start);
      window.dispatchEvent(new Event('resize'));
    };
    const onUp = () => end(false);
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        e.preventDefault();
        e.stopPropagation();
        end(true);
      }
    };
    addEventListener('pointermove', onMove);
    addEventListener('pointerup', onUp);
    addEventListener('pointercancel', onUp);
    addEventListener('keydown', onKey, true);
    raf = requestAnimationFrame(frame);
  };
  // a press on a header while customizing: a mouse or pen drags once it moves a little; a finger holds a moment
  // first (a quick swipe still scrolls the page)
  const press = (id: string, cell: HTMLElement, e: PointerEvent) => {
    if (!host.classList.contains('is-customizing') || e.button !== 0 || lifted ||
      (e.target as HTMLElement).closest('.pmx-hide, .pmx-info, a, input, textarea, select'))
      return;
    const touch = e.pointerType === 'touch';
    const [x0, y0] = [e.clientX, e.clientY];
    let timer = 0;
    const stop = () => {
      clearTimeout(timer);
      removeEventListener('pointermove', onMove);
      removeEventListener('pointerup', stop);
      removeEventListener('pointercancel', stop);
    };
    const go = (ev: PointerEvent) => {
      stop();
      drag(id, cell, ev);
    };
    const onMove = (ev: PointerEvent) => {
      const d = Math.hypot(ev.clientX - x0, ev.clientY - y0);
      if (touch ? d > 8 : d > 4)
        touch ? stop() : go(ev);
    };
    if (touch)
      timer = window.setTimeout(() => go(e), 350);
    addEventListener('pointermove', onMove);
    addEventListener('pointerup', stop);
    addEventListener('pointercancel', stop);
  };
  // a panel's right edge dragged: the panel follows the pointer freely, over its neighbours, while a ghost shows the
  // stop it will snap to (a quarter, a third, a half, two thirds, the row) and names it; released, it settles there on
  // the spring, the others reflowing with it (with reduced motion: the ghost alone, then the snap at once)
  const ghost = div([], 'pmx-snap-ghost');
  const resize = (id: string, cell: HTMLElement, e0: PointerEvent) => {
    e0.preventDefault();
    const follow = !reducedMotion();
    const gap = parseFloat(getComputedStyle(grid).columnGap) || 0;
    const col = () => (grid.getBoundingClientRect().width + gap) / 12;
    let span = Number(cell.dataset.span) as Span;
    const show = () => {
      Object.assign(ghost.style, {left: `${cell.offsetLeft}px`, top: `${cell.offsetTop}px`,
        width: `${col() * span - gap}px`, height: `${cell.offsetHeight}px`});
      ghost.dataset.label = SPAN_WORDS[span][0];
    };
    grid.append(ghost);
    cell.classList.add('is-resizing');
    show();
    const onMove = (e: PointerEvent) => {
      const left = laidOut(cell).left;
      const w = e.clientX - left;
      span = snapSpan((w + gap) / col());
      show();
      if (follow) // (from a little under a quarter to the grid's right edge)
        cell.style.width = `${Math.max(col() * 2 - gap, Math.min(grid.getBoundingClientRect().right - left, w))}px`;
    };
    const onUp = (e: PointerEvent) => {
      removeEventListener('pointermove', onMove);
      removeEventListener('pointerup', onUp);
      removeEventListener('pointercancel', onUp);
      ghost.remove();
      const first = boxesOf(visible());
      cell.style.width = '';
      cell.classList.remove('is-resizing');
      setSpan(id, e.type === 'pointercancel' ? Number(cell.dataset.span) as Span : span, first);
    };
    addEventListener('pointermove', onMove);
    addEventListener('pointerup', onUp);
    addEventListener('pointercancel', onUp);
  };
  // Ctrl+Z while customizing the page on screen (a hidden tab's grid lets the open one take it): ours alone, caught
  // first (capture) and kept from the platform even with nothing left to undo, whose own Ctrl+Z opened a second
  // Overview (QA 0.1.3)
  const onUndo = (e: KeyboardEvent) => {
    // Esc leaves customizing as Cancel does, the session's changes put back; during a drag, the drag's own Esc puts
    // that drag back first
    if (e.key === 'Escape' && !lifted && host.offsetParent !== null &&
      !(e.target as HTMLElement | null)?.closest?.('input, textarea, select, [contenteditable="true"]')) {
      e.preventDefault();
      e.stopImmediatePropagation();
      page.cancel();
      return;
    }
    if (!undoKey(e, host.offsetParent !== null))
      return;
    e.preventDefault();
    e.stopImmediatePropagation();
    page.undo();
  };
  const setFold = (id: string, f: boolean, keep = true) => {
    folds.set(id, f);
    const c = cells.get(id)!;
    c.classList.toggle('is-folded', f);
    c.querySelector('.pmx-block-header')?.setAttribute('aria-expanded', String(!f));
    if (!f) // a chart laid out folded follows its width once shown
      setTimeout(() => window.dispatchEvent(new Event('resize')));
    if (keep && stage) {
      grok.userSettings.add(SETTINGS, foldKey(stage), JSON.stringify(Object.fromEntries(folds)));
      void grok.userSettings.flush();
    }
    changed();
  };
  // a slot changed in place: the page shows it at once, on the spring (kept on Done, taken back by Ctrl+Z)
  const edit = (id: string, f: (s: Slot) => Slot, first?: Map<HTMLElement, Box>) => {
    remember();
    animate(current.map((s) => s.id === id ? f({...s}) : s), false, first);
  };
  // a panel's width set by its edge (drag, keys, double-click); the row's full width toggles back to the one before
  const widthWas = new Map<string, Span>();
  const setSpan = (id: string, span: Span, first?: Map<HTMLElement, Box>) => {
    const was = current.find((s) => s.id === id)!.span;
    if (span === was)
      return first && flip(visible(), first);
    if (span === 12)
      widthWas.set(id, was);
    edit(id, (s) => ({...s, span}), first);
  };
  const move = (id: string, to: number) => {
    const work = [...current];
    const from = work.findIndex((x) => x.id === id);
    const [m] = work.splice(from, 1);
    work.splice(Math.max(0, Math.min(work.length, to)), 0, m);
    remember();
    animate(work);
  };
  const cells = new Map(panels.map((p) => {
    const b = block(p.title, p.hint, ...p.content);
    const cell = div([b], 'pmx-panel');
    cell.dataset.panel = p.id;
    const header = b.querySelector<HTMLElement>('.pmx-block-header')!;
    // what customizing offers on the header: a handle to drag it (arrows move it too), its widths, × hides it
    const handle = button('⠿', 'Drag the header to move this panel (or arrows)', 'pmx-handle', () => {});
    header.addEventListener('pointerdown', (e) => press(p.id, cell, e));
    const edge = div([], 'pmx-resize');
    edge.title = 'Drag to change the width (or ← →); double-click: the full row and back';
    Object.entries({'tabindex': '0', 'role': 'separator', 'aria-orientation': 'vertical', 'aria-valuemin': '3',
      'aria-valuemax': '12', 'aria-label': `Width of ${p.title}`}).forEach(([k, v]) => edge.setAttribute(k, v));
    edge.addEventListener('pointerdown', (e) => host.classList.contains('is-customizing') && e.button === 0 &&
      resize(p.id, cell, e));
    edge.addEventListener('dblclick', () => {
      const now = current.find((s) => s.id === p.id)!.span;
      setSpan(p.id, now === 12 ? widthWas.get(p.id) ?? 6 : 12);
    });
    edge.addEventListener('keydown', (e) => {
      const now = current.find((s) => s.id === p.id)!.span;
      const to = e.key === 'ArrowLeft' ? stepSpan(now, -1) : e.key === 'ArrowRight' ? stepSpan(now, 1) :
        e.key === 'Home' ? 3 : e.key === 'End' ? 12 : 0;
      if (!to)
        return;
      e.preventDefault();
      setSpan(p.id, to);
      edge.focus(); // (the grid laid out again lets go of it)
    });
    cell.append(edge);
    handle.addEventListener('keydown', (e) => {
      const step = e.key === 'ArrowUp' || e.key === 'ArrowLeft' ? -1 :
        e.key === 'ArrowDown' || e.key === 'ArrowRight' ? 1 : 0;
      if (!step)
        return;
      e.preventDefault();
      move(p.id, current.findIndex((s) => s.id === p.id) + step);
      handle.focus();
    });
    const tools = span('', 'pmx-panel-tools');
    tools.append(button('×', 'Hide this panel', 'pmx-hide', () => {
      if (current.find((s) => s.id === p.id)?.hidden)
        return;
      const [h, b] = [host.getBoundingClientRect(), cell.getBoundingClientRect()];
      edit(p.id, (s) => ({...s, hidden: true}));
      fadeOut(cell, host, {left: b.left - h.left, top: b.top - h.top, width: b.width, height: b.height});
    }));
    header.prepend(handle);
    header.append(tools);
    if (foldable.has(p.id)) {
      cell.classList.add('pmx-foldable');
      header.insertBefore(span('', 'pmx-chevron'), handle.nextSibling);
      const sum = span('', 'pmx-panel-summary');
      header.insertBefore(sum, tools);
      if (p.summary && owner)
        owner.effect(() => sum.textContent = p.summary!.value);
      header.tabIndex = 0;
      header.setAttribute('role', 'button');
      // its name is its title: as a button it would otherwise be all its text, the (i) and the summary that repeats
      // it (QA 0.1.3 final designer #29); the summary is its description
      header.setAttribute('aria-label', p.title);
      if (p.summary && owner) {
        owner.effect(() => p.summary!.value ? header.setAttribute('aria-description', p.summary!.value) :
          header.removeAttribute('aria-description'));
      }
      const toggle = () => !grid.closest('.is-customizing') && setFold(p.id, !folds.get(p.id));
      header.addEventListener('click', (e) => {
        if (!(e.target as HTMLElement).closest('button, a, .pmx-info'))
          toggle();
      });
      header.addEventListener('keydown', (e) => {
        if ((e.key === 'Enter' || e.key === ' ') && e.target === header) {
          e.preventDefault();
          toggle();
        }
      });
    }
    return [p.id, cell];
  }));
  folds.forEach((f, id) => setFold(id, f, false));
  // the panels hidden: one quiet line at the end of the page, its "show" listing them
  const hiddenLine = div([], 'pmx-hidden-line pmx-noprint');
  const showPanel = (id: string) => {
    const l = current.map((s) => s.id === id ? {...s, hidden: false} : s);
    if (host.classList.contains('is-customizing')) {
      remember();
      animate(l);
    } else {
      // (grown back in its place, the rest sliding to make room)
      const els = visible();
      const first = boxesOf(els);
      saveLayout(l);
      flip(visible(), first);
    }
  };
  const drawHidden = () => {
    const hidden = current.filter((s) => s.hidden);
    hiddenLine.hidden = !hidden.length;
    const show = document.createElement('a');
    show.textContent = 'show';
    show.className = 'pmx-hint-link';
    show.addEventListener('click', (e) => {
      const menu = DG.Menu.popup();
      for (const s of current.filter((x) => x.hidden))
        menu.item(`+ ${titles.get(s.id) ?? s.id}`, () => showPanel(s.id));
      menu.show({causedBy: e});
    });
    hiddenLine.replaceChildren(`${hidden.length} hidden · `, show);
    hiddenLine.title = hidden.map((s) => titles.get(s.id)).join(', ');
  };
  const saveLayout = (l: Slot[]) => {
    if (stage) {
      // the lab's layout (or the stage's) kept as it is is no layout of one's own: a later lab change reaches it
      if (same(l, labLayout(stage, defaults) ?? defaults))
        grok.userSettings.delete(SETTINGS, keyOf(stage));
      else
        grok.userSettings.add(SETTINGS, keyOf(stage), stored(l));
      // saved now, not at the next sync (every 10 s): a reload right after keeps it
      void grok.userSettings.flush();
    }
    apply(l);
  };
  const host = div([grid, hiddenLine], `pmx-grid-host ${cls}`);
  let before: Slot[] | null = null;
  const page: Page = {
    entries: () => current.map((s) => ({id: s.id, title: titles.get(s.id) ?? s.id, hidden: !!s.hidden,
      folded: foldable.has(s.id) ? !!folds.get(s.id) : null,
      go: () => {
        if (folds.get(s.id))
          setFold(s.id, false);
        // (after the tab it is on is drawn)
        setTimeout(() => cells.get(s.id)!.scrollIntoView({behavior: 'smooth', block: 'start'}));
      },
      show: () => showPanel(s.id)})),
    onChange: (f) => void listeners.push(f),
    customizing: () => host.classList.contains('is-customizing'),
    customize: (on) => {
      if (on === page.customizing())
        return;
      host.classList.toggle('is-customizing', on);
      history.length = 0;
      if (on) {
        before = current.map((s) => ({...s}));
        addEventListener('keydown', onUndo, true);
      } else
        removeEventListener('keydown', onUndo, true);
      if (!on && before) {
        before = null;
        saveLayout(current);
      }
      changed();
    },
    cancel: () => {
      if (before)
        animate(before);
      before = null;
      history.length = 0;
      removeEventListener('keydown', onUndo, true);
      host.classList.remove('is-customizing');
      changed();
    },
    undo: () => {
      const l = history.pop();
      if (l)
        animate(l);
      return !!l;
    },
    reset: () => {
      if (stage) {
        grok.userSettings.delete(SETTINGS, foldKey(stage));
        void grok.userSettings.flush();
      }
      defaultFolds();
      folds.forEach((f, id) => setFold(id, f, false));
      apply((stage ? labLayout(stage, defaults) : null) ?? defaults.map((s) => ({...s})));
    },
    saveForLab: async () => {
      if (!stage)
        return;
      await saveLabSetting(keyOf(stage), stored(current));
      grok.userSettings.delete(SETTINGS, keyOf(stage));
      void grok.userSettings.flush();
    },
  };
  pages.set(host, page);
  apply(stage ? userLayout(stage, defaults) ?? labLayout(stage, defaults) ?? defaults : defaults);
  if (bar && stage)
    host.prepend(customizeBar(page));
  return host;
}

/** What the customize bar says of a page `w` px wide (the grid host's content box), as studio.css's container queries
 *  lay it out: under 980 px ¼ and ⅓ show as ½ and ⅔ as the full row; under 640 every panel is the full row. */
export const narrowNote = (w: number): string => w <= 0 ? '' : w <= 640 ? 'Widths apply when the page is wider' :
  w <= 980 ? 'On a page this narrow, ¼ and ⅓ show as ½, ⅔ as the full row' : '';
const contentWidth = (e: HTMLElement): number => {
  const st = getComputedStyle(e);
  return e.clientWidth - parseFloat(st.paddingLeft || '0') - parseFloat(st.paddingRight || '0');
};

/** The bar over a page being customized: what to do, Reset, Save as the lab's (for whoever may), Cancel, Done. */
function customizeBar(page: Page): HTMLElement {
  const lab = button('Save as the lab\'s', 'Every user\'s layout who has none of their own', 'pmx-bar-button', () => {
    void page.saveForLab().then(() => {
      page.customize(false);
      grok.shell.info('Saved as the lab\'s layout');
    }).catch((e: any) => grok.shell.error(String(e?.message ?? e)));
  });
  lab.hidden = true;
  void canEditLabRules().then((can) => lab.hidden = !can).catch(() => {});
  const help = span('Drag a header to move · drag the right edge to resize · × hides · Ctrl+Z undoes', 'pmx-bar-text');
  // (never narrower than a line of words: with the context panel open it wrapped one word a line beside the buttons;
  // short of room, it takes a line of its own over them, the bar wrapping, QA 0.1.3 browser check, designer U-a)
  Object.assign(help.style, {flex: '1 1 260px', minWidth: '0'});
  // a page too narrow for the widths chosen (studio.css container queries: ¼ and ⅓ show as ½ under 980 px, every panel
  // the full row under 640): said, not a width springing back with no reason (QA round 3 designer B5)
  const narrow = span('', 'pmx-bar-note');
  // (its look inline: the package CSS is capped; the grip clear of the scrollbar and the ghost's label above the
  // panel being resized are studio.css's, QA round 3 designer U7, B5)
  Object.assign(narrow.style, {fontSize: '12px', color: 'var(--dg-grey-5, #717581)'});
  const done = button('Done', 'Keep this layout, yours on every study', 'pmx-bar-button is-primary',
    () => page.customize(false));
  const bar = div([help, narrow, button('Reset', 'Back to the lab\'s layout, or the stage\'s own, and its folds',
    'pmx-bar-button', () => page.reset()), lab, button('Cancel', 'Put the page back as it was', 'pmx-bar-button',
    () => page.cancel()), done], 'pmx-customize-bar pmx-noprint');
  bar.style.flexWrap = 'wrap';
  // (hidden while not customizing: width 0, no note; ponytail: observed for the page's life, one per stage view). The
  // width the container queries read, the grid host's content box, not the bar's: the note said ½ while Home still
  // gave a real ¼ (QA 0.1.3 round 4 designer U-b)
  // (on the next frame: observing and writing the note inside the callback is a "ResizeObserver loop" error)
  const seen = new ResizeObserver(() => requestAnimationFrame(() => {
    const host = bar.parentElement?.classList.contains('pmx-grid-host') ? bar.parentElement :
      [...bar.parentElement?.querySelectorAll<HTMLElement>('.pmx-grid-host') ?? []].find((h) => h.offsetWidth > 0);
    if (host)
      seen.observe(host);
    const w = bar.offsetWidth > 0 && host ? contentWidth(host) : 0;
    narrow.textContent = narrowNote(w);
  }));
  seen.observe(bar);
  return bar;
}

/** A tab: its name and its panels, in their default order and widths. */
export interface Tab {name: string, layout: Slot[]}

export interface Tabbed {
  root: HTMLElement;
  /** Opens the tab `panel` is on. */
  show(panel: string): void;
  /** Runs `f` with every tab shown (so its charts are drawn: the PDF report), then the open tab alone again. */
  all<T>(f: () => Promise<T>): Promise<T>;
}

/** The stage's panels on tabs, each tab a grid of its own with its own layout (`${stage}-${tab}`), customized
 *  together. The tabs are a row of their own, under what stays on top (the study, the figures). */
export function panelTabs(panels: Panel[], tabs: Tab[], cls: string, stage: string, owner?: Owner): Tabbed {
  const grids = tabs.map((t) => panelGrid(panels.filter((p) => t.layout.some((s) => s.id === p.id)), t.layout, cls,
    `${stage}-${t.name.toLowerCase().replace(/\W+/g, '-')}`, owner, false));
  const heads = tabs.map((t, i) => {
    const b = document.createElement('button');
    b.className = 'pmx-tab';
    b.textContent = t.name;
    b.setAttribute('role', 'tab');
    b.onclick = () => openAtTop(i);
    return b;
  });
  // a tab picked by the reader starts at its top: opened with the page scrolled down, it showed the blank space the
  // previous tab's height left, its panels above out of sight
  const openAtTop = (i: number) => {
    open(i);
    const scroller = scrollerOf(root);
    const above = strip.getBoundingClientRect().top - (scroller === document.scrollingElement ? 0 :
      scroller.getBoundingClientRect().top);
    if (above < 0)
      scroller.scrollTop += above;
  };
  const strip = div(heads, 'pmx-tabs pmx-noprint');
  strip.setAttribute('role', 'tablist');
  const root = div([strip, ...grids], 'pmx-tabbed');
  let current = 0;
  const open = (i: number) => {
    // a shorter tab opened while the page is scrolled down would shrink the page under the reader, and the browser
    // pulls everything down to keep the scroll in range: the tabs keep the height they had, the new tab's space below
    // ponytail: the height only grows while the stage is open; reset it if a short tab's blank tail ever bothers
    root.style.minHeight = `${root.offsetHeight}px`;
    current = i;
    grids.forEach((g, j) => g.hidden = j !== i);
    heads.forEach((h, j) => {
      h.classList.toggle('is-open', j === i);
      h.setAttribute('aria-selected', String(j === i));
      h.tabIndex = j === i ? 0 : -1;
    });
    window.dispatchEvent(new Event('resize'));
  };
  // arrows move between the tabs, as in any tab list
  strip.addEventListener('keydown', (e) => {
    const step = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
    if (!step)
      return;
    const i = (current + step + tabs.length) % tabs.length;
    openAtTop(i);
    heads[i].focus();
  });
  open(0);
  const each = grids.map((g) => pages.get(g)!);
  const showTab = (panel: string) => {
    const i = tabs.findIndex((t) => t.layout.some((s) => s.id === panel));
    if (i >= 0 && i !== current)
      open(i);
  };
  const page: Page = {
    entries: () => each.flatMap((p, i) => p.entries().map((e) => ({...e, tab: tabs[i].name,
      go: () => {
        showTab(e.id);
        e.go();
      }}))),
    onChange: (f) => each.forEach((p) => p.onChange(f)),
    customizing: () => each[0].customizing(),
    customize: (on) => {
      each.forEach((p) => p.customize(on));
      root.classList.toggle('is-customizing', on);
    },
    cancel: () => {
      each.forEach((p) => p.cancel());
      root.classList.remove('is-customizing');
    },
    undo: () => each[current].undo(),
    reset: () => each.forEach((p) => p.reset()),
    saveForLab: async () => {
      for (const p of each)
        await p.saveForLab();
    },
  };
  pages.set(root, page);
  root.prepend(customizeBar(page));
  return {root, show: showTab,
    all: async (f) => {
      const was = current;
      grids.forEach((g) => g.hidden = false);
      window.dispatchEvent(new Event('resize'));
      // the charts redraw at their width before they are pictured
      await new Promise((r) => setTimeout(r, 300));
      try {
        return await f();
      } finally {
        open(was);
        // (not the height of every tab at once)
        root.style.minHeight = '';
      }
    }};
}

/** The page's index for the toolbox: each panel a link (its tab first, on a tabbed page) that opens its tab,
 *  unfolds it and scrolls to it; a hidden one greyed, with "+ show". It follows the page's changes. */
export function pageIndex(root: HTMLElement): HTMLElement {
  const page = pages.get(root);
  const box = div([], 'pmx-index');
  if (!page)
    return box;
  const draw = () => {
    const rows: HTMLElement[] = [];
    let tab: string | undefined;
    for (const e of page.entries()) {
      if (e.tab !== undefined && e.tab !== tab) {
        tab = e.tab;
        const first = e;
        const t = document.createElement('a');
        t.className = 'pmx-index-tab';
        t.textContent = tab;
        t.addEventListener('click', () => first.go());
        rows.push(t);
      }
      // (a panel named as its tab is the tab's link)
      if (e.tab !== undefined && e.title === e.tab && !e.hidden)
        continue;
      const row = div([], `pmx-index-item${e.tab !== undefined ? ' is-nested' : ''}${e.hidden ? ' is-hidden' : ''}`);
      if (e.hidden) {
        const plus = document.createElement('a');
        plus.textContent = '+ show';
        plus.className = 'pmx-index-show';
        plus.addEventListener('click', () => e.show());
        row.append(span(e.title), plus);
        row.title = 'Hidden on this page';
      } else {
        const a = document.createElement('a');
        a.textContent = e.title;
        a.addEventListener('click', () => e.go());
        if (e.folded)
          a.title = 'Folded: opens it';
        row.append(a);
      }
      rows.push(row);
    }
    box.replaceChildren(...rows);
  };
  page.onChange(draw);
  draw();
  return box;
}
