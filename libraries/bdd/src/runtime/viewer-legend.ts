/* The legend a viewer hosts — read from the `data-legend-*` attributes it publishes on every
   commit, never from its render path — and the row tooltip, which is one element hidden between
   hovers (so every read matches the visible one only). */
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from './patience.js';
import type {ElementRef} from './args.js';
import {withKeys} from './gestures.js';
import {exactText} from './locate.js';
import {LegendState} from './viewer-runtime.js';
import {near, parseHex} from './viewer-pixels.js';
import {onViewer, settle, snapshot, viewerLocator} from './viewers.js';

const TOOLTIP_ROWS = '.d4-tooltip:visible table.d4-row-tooltip-table tr';
const TOOLTIP_COLUMNS = TOOLTIP_ROWS + ' td:first-child';
const LEGEND_MODES: Record<string, string> = {docked: 'docked', corner: 'corner', 'mini icon': 'miniIcon', tooltip: 'tooltip', hidden: 'hidden'};

/** The legend the viewer hosts, wherever it sits now (in the viewer, or re-parented into the
 * tooltip when collapsed to the mini icon). */
function legendRoot(page: Page, viewer: Locator): Locator {
  return viewer.locator('[name="legend"]').or(page.locator('.d4-tooltip [name="legend"]')).filter({visible: true}).first();
}

/** A legend item by the label it shows (`aria-label`, else the label text). */
export async function legendItem(page: Page, target: ElementRef, label: string): Promise<Locator> {
  const legend = legendRoot(page, await viewerLocator(page, target));
  const items = legend.locator('[name="legend-item"]');
  const byAria = items.filter({has: page.locator(`:scope[aria-label="${label.replace(/"/g, '\\"')}" i]`)});
  const byText = items.filter({has: page.locator('.d4-legend-value', {hasText: exactText(label)})});
  return byAria.or(byText).first();
}

export function legendStateOf(page: Page, target: ElementRef): Promise<LegendState | undefined> {
  return onViewer(page, target, (el) => (window as any).__bdd.legendState(el), undefined);
}

function describeLegend(s: LegendState | undefined): string {
  return s ? `mode ${s.mode || 'none'}, slot ${s.slot || 'none'}, ${s.items} items (${s.keys.length} rendered), ${s.width}x${s.height}` : 'no legend';
}

/** The legend sits on the side the viewer's own slot class names (`d4-legend-left`, …). */
export async function expectLegendSide(page: Page, target: ElementRef, side: string): Promise<void> {
  const want = side.toLowerCase();
  const legend = (await viewerLocator(page, target)).locator('[name="legend"]').filter({visible: true}).first();
  await expect(legend, `the legend of ${target.phrase}`).toBeVisible();
  const classes = (await legend.getAttribute('class') ?? '').split(/\s+/);
  const sides = classes.filter((c) => /^d4-legend-(left|right|top|bottom)$/.test(c)).map((c) => c.slice('d4-legend-'.length));
  expect(sides, `the legend of ${target.phrase} is on the ${sides.join(', ') || 'unknown'} side, not the ${want}`).toContain(want);
}

/** The legend's item total (`data-legend-items`, every section, rendered or not). */
export async function expectLegendItems(page: Page, target: ElementRef, count: number): Promise<void> {
  let last: LegendState | undefined;
  try {
    await expect.poll(async () => (last = await legendStateOf(page, target))?.items, {timeout: pollMs(5000)}).toBe(count);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} lists ${last?.items ?? 'no'} items, not ${count} (${describeLegend(last)})`);
  }
}

export async function expectLegendItemsChange(page: Page, target: ElementRef, compare: 'fewer' | 'same'): Promise<void> {
  let last: {before?: LegendState; now?: LegendState} = {};
  const holds = async (): Promise<boolean | string> => {
    last = await onViewer(page, target, (el) => (window as any).__bdd.legendChange(el), undefined);
    if (!last.before)
      return 'no legend at the snapshot';
    if (!last.now)
      return 'no legend now';
    return compare === 'fewer' ? last.now.items < last.before.items :
      last.now.items === last.before.items && last.now.keys.join('|') === last.before.keys.join('|');
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} does not list ${compare === 'fewer' ? 'fewer items than' : 'the same items as'} before (before ${describeLegend(last.before)}; now ${describeLegend(last.now)})`);
  }
}

export async function expectLegendMode(page: Page, target: ElementRef, mode: string): Promise<void> {
  const want = LEGEND_MODES[mode];
  if (!want)
    throw new Error(`a legend is docked, in a corner, collapsed to the mini icon, shown in the tooltip or hidden — not "${mode}"`);
  let last: LegendState | undefined;
  try {
    await expect.poll(async () => (last = await legendStateOf(page, target))?.mode, {timeout: pollMs(5000)}).toBe(want);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} is not ${mode} (${describeLegend(last)})`);
  }
}

/** In a corner by the mode it publishes, and put there: the legend is laid over the viewer
 * (absolutely positioned inside its root) and anchored by the two edges its corner slot names and
 * by no other, so a slot the placement decided but did not apply fails. The anchor, not the box's
 * geometry: a viewer insets a corner legend past its axes and strips, which can put a
 * bottom-anchored legend above the middle. */
export async function expectLegendInCorner(page: Page, target: ElementRef): Promise<void> {
  await expectLegendMode(page, target, 'corner');
  let last = '';
  const holds = async (): Promise<boolean> => {
    const r: {slot: string; position: string; inRoot: boolean; top: string; bottom: string; left: string; right: string} | undefined =
      await onViewer(page, target, (el) => {
        const b = (window as any).__bdd;
        const v = b.viewerOf(el);
        const l: HTMLElement | null = v.root.querySelector('[name="legend"]');
        return l ? {slot: l.dataset.legendSlot ?? '', position: l.style.position, inRoot: l.parentElement === v.root,
          top: l.style.top, bottom: l.style.bottom, left: l.style.left, right: l.style.right} : undefined;
      }, undefined);
    if (!r)
      return false;
    const anchored = (['top', 'bottom', 'left', 'right'] as const).filter((e) => r[e] !== '');
    last = `slot ${r.slot}, ${r.position || 'static'}${r.inRoot ? ' over the viewer' : ' inside a layout part'}, anchored at ${anchored.join(' and ') || 'no edge'}`;
    const want = {leftTop: 'top,left', leftBottom: 'bottom,left', rightTop: 'top,right', rightBottom: 'bottom,right'}[r.slot];
    return r.position === 'absolute' && r.inRoot && want !== undefined && [...anchored].sort().join(',') === want.split(',').sort().join(',');
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} is not put in the corner it names (${last || 'no legend'})`);
  }
}

export async function expectLegendSlot(page: Page, target: ElementRef, slot: string): Promise<void> {
  const norm = (s: string) => s.toLowerCase().replace(/[^a-z]/g, '');
  let last: LegendState | undefined;
  try {
    await expect.poll(async () => norm((last = await legendStateOf(page, target))?.slot ?? ''), {timeout: pollMs(5000)}).toBe(norm(slot));
  }
  catch {
    throw new Error(`the legend of ${target.phrase} is not in the ${slot} slot (${describeLegend(last)})`);
  }
}

/** Mode, slot and size as at the snapshot before the last change, read once the viewer is quiet:
 * a legend that keeps its slot but jumps or resizes inside it has moved. */
export async function expectLegendPlacedAsBefore(page: Page, target: ElementRef): Promise<void> {
  const r: {before?: LegendState; now?: LegendState} = await onViewer(page, target, async (el) => {
    const b = (window as any).__bdd;
    await b.quiet(b.viewerOf(el));
    return b.legendChange(el);
  }, undefined);
  if (!r.before || !r.now)
    throw new Error(`the legend of ${target.phrase}: ${!r.before ? 'no legend at the snapshot' : 'no legend now'}`);
  const where = (l: LegendState) => `${l.mode}/${l.slot}/${Math.round(l.width)}x${Math.round(l.height)}`;
  expect(where(r.now), `the legend of ${target.phrase} moved (before ${describeLegend(r.before)}; now ${describeLegend(r.now)})`).toBe(where(r.before));
}

/** Drags the splitter between a docked legend and the plot: the baseline is taken first, so "wider
 * / narrower / taller / shorter than before" compare with the legend as it was before the drag. */
export async function dragLegendSplitter(page: Page, target: ElementRef, px: number, direction: string): Promise<void> {
  const d = direction.toLowerCase();
  if (!['left', 'right', 'up', 'down'].includes(d))
    throw new Error(`a splitter is dragged left, right, up or down, not "${direction}"`);
  await settle(page, target);
  await snapshot(page, target);
  const splitter = (await viewerLocator(page, target)).locator('[name="legend-splitter"]').filter({visible: true}).first();
  await splitter.waitFor({state: 'visible', timeout: 5000}).catch(async () => {
    throw new Error(`the legend of ${target.phrase} has no splitter to drag (${describeLegend(await legendStateOf(page, target))}) — only a docked legend has one`);
  });
  const box = await splitter.boundingBox();
  if (!box)
    throw new Error(`the legend splitter of ${target.phrase} has no box`);
  const x = box.x + box.width / 2;
  const y = box.y + box.height / 2;
  const dx = d === 'left' ? -px : d === 'right' ? px : 0;
  const dy = d === 'up' ? -px : d === 'down' ? px : 0;
  await page.mouse.move(x, y);
  await page.mouse.down();
  await page.mouse.move(x + dx / 2, y + dy / 2);
  await page.mouse.move(x + dx, y + dy);
  await page.mouse.up();
  await settle(page, target);
}

export type LegendSizeChange = 'wider' | 'narrower' | 'taller' | 'shorter';

/** The legend's own box against the snapshot before the last change, by more than a rounding
 * pixel, read once the viewer is quiet. */
export async function expectLegendSize(page: Page, target: ElementRef, change: LegendSizeChange): Promise<void> {
  let last: {before?: LegendState; now?: LegendState} = {};
  const holds = async (): Promise<boolean> => {
    last = await onViewer(page, target, async (el) => {
      const b = (window as any).__bdd;
      await b.quiet(b.viewerOf(el));
      return b.legendChange(el);
    }, undefined);
    if (!last.before || !last.now)
      return false;
    const delta = change === 'wider' || change === 'narrower' ? last.now.width - last.before.width : last.now.height - last.before.height;
    return change === 'wider' || change === 'taller' ? delta > 2 : delta < -2;
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} is not ${change} than before (before ${describeLegend(last.before)}; now ${describeLegend(last.now)})`);
  }
}

interface DrawnItem {key: string; empty: boolean; canvas: boolean; ink: number; shaped: boolean; text: boolean}

/** How every item of the legend draws its category — the lists are virtualised, so each section is
 * scrolled through and put back. A structure is a renderer's canvas (`d4-legend-canvas-item`) whose
 * painted pixels span a figure, not one line of text (blank only for the empty category, and not
 * every item blank); text is a label span. */
export async function expectLegendItemsDrawnAs(page: Page, target: ElementRef, as: 'structure' | 'text'): Promise<void> {
  let last = '';
  const holds = async (): Promise<boolean> => {
    const r: {total: number; items: DrawnItem[]} = await onViewer(page, target, async (el) => {
      const b = (window as any).__bdd;
      const v = b.viewerOf(el);
      const l: HTMLElement | null = v.root.querySelector('[name="legend"]') ?? document.querySelector('.d4-tooltip [name="legend"]');
      const seen = new Map<string, DrawnItem>();
      const read = (i: Element): void => {
        const key = i.getAttribute('data-item-key') ?? i.getAttribute('aria-label') ?? '';
        const cv = i.querySelector('canvas.d4-legend-value') as HTMLCanvasElement | null;
        let ink = 0;
        let top = Infinity, bottom = -1, left = Infinity, right = -1;
        if (cv && cv.width > 0 && cv.height > 0) {
          const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
          for (let p = 3; p < data.length; p += 4) {
            if (data[p] === 0 || (data[p - 3] > 245 && data[p - 2] > 245 && data[p - 1] > 245))
              continue;
            ink++;
            const px = (p - 3) / 4;
            const x = px % cv.width, y = Math.floor(px / cv.width);
            top = Math.min(top, y); bottom = Math.max(bottom, y); left = Math.min(left, x); right = Math.max(right, x);
          }
        }
        // one line of text is a band a fifth of the canvas tall; a drawn figure spans most of it
        const shaped = cv !== null && ink > 20 && bottom - top >= cv.height * 0.35 && right - left >= cv.width * 0.25;
        seen.set(key, {key, empty: (i.getAttribute('aria-label') ?? '') === '', ink, shaped,
          canvas: i.classList.contains('d4-legend-canvas-item') && cv !== null, text: i.classList.contains('d4-legend-text-item')});
      };
      const frame = (): Promise<void> => new Promise((res) => requestAnimationFrame(() => requestAnimationFrame(() => res())));
      for (const list of Array.from(l?.querySelectorAll('[name^="legend-list-"]') ?? []) as HTMLElement[]) {
        const start = list.scrollTop;
        const step = Math.max(1, list.clientHeight);
        for (let y = 0; ; y += step) {
          list.scrollTop = y;
          await frame();
          list.querySelectorAll('[name="legend-item"]').forEach(read);
          if (y + step >= list.scrollHeight)
            break;
        }
        list.scrollTop = start;
      }
      return {total: Number(l?.dataset.legendItems ?? 0), items: [...seen.values()]};
    }, undefined);
    const items = r.items;
    last = `${items.length} of ${r.total} items read: ` + items.map((i) => `${i.key}${i.empty ? ' (empty)' : ''}: ` +
      (i.canvas ? `canvas, ${i.ink} px painted${i.shaped ? '' : ' in no figure'}` : i.text ? 'text' : 'neither')).join('; ');
    if (items.length === 0 || items.length !== r.total)
      return false;
    // the empty category goes through the renderer too, and a renderer draws nothing for it
    return as === 'structure' ? items.every((i) => i.canvas && (i.shaped || i.empty && i.ink === 0)) && items.some((i) => i.shaped) :
      items.every((i) => i.text && !i.canvas);
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`not every item in the legend of ${target.phrase} is drawn as ${as === 'structure' ? 'a structure' : 'text'}: ${last || 'no item rendered'}`);
  }
}

/** A click on a legend item (the category filters the viewer): the baseline is taken first and
 * the viewer settles after, so the checks that follow read the state after the repaint. */
export async function clickLegendItem(page: Page, target: ElementRef, label: string, options: {key?: string; cross?: boolean} = {}): Promise<void> {
  await snapshot(page, target);
  const item = await legendItem(page, target, label);
  await item.waitFor({state: 'visible', timeout: 5000}).catch(async () => {
    const shown = await legendRoot(page, await viewerLocator(page, target)).locator('.d4-legend-value').allTextContents();
    throw new Error(`no "${label}" item in the legend of ${target.phrase}; it shows: ${shown.map((s) => s.trim()).filter(Boolean).join(' | ') || 'nothing'}`);
  });
  const what = options.cross ? item.locator('.d4-legend-cross').first() : item.locator('.d4-legend-value').first();
  if (options.cross)
    await item.hover();
  await withKeys(page, options.key ? options.key.split('+') : [], () => what.click());
  await settle(page, target);
}

function cssToHex(color: string): string {
  const m = /rgba?\((\d+),\s*(\d+),\s*(\d+)/.exec(color);
  if (m)
    return '#' + [m[1], m[2], m[3]].map((x) => Number(x).toString(16).padStart(2, '0')).join('').toUpperCase();
  const h = /^#([0-9a-f]{6})$/i.exec(color.trim());
  return h ? '#' + h[1].toUpperCase() : '';
}

/** The color a legend item is drawn in (its inline color, the category's color). */
export async function legendItemColor(page: Page, target: ElementRef, label: string): Promise<string> {
  const item = await legendItem(page, target, label);
  await item.waitFor({state: 'visible', timeout: 5000});
  return cssToHex(await item.evaluate((e) => (e as HTMLElement).style.color || getComputedStyle(e).color));
}

export async function expectLegendItemColor(page: Page, target: ElementRef, label: string, color: string): Promise<void> {
  const want = parseHex(color);
  let last = '';
  try {
    await expect.poll(async () => near(last = await legendItemColor(page, target, label), want), {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${label}" item in the legend of ${target.phrase} is not colored ${want}; it is ${last || 'colored nothing'}`);
  }
}

export async function expectLegendItemsDiffer(page: Page, target: ElementRef, a: string, b: string): Promise<void> {
  let ca = '';
  let cb = '';
  try {
    await expect.poll(async () => {
      ca = await legendItemColor(page, target, a);
      cb = await legendItemColor(page, target, b);
      return ca !== '' && cb !== '' && !near(ca, cb);
    }, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${a}" and "${b}" items in the legend of ${target.phrase} are colored alike (${ca || 'nothing'} and ${cb || 'nothing'})`);
  }
}

// --- the row tooltip ----------------------------------------------------------------------------------

/** The value the row tooltip shows for a column name. */
export async function tooltipValue(page: Page, column: string): Promise<string | undefined> {
  const rows = await page.locator(TOOLTIP_ROWS).evaluateAll((trs) =>
    trs.map((tr) => Array.from(tr.querySelectorAll('td')).map((td) => (td.textContent ?? '').trim())));
  const hit = rows.find((cells) => cells[0]?.toUpperCase() === column.toUpperCase());
  return hit ? hit.slice(1).join(' ').trim() : undefined;
}

export async function expectTooltipValue(page: Page, column: string, value: string): Promise<void> {
  await expect.poll(() => tooltipValue(page, column), {timeout: pollMs(5000), message: `"${column}" in the row tooltip`}).toBe(value);
}

export async function tooltipColumns(page: Page): Promise<string[]> {
  const cells = await page.locator(TOOLTIP_COLUMNS).allTextContents();
  return [...new Set(cells.map((c) => c.trim().toUpperCase()).filter((c) => c.length > 0))].sort();
}

export async function expectTooltipColumns(page: Page, list: string, negate = false): Promise<void> {
  const want = [...new Set(list.split(/\s*,\s*/).map((c) => c.trim().toUpperCase()).filter((c) => c.length > 0))].sort();
  const poll = expect.poll(() => tooltipColumns(page), {timeout: pollMs(5000), message: negate ? 'the tooltip shows exactly these columns' : 'the tooltip does not show these columns'});
  await (negate ? poll.not : poll).toEqual(want);
}
