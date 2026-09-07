import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

declare const grok: any;

export interface Pt { x: number; y: number; }
export interface Rect { x: number; y: number; w: number; h: number; }
export interface Viewport { top: number; bottom: number; height: number; }

export const BOX = '[name="viewer-Box-plot"]';

export function bpProp(page: Page, prop: string): Promise<any> {
  return page.evaluate((p) => grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot')?.props?.[p], prop);
}

export async function setBpProp(page: Page, prop: string, value: any, capMs = 900): Promise<void> {
  await v.setViewerProps(page, 'Box plot', [{set: {[prop]: value}, wait: capMs}]);
}

/**
 * Sets several properties in one go and settles once.
 *
 * For a run of sets whose intermediate states nothing asserts — a step's restore tail, a
 * value/category pair — this is one round trip and one render wait instead of N of each.
 */
export async function setBpProps(page: Page, sets: Record<string, any>, capMs = 900): Promise<void> {
  await v.setViewerProps(page, 'Box plot', [{set: sets, wait: capMs}]);
}

export function canvasRect(page: Page): Promise<Rect> {
  return page.evaluate((sel) => {
    const c = document.querySelector(`${sel} canvas[name="canvas"]`)!.getBoundingClientRect();
    return {x: c.x, y: c.y, w: c.width, h: c.height};
  }, BOX);
}

export function viewportRect(page: Page): Promise<Viewport> {
  return page.evaluate(() => {
    const vp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot').viewport;
    return {top: vp.top, bottom: vp.bottom, height: vp.height};
  });
}

export function readLadder(page: Page): Promise<Record<string, any>> {
  return page.evaluate(() => {
    const p = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot').props;
    return {
      valueColumnName: p.valueColumnName,
      category1ColumnName: p.category1ColumnName,
      category2ColumnName: p.category2ColumnName,
      showMinorCategories: p.showMinorCategories,
      showAllCategories: p.showAllCategories,
      markerColorColumnName: p.markerColorColumnName,
      invertColorScheme: p.invertColorScheme,
      colorMin: p.colorMin,
      colorMax: p.colorMax,
      axisType: p.axisType,
      invertYAxis: p.invertYAxis,
      plotStyle: p.plotStyle,
    };
  });
}

// The viewer parks both toggle icons on every draw, shown or not (box_plot_group_comparison.dart
// _positionToggleIcons): the reveal icon 19px left of the bare p-value text, the cross at the right
// end of the comparison strip. Their stored offsets are the only JS-readable trace of where those
// canvas-drawn regions are.
function iconOrigin(page: Page, iconName: string): Promise<Pt | null> {
  return page.evaluate(({sel, name}) => {
    const icon = document.querySelector(`${sel} [name="${name}"]`) as HTMLElement | null;
    const host = icon?.offsetParent as HTMLElement | null;
    if (!icon || !host || icon.style.left === '') return null;
    const h = host.getBoundingClientRect();
    return {x: h.x + parseFloat(icon.style.left), y: h.y + parseFloat(icon.style.top)};
  }, {sel: BOX, name: iconName});
}

/**
 * A point on the p-value text: the reveal hover area when comparison is off, the strip menu region
 * when on. The context-menu hit test (box_plot_t_test.dart pValueHitTest) starts 10px right of the
 * plot edge, so the point sits 20px in.
 */
export async function pValuePoint(page: Page): Promise<Pt | null> {
  const o = await iconOrigin(page, 'show-group-stats');
  return o ? {x: o.x + 39, y: o.y + 7} : null;
}

/**
 * Hovers the region that reveals a group-comparison toggle icon and returns the icon's centre, or
 * null when it never became visible. Only an ENTER into the region shows the icon
 * (htmlMouseOverVisibility keeps an edge state), so the pointer is parked outside the viewer first:
 * a pointer that never left since the last reveal raises no enter, which is how the Scenario 2
 * setup of the group-comparison spec read false after its click.
 */
export async function revealToggleIcon(page: Page, iconName: 'show-group-stats' | 'close-group-stats'): Promise<Pt | null> {
  await page.mouse.move(0, 0);
  const r = await canvasRect(page);
  const own = await iconOrigin(page, iconName);
  // the icon's own slot, not the p-value text: the text carries the test-name tooltip
  const spots: Pt[] = own
    ? [iconName === 'show-group-stats' ? {x: own.x + 5, y: own.y + 7} : {x: own.x - 20, y: own.y + 7}]
    : [];
  if (spots.length === 0)
    for (const [dx, dy] of [[35, 15], [40, 17], [30, 14], [45, 16]]) spots.push({x: r.x + dx, y: r.y + dy});
  const centre = () => page.evaluate(({sel, name}) => {
    const el = document.querySelector(`${sel} [name="${name}"]`) as HTMLElement | null;
    const b = el?.getBoundingClientRect();
    return el && b && b.width > 0 && getComputedStyle(el).visibility === 'visible'
      ? {x: b.x + b.width / 2, y: b.y + b.height / 2} : null;
  }, {sel: BOX, name: iconName});
  for (const s of spots) {
    await page.mouse.move(s.x, s.y);
    const pt = await v.pollValue(centre, (p) => p !== null, 500, 50);
    if (pt) return pt;
  }
  return null;
}

export async function clickToggleIcon(page: Page, iconName: 'show-group-stats' | 'close-group-stats'): Promise<void> {
  const pt = await revealToggleIcon(page, iconName);
  if (!pt) throw new Error(`${iconName} icon did not appear on hover`);
  await page.mouse.click(pt.x, pt.y);
  const want = iconName === 'show-group-stats';
  await v.pollValue(() => bpProp(page, 'showGroupComparison'), (on) => on === want, 1500, 50);
}

interface SliderRead { rect: Rect; top: Pt; bottom: Pt; }

/** Handle centres of the value-axis range slider; the slider lays out on hovering the axis strip it sits on. */
export async function verticalSliderHandles(page: Page): Promise<{top: Pt; bottom: Pt}> {
  const read = (): Promise<SliderRead | null> => page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const slider = Array.from(bp.root.querySelectorAll('svg[type="range-slider"]'))
      .find((s: any) => { const r = s.getBoundingClientRect(); return r.height > r.width; }) as SVGElement | undefined;
    if (!slider) return null;
    const sr = slider.getBoundingClientRect();
    const centres = Array.from(slider.querySelectorAll('circle')).slice(0, 2)
      .map((c) => { const r = c.getBoundingClientRect(); return {x: r.x + r.width / 2, y: r.y + r.height / 2}; })
      .sort((a, b) => a.y - b.y);
    if (centres.length < 2) return null;
    return {rect: {x: sr.x, y: sr.y, w: sr.width, h: sr.height}, top: centres[0], bottom: centres[1]};
  });
  const laidOut = (s: SliderRead | null) => !!s && s.bottom.y - s.top.y > 50;
  let s = await read();
  if (!laidOut(s)) {
    const r = await canvasRect(page);
    const hover = s ? {x: s.rect.x + s.rect.w / 2, y: s.rect.y + s.rect.h / 2} : {x: r.x + 12, y: r.y + r.h / 2};
    await page.mouse.move(hover.x, hover.y);
    s = await v.pollValue(read, laidOut, 1500, 50);
  }
  if (!laidOut(s)) throw new Error('vertical range slider handles did not lay out to a usable span');
  return {top: s!.top, bottom: s!.bottom};
}

/**
 * Drags the top value-axis handle down by `frac` of the handle span, narrowing the viewport.
 *
 * Three intermediate moves, not sixteen: the slider tracks the pointer on every mousemove, so the
 * extra thirteen only bought thirteen more actionability round trips (2.1s, measured 2026-09-04).
 */
export async function dragTopHandle(page: Page, frac: number): Promise<void> {
  const h = await verticalSliderHandles(page);
  await page.mouse.move(h.top.x, h.top.y);
  await page.mouse.down();
  await page.mouse.move(h.top.x, h.top.y + (h.bottom.y - h.top.y) * frac, {steps: 3});
  await page.mouse.up();
}

/**
 * Resolves once the box plot has painted, and installs the render stamp while doing it.
 *
 * The plain waitForViewerRendered subscribes when it is called, so the very first settle after
 * addViewer has already missed the paint it waits for and burns its whole cap.
 */
export async function bpPainted(page: Page, capMs = 3000): Promise<number> {
  await v.waitForViewerRendered(page, 'Box plot', 0);
  return v.pollValue(() => bpCanvasInk(page), (n) => n > 0, capMs, 50);
}

/** Non-white, non-transparent pixels on the box plot canvas. */
export function bpCanvasInk(page: Page): Promise<number> {
  return page.evaluate(() => {
    const bp = grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot');
    const cv = bp?.root?.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
    if (!cv) return 0;
    const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      const r = data[i], g = data[i + 1], b = data[i + 2], a = data[i + 3];
      if (a !== 0 && !(r >= 250 && g >= 250 && b >= 250)) n++;
    }
    return n;
  });
}

