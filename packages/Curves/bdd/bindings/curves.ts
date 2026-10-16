/* Curves' own screen parts. The datasets are the package's files on the stand (`System:AppData/Curves/...`,
   published with the package). The legend of a fitted-curve chart is drawn on its canvas, and the viewer
   reports no hit areas for it: its rows are found the way a person finds them, by where the chart's own
   hover answers - a synthetic pointer move over the canvas raises the tooltip the chart gives a shortened
   name (one row) or `+N more` (the list), and the rows found are reported as the viewer's hit areas
   `shortened legend name` and `more curves row` until the next measure. Synthetic moves are not filmed. */
import {Page} from '@playwright/test';
import {dataset, Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, expect, locate, silent, viewers} from '@datagrok-libraries/bdd/runtime';
import {hoverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';

declare const grok: any;
declare const DG: any;

dataset('assay-curves', {path: 'System:AppData/Curves/demo/Assay_Curves.csv', aliases: ['Assay Curves'],
  description: '43 dose-response curves in "Fitted Curve", one per row, with "Compound|Assay|Target" naming each'});

/* The legend needs a plot of 325 x 275 px at least (`FitConstants.MIN_LEGEND_PX_*`): beside the grid in a
   help GIF's 800 px window the viewer docked at its default half would draw none. It overlays the
   selected rows only, as in the Assay Curves demo: the current row's curve is one of them, and would be
   named twice. Show Legend is set as the legend is drawn: left alone, the property starts empty and its
   checkbox reads off over a legend that is on, so the first click would change nothing on the plot. */
export const addWideViewer = Given('user adds a MultiCurveViewer viewer of the selected rows taking {int} percent of the view width', async (page: Page, percent: number) => {
  silent(page);
  await page.evaluate(async (ratio) => {
    const tv = grok.shell.tv;
    const v = await tv.dataFrame.plot.fromType('MultiCurveViewer', {showCurrentRowCurve: false, showMouseOverRowCurve: false, showLegend: true});
    tv.dockManager.dock(v, DG.DOCK_TYPE.RIGHT, null, 'Multi curve viewer', ratio);
  }, percent / 100);
  await page.locator('[name="viewer-MultiCurveViewer"]').filter({visible: true}).first().waitFor();
}, {tier: 'api', description: 'a multi curve viewer on the current table, docked to the right of the grid at that share of the view; silent in a guide'});

/* A molecule column is drawn as text until Chem's renderer is up, which on a busy stand happens in the
   middle of a video: a guide that is not about it leaves it out of the grid. */
export const hideGridColumn = Given('the grid of the current view does not show the {string} column', async (page: Page, column: string) => {
  silent(page);
  await page.evaluate((name) => { grok.shell.tv.grid.columns.byName(name).visible = false; }, column);
}, {tier: 'api', description: 'hides the column in the grid of the current table view (the table keeps it); silent in a guide'});

export const SHORTENED = 'shortened legend name';
export const MORE = 'more curves row';

interface LegendFound {areas: Record<string, {x: number; y: number; width: number; height: number}>; name?: string; listed: number;
  width: number; height: number}

/* Page source, not a function: tsx would wrap a named function in a helper the page lacks. */
const MEASURE = `(root) => {
  const v = grok.shell.tv.viewers.find((x) => x.root === root || root.contains(x.root) || x.root.contains(root));
  const cv = [...v.root.querySelectorAll('canvas')].sort((a, b) => b.clientWidth * b.clientHeight - a.clientWidth * a.clientHeight)[0];
  const box = cv.getBoundingClientRect();
  const STEP_X = 4, STEP_Y = 3, ROW_PX = 18;
  const hits = {};
  let name, listed = 0;
  for (let y = 1; y < box.height; y += STEP_Y) {
    for (let x = 1; x < box.width; x += STEP_X) {
      ui.tooltip.hide();
      cv.dispatchEvent(new MouseEvent('mousemove', {clientX: box.x + x, clientY: box.y + y, bubbles: true}));
      const tip = [...document.querySelectorAll('.d4-tooltip')].find((t) => getComputedStyle(t).display !== 'none' && t.querySelector('canvas'));
      if (!tip)
        continue;
      const rows = tip.querySelectorAll('canvas').length;
      const key = rows > 1 ? '${MORE}' : '${SHORTENED}';
      if (rows > 1)
        listed = rows;
      else if (!hits[key])
        name = tip.textContent.trim();
      else if (tip.textContent.trim() !== name || y > hits[key].y0 + ROW_PX - STEP_Y)
        continue;
      const h = hits[key] ??= {x0: x, y0: y, x1: x, y1: y};
      h.x0 = Math.min(h.x0, x); h.y0 = Math.min(h.y0, y); h.x1 = Math.max(h.x1, x); h.y1 = Math.max(h.y1, y);
    }
  }
  ui.tooltip.hide();
  const areas = {};
  for (const [k, h] of Object.entries(hits))
    areas[k] = {x: h.x0 - STEP_X / 2, y: h.y0 - STEP_Y / 2, width: h.x1 - h.x0 + STEP_X, height: h.y1 - h.y0 + STEP_Y};
  v.addStatusProvider('curves-legend', () => ({hitAreas: areas}));
  return {areas, name, listed, width: box.width, height: box.height};
}`;

const lastFound = new WeakMap<Page, LegendFound>();

async function measureLegend(page: Page, target: ElementRef): Promise<LegendFound> {
  await viewers.settle(page, target);
  const root = await locate(page, target);
  const found = await root.evaluate((el, src) => (new Function(`return (${src});`)())(el), MEASURE) as LegendFound;
  lastFound.set(page, found);
  return found;
}

export const hoverShortened = When('user hovers over a shortened name in the curves legend of {widget}', async (page: Page, target: ElementRef) => {
  const found = await measureLegend(page, target);
  expect(found.areas[SHORTENED], `a shortened name in the legend of ${target.phrase}`).toBeTruthy();
  await hoverArea(page, SHORTENED, target);
}, {tier: 'ui', description: 'finds the first legend row whose name the legend shortened, and moves the pointer onto it'});

export const hoverMore = When('user hovers over the "+N more" row of the curves legend of {widget}', async (page: Page, target: ElementRef) => {
  const found = await measureLegend(page, target);
  expect(found.areas[MORE], `the "+N more" row in the legend of ${target.phrase}`).toBeTruthy();
  await hoverArea(page, MORE, target);
}, {tier: 'ui', description: 'finds the row counting the curves the legend had no room for, and moves the pointer onto it'});

/* Just below the legend, where the curves leave the corner free: the tooltip closes as the pointer leaves
   the row, as it does for a person heading elsewhere. */
export const pointerOffLegend = When('user moves the pointer off the curves legend of {widget}', async (page: Page, target: ElementRef) => {
  const row = (await viewers.hitAreas(page, target))[MORE];
  expect(row, `the "+N more" row in the legend of ${target.phrase}, measured before`).toBeTruthy();
  await page.mouse.move(row.x + row.width / 2, row.y + row.height + 30);
  await expect(page.locator('.d4-tooltip').filter({visible: true})).toHaveCount(0);
}, {tier: 'ui', description: 'moves the pointer 30 px below the "+N more" row the last hover measured, and waits for the tooltip to close'});

const visibleTooltip =(page: Page) => page.locator('.d4-tooltip').filter({visible: true}).first();

export const tooltipWholeName = Then('the tooltip should show the whole of the shortened name', async (page: Page) => {
  const name = lastFound.get(page)?.name;
  expect(name, 'the name the legend shortened').toBeTruthy();
  await expect(visibleTooltip(page)).toHaveText(name!);
  expect(name!.endsWith('…'), 'the tooltip says the name without the ellipsis').toBe(false);
});

export const tooltipListsCurves = Then('the tooltip should list the curves of the legend, each in its own color', async (page: Page) => {
  const listed = lastFound.get(page)!.listed;
  expect(listed, 'curves the "+N more" row stands for').toBeGreaterThan(1);
  const tip = visibleTooltip(page);
  await expect(tip.locator('canvas')).toHaveCount(listed);
  const colors = await tip.locator('canvas + div').evaluateAll((els) => els.map((e) => (e as HTMLElement).style.color));
  expect(new Set(colors).size, `the colors of the listed names (${colors.join(', ')})`).toBeGreaterThan(Math.min(listed, 8) - 1);
});

/* The corner is the one the "+N more" row is in: the legend's last row, so the side it hangs from. */
export const legendCorner = Then('the curves legend of {widget} should be in the {word} {word} corner', async (page: Page, target: ElementRef, vertical: string, horizontal: string) => {
  const found = await measureLegend(page, target);
  const row = found.areas[MORE];
  expect(row, `the "+N more" row in the legend of ${target.phrase}`).toBeTruthy();
  const corner = `${row.y + row.height / 2 < found.height / 2 ? 'top' : 'bottom'} ${row.x + row.width / 2 < found.width / 2 ? 'left' : 'right'}`;
  expect(corner, `the corner of ${target.phrase} its legend is in`).toBe(`${vertical} ${horizontal}`);
});

export const legendDrawn = Then('the curves legend of {widget} should be drawn', async (page: Page, target: ElementRef) => {
  await expect.poll(async () => Object.keys((await measureLegend(page, target)).areas).length, {message: `the legend of ${target.phrase}`})
    .toBeGreaterThan(0);
});

export const legendNotDrawn = Then('the curves legend of {widget} should not be drawn', async (page: Page, target: ElementRef) => {
  await expect.poll(async () => Object.keys((await measureLegend(page, target)).areas).length, {message: `the legend of ${target.phrase}`})
    .toBe(0);
});

/* A grid cell draws a legend only from 325 x 275 px inside its 6 px margin, more than a grid gives a row
   by default; the header is kept one line high so the rows get the height a help GIF's window has. */
export const sizeCurveCells = Given('the grid of the current view has rows {int} px high and the {string} column {int} px wide', async (page: Page, height: number, column: string, width: number) => {
  silent(page);
  await page.evaluate(([h, name, w]) => {
    const grid = grok.shell.tv.grid;
    grid.props.rowHeight = h;
    grid.props.colHeaderHeight = 30;
    grid.columns.byName(name).width = w;
  }, [height, column, width] as [number, string, number]);
}, {tier: 'api', description: 'sets the row height of the current grid, a one-line header and the width of one of its columns; silent in a guide'});

/* The legend is drawn on the grid's canvas, so a cell is read by what it writes there: one repaint with
   fillText recorded, each text placed in the cell its point falls in. A cell draws a legend when it
   writes one of its curves' names, a name shortened with an ellipsis, or the "+N more" row. Rows of
   the column the grid paints - in view, even partly - are the ones read. */
const LEGEND_CELLS = `async (name) => {
  const grid = grok.shell.tv.grid;
  const canvas = grid.canvas;
  const drawn = [];
  const original = CanvasRenderingContext2D.prototype.fillText;
  CanvasRenderingContext2D.prototype.fillText = function (text, x, y, ...rest) {
    if (this.canvas === canvas) {
      const m = this.getTransform();
      const k = canvas.width / canvas.clientWidth;
      drawn.push({text: String(text), x: (m.a * x + m.c * y + m.e) / k, y: (m.b * x + m.d * y + m.f) / k});
    }
    return original.call(this, text, x, y, ...rest);
  };
  try {
    grid.invalidate();
    for (let i = 0; i < 3; i++)
      await new Promise((r) => requestAnimationFrame(r));
  }
  finally {
    CanvasRenderingContext2D.prototype.fillText = original;
  }
  const col = grid.columns.byName(name);
  const rows = [];
  for (let r = 0; r < grid.table.filter.trueCount; r++) {
    const cell = grid.cell(name, r);
    const b = cell.bounds;
    if (b.y + b.height <= 0 || b.y >= canvas.clientHeight || b.x + b.width <= 0 || b.x >= canvas.clientWidth)
      continue;
    let names = [];
    try {
      names = (JSON.parse(cell.cell.value || '{}').series ?? []).map((s) => s.name).filter((n) => n);
    }
    catch (_) { }
    const inside = drawn.filter((d) => d.x >= b.x && d.x <= b.x + b.width && d.y >= b.y && d.y <= b.y + b.height);
    const legend = inside.filter((d) => names.includes(d.text) || /^\\+\\d+ more$/.test(d.text) ||
      (d.text.endsWith('…') && names.some((n) => n.startsWith(d.text.slice(0, -1)))));
    rows.push({row: cell.tableRowIndex, legend: legend.map((d) => d.text)});
  }
  return {rows, current: grid.table.currentRowIdx, column: col ? name : null};
}`;

interface LegendCells {rows: {row: number; legend: string[]}[]; current: number; column: string | null}

async function legendCells(page: Page, column: string): Promise<LegendCells> {
  const found = await page.evaluate((args) => (new Function(`return (${args[0]});`)())(args[1]), [LEGEND_CELLS, column]) as LegendCells;
  expect(found.column, `the "${column}" column of the current grid`).toBe(column);
  expect(found.rows.length, `the rows of "${column}" in view`).toBeGreaterThan(1);
  return found;
}

const describe = (c: LegendCells) => c.rows.map((r) => `row ${r.row}${r.row === c.current ? ' (current)' : ''}: ${r.legend.length ? r.legend.join(', ') : 'no legend'}`).join('; ');

export const noCellLegend = Then('no curve cell in view of {string} column should draw a legend', async (page: Page, column: string) => {
  await expect.poll(async () => {
    const c = await legendCells(page, column);
    return c.rows.filter((r) => r.legend.length > 0).length === 0 ? 'none' : describe(c);
  }, {message: `the legends drawn in "${column}"`}).toBe('none');
});

export const everyCellLegend = Then('every curve cell in view of {string} column should draw a legend', async (page: Page, column: string) => {
  await expect.poll(async () => {
    const c = await legendCells(page, column);
    return c.rows.every((r) => r.legend.length > 0) ? 'all' : describe(c);
  }, {message: `the legends drawn in "${column}"`}).toBe('all');
});

export const onlyCurrentCellLegend = Then('only the current curve cell of {string} column should draw a legend', async (page: Page, column: string) => {
  await expect.poll(async () => {
    const c = await legendCells(page, column);
    const ok = c.rows.some((r) => r.row === c.current) && c.rows.every((r) => (r.legend.length > 0) === (r.row === c.current));
    return ok ? 'only current' : describe(c);
  }, {message: `the legends drawn in "${column}"`}).toBe('only current');
});

/* The Options pane is taller than a help GIF's window: a person turns the wheel over the panel until the
   input shows, either way, rather than the page jumping to it. */
export const wheelIntoView = When('user scrolls {element} into view', async (page: Page, target: ElementRef) => {
  const loc = (await locate(page, target)).first();
  const aim = await loc.evaluate((e) => {
    let list = e.parentElement;
    while (list && !(/auto|scroll/.test(getComputedStyle(list).overflowY) && list.scrollHeight > list.clientHeight))
      list = list.parentElement;
    if (!list)
      return null;
    const r = e.getBoundingClientRect();
    const l = list.getBoundingClientRect();
    const MARGIN = 24;
    const dy = r.top < l.top + MARGIN ? r.top - l.top - MARGIN : r.bottom > l.bottom - MARGIN ? r.bottom - l.bottom + MARGIN : 0;
    return {x: l.x + l.width / 2, y: l.y + l.height / 2, dy, top: l.top, bottom: l.bottom};
  });
  if (aim && Math.abs(aim.dy) > 1) {
    await page.mouse.move(aim.x, aim.y, {steps: 6});
    await page.mouse.wheel(0, aim.dy);
  }
  await expect.poll(async () => {
    const box = await loc.boundingBox();
    return !!box && (!aim || (box.y >= aim.top && box.y + box.height <= aim.bottom));
  }, {message: `${target.phrase} within its scrolled panel`}).toBe(true);
}, {tier: 'ui', description: 'the mouse wheel over the panel the element scrolls in until the whole of it shows, up or down'});
