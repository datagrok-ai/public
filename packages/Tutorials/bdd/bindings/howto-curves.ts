/* The Dose-response curves how-to clips (features/howto/16-curves): the Curves package's Assay Curves table, and a point
   of a curve drawn in a grid cell pointed at or clicked. The Curves renderer reports no hit areas for its points (what
   a renderer draws is the package's test), so a point is placed here the way the renderer's own hit test places it
   (Curves src/fit/fit-interaction.ts `pointViewport`, fit-layout.ts `layoutChart`, the statistics library's
   `getChartBounds`): the cell's chart options, its points' bounds padded by 5% (log axes by a factor of 1.05), the
   plot box inside the cell's margins. Everything else is library vocabulary. */
import type {Page} from '@playwright/test';
import {Given, Then, When, dataset} from '@datagrok-libraries/bdd';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {exactText, expect, guideHop, locate, pollMs, silent, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

dataset('assay-curves', {path: 'System:AppData/Curves/demo/Assay_Curves.csv', aliases: ['Assay Curves'],
  description: '43 dose-response curves (4PL, log concentration) in "Fitted Curve", one per compound, assay, target and run, with Molecule, IC50, Hill, Min, Max, AUC'});

type Offset = {x: number; y: number};

/** What the "+" next to a statistic in Context Panel > Fit does (Curves fit-grid-cell-handler.ts `addStatisticColumn`):
 * the calculated column `<curve column> 1 <statistic>` of the first series, put right after the curve column. */
export const statisticColumnAdded = Given('the {string} statistic of the {string} column is a column of the table',
  async (page: Page, statistic: string, column: string) => {
    silent(page);
    await page.evaluate(async ([stat, name]) => {
      const DG = (window as any).DG;
      const tv = grok.shell.tv;
      const df = tv.dataFrame;
      const func = DG.Func.find({package: 'Curves', name: 'curveStatistic'})[0];
      await func.prepare({table: df, curveColumn: df.col(name), propName: stat, seriesNumber: 0})
        .call(false, undefined, {processed: false});
      const added = `${name} 1 ${stat}`;
      const order = (names: string[]) => {
        const rest = names.filter((n) => n !== added);
        const at = rest.indexOf(name) + 1;
        return [...rest.slice(0, at), added, ...rest.slice(at)];
      };
      df.columns.setOrder(order(df.columns.names()));
      const gridNames: string[] = [];
      for (let i = 0; i < tv.grid.columns.length; i++) {
        const n = tv.grid.columns.byIndex(i)?.name;
        if (n)
          gridNames.push(n);
      }
      tv.grid.columns.setOrder(order(gridNames));
    }, [statistic, column] as [string, string]);
    await expect.poll(() => page.evaluate((n) => grok.shell.tv.dataFrame.col(n) != null, `${column} 1 ${statistic}`),
      {message: `the "${column} 1 ${statistic}" column`}).toBe(true);
    await viewers.settleAll(page);
  }, {tier: 'api', description: 'set-up a clip leaves out (another clip shows the + of the Fit pane): the statistic of the curves\' first series as a calculated column after the curve column'});

/** The starting state of a clip about one curve further down the table: that row at the top of the grid. */
export const gridScrolledToRow = Given('the grid of the current view starts at row {int}', async (page: Page, row: number) => {
  silent(page);
  await page.evaluate((r) => grok.shell.tv.grid.vertScroll.scrollTo(r - 1), row);
  await expect.poll(() => page.evaluate(() => Math.round(grok.shell.tv.grid.vertScroll.min) + 1),
    {message: 'the first row the grid shows'}).toBe(row);
  await viewers.settleAll(page);
}, {tier: 'api', description: 'set-up a clip leaves out: the grid scrolled so that the row (from 1) is the first it shows'});

/** Where point [point] (from 1, in the order the cell stores them) of the cell's first curve is drawn, relative to the
 * cell's top left corner, as the renderer's hit test places it. */
async function pointOffset(page: Page, target: ElementRef, area: string, point: number): Promise<Offset> {
  const cellOf = /^cell (\d+) of (.+)$/.exec(area);
  if (!cellOf)
    throw new Error(`"${area}" is no grid cell: a point is aimed at in "cell <row> of <column>"`);
  const found = await viewers.onViewer(page, target, (el, [row, column, index]) => {
    const grid = (window as any).__bdd.viewerOf(el);
    const gc = grid.cell(column, grid.tableRowToGrid(row - 1));
    const data = JSON.parse(gc.cell.value);
    const tag = gc.cell.column.getTag('.%fit') ?? gc.cell.column.getTag('.fit');
    const o = {...(data.chartOptions ?? {}), ...(tag ? JSON.parse(tag).chartOptions ?? {} : {})};
    const b = gc.bounds;
    const sb = b.inflate(b.width < 70 ? -2 : -6, b.height < 45 ? -2 : -6);
    const roomy = sb.width >= 180 && sb.height >= 140;
    const showTitle = sb.width >= 275 && sb.height >= 225 && !!o.title;
    const box = sb.cutLeft(roomy && o.yAxisName ? 45 : 38).cutBottom(roomy && o.xAxisName ? 30 : 15)
      .cutTop(showTitle ? 25 : 5).cutRight(18);
    const xs: number[] = [];
    const ys: number[] = [];
    for (const s of data.series)
      for (const p of s.points) { xs.push(p.x); ys.push(p.y); }
    let [minX, maxX, minY, maxY] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
    if (o.logX) { minX = minX > 0 ? minX / 1.05 : minX; maxX *= 1.05; } else { const d = (maxX - minX) * 0.05; minX -= d; maxX += d; }
    if (o.logY) { minY = minY > 0 ? minY / 1.05 : minY; maxY *= 1.05; } else { const d = (maxY - minY) * 0.05; minY -= d; maxY += d; }
    if (o.minX != null) minX = o.minX;
    if (o.maxX != null) maxX = o.maxX;
    if (o.minY != null) minY = o.minY;
    if (o.maxY != null) maxY = o.maxY;
    const p = data.series[0].points[index - 1];
    if (!p)
      return `the curve has ${data.series[0].points.length} points`;
    const tx = o.logX ? Math.log(p.x / minX) / Math.log(maxX / minX) : (p.x - minX) / (maxX - minX);
    const ty = o.logY ? Math.log(p.y / minY) / Math.log(maxY / minY) : (p.y - minY) / (maxY - minY);
    return {x: box.left + tx * box.width - b.left, y: box.top + (1 - ty) * box.height - b.top};
  }, [Number(cellOf[1]), cellOf[2], point] as [number, string, number]);
  if (typeof found === 'string')
    throw new Error(`${area} of ${target.phrase} has no point ${point}: ${found}`);
  return found as Offset;
}

export const hoverCurvePoint = When('user hovers over point {int} of the curve in the {string} area of {widget}',
  async (page: Page, point: number, area: string, target: ElementRef) => {
    const cell = await viewers.hitArea(page, target, area, true);
    const at = await pointOffset(page, target, area, point);
    await page.mouse.move(cell.x + at.x - 3, cell.y + at.y - 3);
    await page.mouse.move(cell.x + at.x, cell.y + at.y);
    await expect(page.locator('.d4-tooltip').filter({visible: true}), `the tooltip of point ${point}`).toHaveCount(1);
  }, {tier: 'ui', description: 'a point of the cell\'s first curve, counted from 1 in the order the cell stores them; done when its tooltip (its x and y) shows'});

export const clickCurvePoint = When('user clicks on point {int} of the curve in the {string} area of {widget}',
  async (page: Page, point: number, area: string, target: ElementRef) => {
    const cell = await viewers.hitArea(page, target, area, true);
    const at = await pointOffset(page, target, area, point);
    await page.mouse.move(cell.x + at.x - 3, cell.y + at.y - 3);
    await page.mouse.click(cell.x + at.x, cell.y + at.y);
    await viewers.settle(page, target);
  }, {tier: 'ui', description: 'a click on a point of the cell\'s first curve (from 1, in stored order): a curve that allows it marks the point an outlier, or takes it back'});

/** A click in the middle of a curve cell can land on a point and leave it out of the fit (Click To Toggle): this one
 * goes to the strip along the cell's left edge, where the y axis is named, which no point reaches. */
export const clickBesideCurve = When('user clicks on the {string} area of {widget} beside its curve',
  async (page: Page, area: string, target: ElementRef) => {
    const cell = await viewers.hitArea(page, target, area, true);
    const at = {x: cell.x + 14, y: cell.y + cell.height / 2};
    await page.mouse.move(at.x + 4, at.y - 4);
    await page.mouse.click(at.x, at.y);
    await viewers.settle(page, target);
  }, {tier: 'ui', description: 'a click inside the cell, 14 px from its left edge (the y axis name), so no point of the curve is toggled'});

export const doubleClickBesideCurve = When('user double-clicks on the {string} area of {widget} beside its curve',
  async (page: Page, area: string, target: ElementRef) => {
    const cell = await viewers.hitArea(page, target, area, true);
    const at = {x: cell.x + 14, y: cell.y + cell.height / 2};
    await page.mouse.move(at.x + 4, at.y - 4);
    await page.mouse.dblclick(at.x, at.y);
  }, {tier: 'ui', description: 'a double-click inside the cell, 14 px from its left edge, so neither click toggles a point of the curve'});

async function isOutlier(page: Page, row: number, column: string, point: number): Promise<boolean | string> {
  return page.evaluate(([r, c, i]) => {
    const col = grok.shell.tv?.dataFrame?.col(c);
    if (!col)
      return `the table has no "${c}" column`;
    const p = JSON.parse(col.get(r - 1)).series[0].points[i - 1];
    return p ? !!p.outlier : `row ${r} has no point ${i}`;
  }, [row, column, point] as [number, string, number]);
}

export const pointLeftOut = Then('point {int} of the curve in row {int} of the {string} column should be left out',
  async (page: Page, point: number, row: number, column: string) => {
    await expect.poll(() => isOutlier(page, row, column, point), {timeout: pollMs(5000),
      message: `point ${point} of row ${row} of "${column}" marked an outlier`}).toBe(true);
  }, {description: 'the point is stored as an outlier in the cell (`outlier: true`), so the fit leaves it out'});

export const noPointLeftOut = Then('no point of the curve in row {int} of the {string} column should be left out',
  async (page: Page, row: number, column: string) => {
    const outliers = await page.evaluate(([r, c]) => {
      const col = grok.shell.tv?.dataFrame?.col(c);
      if (!col)
        return `the table has no "${c}" column`;
      return JSON.parse(col.get(r - 1)).series.flatMap((s: any) => s.points).filter((p: any) => p.outlier).length;
    }, [row, column] as [number, string]);
    expect(outliers, `points of row ${row} of "${column}" marked outliers`).toBe(0);
  }, {description: 'a click on the curve cell (away from its points) marked no point an outlier'});

export const pointKept = Then('point {int} of the curve in row {int} of the {string} column should be in the fit',
  async (page: Page, point: number, row: number, column: string) => {
    await expect.poll(() => isOutlier(page, row, column, point), {timeout: pollMs(5000),
      message: `point ${point} of row ${row} of "${column}" not an outlier`}).toBe(false);
  }, {description: 'the point is not stored as an outlier in the cell'});

/* The Options pane's Statistics and Droplines are Dart multi-choice inputs: a column of bare checkboxes, each with a
   <label for> and no name, which no library noun reaches. The label is what a person clicks. */
async function choiceBox(page: Page, option: string, target: ElementRef) {
  const host = (await locate(page, target)).filter({visible: true}).first();
  await expect(host, target.phrase).toBeVisible({timeout: pollMs(10000)});
  const item = host.locator('.ui-input-multi-choice-checks > div').filter({has: page.locator('label', {hasText: exactText(option)})}).first();
  await expect(item, `"${option}" in ${target.phrase}`).toBeVisible({timeout: pollMs(5000)});
  return {label: item.locator('label').first(), box: item.locator('input[type="checkbox"]').first()};
}

export const tickChoice = When('user ticks {string} in {element}', async (page: Page, option: string, target: ElementRef) => {
  const {label, box} = await choiceBox(page, option, target);
  await label.scrollIntoViewIfNeeded();
  await guideHop(page, label);
  await label.click();
  await expect(box, `"${option}" in ${target.phrase} after the click`).toBeChecked({checked: true, timeout: pollMs(5000)});
}, {tier: 'ui', description: 'a click on the label of an unticked option of a Dart multi-choice input (Statistics, Droplines): the option is ticked'});

export const untickChoice = When('user unticks {string} in {element}', async (page: Page, option: string, target: ElementRef) => {
  const {label, box} = await choiceBox(page, option, target);
  await label.scrollIntoViewIfNeeded();
  await guideHop(page, label);
  await label.click();
  await expect(box, `"${option}" in ${target.phrase} after the click`).toBeChecked({checked: false, timeout: pollMs(5000)});
}, {tier: 'ui', description: 'a click on the label of a ticked option of a Dart multi-choice input: the option is unticked'});

export const choiceTicked = Then('{string} should be ticked in {element}', async (page: Page, option: string, target: ElementRef) => {
  const {box} = await choiceBox(page, option, target);
  await expect(box, `"${option}" in ${target.phrase}`).toBeChecked({timeout: pollMs(5000)});
}, {description: 'the option\'s checkbox in a Dart multi-choice input is checked'});
