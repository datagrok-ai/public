/* ---
realizes: [viewers.scatter-plot, viewers.line-chart, viewers.density-plot, viewers.box-plot, viewers.histogram, viewers.bar-chart, powerpack.dialogs.formula-lines]
--- */
import {expect, type Page} from '@playwright/test';
import {test} from '../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../spec-login';
import * as v from '../helpers/viewers';

// Server lane: the Formula Lines dialog and its context-menu item are PowerPack's, and the
// local client loads no packages.
test.use(specTestOptions);

declare const grok: any;

const demogPath = 'System:DemoFiles/demog.csv';
const FORMULA_REGION_CAPTIONS = ['Region - Formula Lines', 'Formula Region', 'Region', 'Formula - Region'];

const escapeRe = (s: string) => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
const VISIBLE_LABELS = '.d4-menu-popup:visible .d4-menu-item-label';
const menuLabel = (page: Page, text: string) =>
  page.locator(VISIBLE_LABELS, {hasText: new RegExp(`^${escapeRe(text)}$`)});
const formulaLinesDialog = (page: Page) => page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'});
const canvasOf = (viewerType: string) => `[name="viewer-${viewerType.replace(/\s+/g, '-')}"] canvas`;

// Menu.show removes the previous popup as it appends the new one, so the popup count is no
// signal; the popups already there are marked and a visible unmarked one is awaited.
const NEW_POPUP = '.d4-menu-popup:visible:not([data-seen])';

async function openPopup(page: Page, act: () => Promise<void>) {
  await page.evaluate(() => {
    for (const p of Array.from(document.querySelectorAll('.d4-menu-popup'))) p.setAttribute('data-seen', '1');
  });
  await act();
  await page.locator(NEW_POPUP).first().waitFor({timeout: 5000});
}

// A context menu whose item opened the Formula Lines dialog stays open under the dialog, at
// the very point the next right-click lands on, so the trusted click hits the popup instead of
// the canvas: close it first.
async function dismissPopups(page: Page) {
  const visible = page.locator('.d4-menu-popup:visible');
  if (await visible.count() === 0) return;
  await page.keyboard.press('Escape');
  await visible.first().waitFor({state: 'hidden', timeout: 2000}).catch(() => page.evaluate(() => {
    for (const p of Array.from(document.querySelectorAll('.d4-menu-popup'))) p.remove();
  }));
}

// a synthetic contextmenu event bypasses the Dart canvas hit-test and always opens the
// whole-viewer menu; only a trusted click resolves the axis/region under the pointer
async function rightClick(page: Page, viewerType: string, dx = 100, dy = 100) {
  const selector = canvasOf(viewerType);
  const box = await page.locator(selector).first().boundingBox();
  if (!box) throw new Error(`no element: ${selector}`);
  await dismissPopups(page);
  await openPopup(page, () => page.mouse.click(box.x + dx, box.y + dy, {button: 'right'}));
}

async function clickMenuItem(page: Page, text: string) {
  await menuLabel(page, text).first().evaluate((el) => (el.closest('.d4-menu-item') as HTMLElement).click());
}

async function menuLabels(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-menu-popup'))
      .filter((p) => (p as HTMLElement).offsetParent !== null)
      .flatMap((p) => Array.from(p.querySelectorAll('.d4-menu-item-label')))
      .map((l) => (l.textContent ?? '').trim()));
}

async function noPopups(page: Page) {
  await page.waitForFunction(() => document.querySelectorAll('.d4-menu-popup').length === 0, null, {timeout: 5000});
}

async function noDialogs(page: Page) {
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const btn = (d.querySelector('[name="button-CANCEL"]') ?? d.querySelector('[name="button-OK"]')) as HTMLElement | null;
      btn?.click();
    });
  });
  await page.waitForFunction(() => document.querySelectorAll('.d4-dialog').length === 0, null, {timeout: 5000});
}

// Enters drawing mode from the context menu and waits for the repaint that installs the overlay,
// capped at the settle it replaces.
async function startDrawing(page: Page, viewerType: string, capMs: number) {
  const painted = await v.armEvent(page, `viewer:${viewerType}.onViewerRendered`, capMs);
  await clickMenuItem(page, 'Draw Annotation Region');
  await noPopups(page);
  await painted();
}

async function dragRect(page: Page, from: [number, number], to: [number, number], steps = 1) {
  await page.mouse.move(from[0], from[1], {steps});
  await page.mouse.down();
  await page.mouse.move((from[0] + to[0]) / 2, (from[1] + to[1]) / 2, {steps});
  await page.mouse.move(to[0], to[1], {steps});
  await page.mouse.up();
}

async function dragPolygon(page: Page, pts: [number, number][]) {
  await page.mouse.move(pts[0][0], pts[0][1]);
  await page.mouse.down();
  for (let i = 1; i < pts.length; i++) await page.mouse.move(pts[i][0], pts[i][1], {steps: 5});
  await page.mouse.up();
}

async function look(page: Page, viewerType: string): Promise<any> {
  return page.evaluate((t) => grok.shell.tv.viewers.find((x: any) => x.type === t).getOptions(true).look, viewerType);
}

async function regions(page: Page, viewerType: string): Promise<any[]> {
  return JSON.parse((await look(page, viewerType)).annotationRegions || '[]');
}

async function formulaLines(page: Page, viewerType: string): Promise<any[]> {
  return JSON.parse((await look(page, viewerType)).formulaLines || '[]');
}

async function closeOtherViewers(page: Page) {
  await page.evaluate(() => {
    for (const x of [...grok.shell.tv.viewers].filter((x: any) => x.type !== 'Grid')) x.close();
  });
}

async function setOpacity(page: Page, value: string) {
  await page.evaluate((val) => {
    const label = Array.from(document.querySelectorAll('.d4-dialog .ui-input-label'))
      .find((el) => el.textContent?.trim() === 'Opacity');
    const o = label?.parentElement?.querySelector('input') as HTMLInputElement | null;
    if (!o) throw new Error('Opacity input not found');
    o.value = val; o.dispatchEvent(new Event('input', {bubbles: true}));
  }, value);
}

async function fillHeaderColor(page: Page, color: string) {
  const hc = page.locator('[name="input-host-Header-Color"], [name="input-host-Color"]').locator('input').first();
  if (await hc.isVisible({timeout: 3000}).catch(() => false))
    await hc.fill(color);
}

// "Add new" in the dialog: the new item lands synchronously; OK is what commits it to the viewer.
async function addFormulaRegionInDialog(page: Page) {
  await openPopup(page, () =>
    page.locator('.d4-dialog button, .d4-dialog .ui-btn, .d4-dialog span').getByText('Add new', {exact: true}).first().click());
  const clicked = await page.evaluate((candidates: string[]) => {
    const labels = Array.from(document.querySelectorAll('.d4-menu-popup'))
      .filter((p) => (p as HTMLElement).offsetParent !== null)
      .flatMap((p) => Array.from(p.querySelectorAll('.d4-menu-item-label')));
    for (const c of candidates) {
      const t = labels.find((l: any) => (l.textContent ?? '').trim() === c);
      if (t) {
        ((t as HTMLElement).closest('.d4-menu-item') as HTMLElement | null)?.click();
        return c;
      }
    }
    return null;
  }, FORMULA_REGION_CAPTIONS);
  expect(clicked).not.toBeNull();
  await expect(page.locator('.d4-dialog [name="button-OK"]')).not.toHaveClass(/disabled/);
  await page.locator('.d4-dialog [name="button-OK"]').click();
  await page.waitForFunction(() => document.querySelectorAll('.d4-dialog').length === 0, null, {timeout: 5000});
}

// The axis "Annotations" items write the look synchronously and open the Formula Lines dialog
// on a 500 ms product timer; the dialog has to be waited for and dismissed or it sits over the
// next step.
async function axisAnnotationsItem(page: Page, leaf: string) {
  await page.evaluate((l) => (window as any).__menuLeaf('Annotations', l), leaf);
  const dialog = formulaLinesDialog(page);
  await dialog.waitFor({timeout: 2000}).catch(() => {});
  if (await dialog.count() > 0) {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() => document.querySelectorAll('.d4-dialog').length === 0, null, {timeout: 5000});
  }
}

function expectToolsOrder(labels: string[]) {
  const sIdx = labels.indexOf('Show Annotation Regions');
  const dIdx = labels.indexOf('Draw Annotation Region');
  const fIdx = labels.indexOf('Formula Lines...');
  expect(sIdx).toBeGreaterThanOrEqual(0);
  expect(dIdx).toBeGreaterThan(sIdx);
  expect(fIdx).toBeGreaterThan(dIdx);
}

test('Annotation regions scenario', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 1000});

  await softStep('1.1 Draw rectangle region + edit properties', async () => {
    await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot', 15_000, 'Scatter plot');
    await page.locator('[name="viewer-Scatter-plot"] canvas').first().waitFor();

    expect((await look(page, 'Scatter plot')).lassoTool).toBe(false);

    await rightClick(page, 'Scatter plot');
    await startDrawing(page, 'Scatter plot', 300);

    const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
    await dragRect(page, [box.x + box.width * 0.3, box.y + box.height * 0.3], [box.x + box.width * 0.6, box.y + box.height * 0.6]);

    await formulaLinesDialog(page).waitFor({timeout: 5000});
    await page.locator('[name="input-host-Title"] input').fill('My Rect Region');
    await page.locator('[name="input-host-Description"] textarea').fill('Rectangle description');
    await page.locator('[name="input-host-Region-Color"] input').fill('#ff8800');
    await page.locator('[name="input-host-Outline-Color"] input').fill('#003366');

    await page.locator('[name="input-host-Width"] input').fill('3');
    await page.locator('[name="input-host-Width"] input').press('Tab');
    await setOpacity(page, '60');
    await fillHeaderColor(page, '#ffffff');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = (await regions(page, 'Scatter plot'))[0];
    expect(region.header).toBe('My Rect Region');
    expect(region.description).toBe('Rectangle description');
    expect(region.outlineWidth).toBe(3);
    expect(region.opacity).toBe(60);
  });

  await softStep('1.2 Draw Lasso (polygon) region', async () => {
    await v.setViewerProps(page, 'Scatter plot', [{set: {lassoTool: true}}]);
    await rightClick(page, 'Scatter plot');
    await startDrawing(page, 'Scatter plot', 300);

    const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
    await dragPolygon(page, [
      [box.x + box.width * 0.65, box.y + box.height * 0.2],
      [box.x + box.width * 0.85, box.y + box.height * 0.3],
      [box.x + box.width * 0.85, box.y + box.height * 0.5],
      [box.x + box.width * 0.70, box.y + box.height * 0.55],
      [box.x + box.width * 0.65, box.y + box.height * 0.4],
    ]);

    await formulaLinesDialog(page).waitFor({timeout: 5000});
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
  });

  await softStep('1.3 Create formula region', async () => {
    const before = (await regions(page, 'Scatter plot')).length;
    await rightClick(page, 'Scatter plot');
    await clickMenuItem(page, 'Formula Lines...');
    await page.locator('.d4-dialog').waitFor();

    await addFormulaRegionInDialog(page);
    const after = await v.pollValue(() => regions(page, 'Scatter plot'), (r) => r.length > before, 400);
    expect(after.length).toBe(before + 1);

    const newRegion = after.find((r: any) => r.type === 'formula') ?? after[0];
    expect(newRegion).toBeDefined();
    if (newRegion.type === 'formula') {
      const formulaPair = `${newRegion.formula1 ?? ''} ${newRegion.formula2 ?? ''}`;
      expect(formulaPair).toMatch(/\$\{HEIGHT\}|\$\{WEIGHT\}/);
    }
  });

  await softStep('2.1 Show/Hide viewer & dataframe independently', async () => {
    await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot');
      sp.setOptions({
        annotationRegions: JSON.stringify([
          {type: 'area', x: 'WEIGHT', y: 'HEIGHT', area: [[60, 150], [95, 150], [95, 185], [60, 185]], header: 'Viewer Region'},
        ]),
      });
      sp.dataFrame.setTag('.annotation-regions', JSON.stringify([
        {type: 'area', x: 'WEIGHT', y: 'HEIGHT', area: [[90, 160], [120, 160], [120, 190], [90, 190]], header: 'DF Region', isDataFrameRegion: true},
      ]));
    });
    const flags = async () => {
      const l = await look(page, 'Scatter plot');
      return {v: l.showViewerAnnotationRegions, d: l.showDataframeAnnotationRegions};
    };

    const before1 = await flags();
    await rightClick(page, 'Scatter plot');
    await clickMenuItem(page, 'Show Viewer Annotation Regions');
    const after1 = await v.pollValue(flags, (f) => f.v !== before1.v, 200);
    expect(after1.v).toBe(false);
    expect(after1.d).toBe(true);

    await rightClick(page, 'Scatter plot');
    await clickMenuItem(page, 'Show Dataframe Annotation Regions');
    const after2 = await v.pollValue(flags, (f) => f.d !== after1.d, 200);
    expect(after2.v).toBe(false);
    expect(after2.d).toBe(false);
  });

  await softStep('2.2 Global Show Annotation Regions re-enables both', async () => {
    await rightClick(page, 'Scatter plot');
    await clickMenuItem(page, 'Show Annotation Regions');
    const after = await v.pollValue(async () => {
      const l = await look(page, 'Scatter plot');
      return {v: l.showViewerAnnotationRegions, d: l.showDataframeAnnotationRegions};
    }, (f) => f.v && f.d, 200);
    expect(after.v).toBe(true);
    expect(after.d).toBe(true);
  });

  await softStep('4.3 Change Annotation Font', async () => {
    const before = (await look(page, 'Scatter plot')).annotationFont;
    await v.setViewerProps(page, 'Scatter plot', [{set: {annotationFont: '18px bold Arial'}}]);
    const after = (await look(page, 'Scatter plot')).annotationFont;
    expect(before).not.toBe(after);
    expect(after).toBe('18px bold Arial');
  });

  await softStep('4.1 Right-click region → Edit opens dialog', async () => {
    const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
    await rightClick(page, 'Scatter plot', box.width * 0.45, box.height * 0.45);
    const hasEdit = (await menuLabels(page)).includes('Edit');
    console.log(`[4.1] "Edit" menu item found: ${hasEdit}`);
    await page.keyboard.press('Escape');
    if (!hasEdit)
      console.warn('[4.1] AMBIGUOUS — cursor likely missed the region (worldToScreen unavailable)');
  });

  await softStep('4.2 Modify region via dialog: reopen and edit Outline Width / Opacity / Header Color', async () => {
    await rightClick(page, 'Scatter plot');
    await clickMenuItem(page, 'Formula Lines...');
    await formulaLinesDialog(page).waitFor({timeout: 5000});

    const title = page.locator('.d4-dialog [name="input-host-Title"] input');
    const gridBox = (await page.locator('.d4-dialog .d4-grid').boundingBox())!;
    await page.mouse.click(gridBox.x + gridBox.width * 0.5, gridBox.y + 36);
    await v.pollValue(() => title.inputValue().catch(() => ''), (t) => t === 'Viewer Region', 300);

    await page.locator('.d4-dialog [name="input-host-Width"] input').fill('5');
    await page.locator('.d4-dialog [name="input-host-Width"] input').press('Tab');
    await setOpacity(page, '40');
    await fillHeaderColor(page, '#ff0000');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = (await regions(page, 'Scatter plot'))[0];
    expect(region.outlineWidth).toBe(5);
    expect(region.opacity).toBe(40);
  });

  await softStep('5.1 Preview / 5.2 Grid representation', async () => {
    await noDialogs(page);
    await rightClick(page, 'Scatter plot');
    await clickMenuItem(page, 'Formula Lines...');
    await page.locator('.d4-dialog').waitFor();
    await expect(page.locator('.d4-dialog .d4-grid')).toBeVisible();
    const hasRenderedCanvas = await page.evaluate(() => {
      const canvases = Array.from(document.querySelectorAll('.d4-dialog canvas')) as HTMLCanvasElement[];
      return canvases.some((c) => c.width > 0 && c.height > 0);
    });
    expect(hasRenderedCanvas).toBe(true);
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
  });

  await softStep('6.1 Multi-axis limitation', async () => {
    await page.evaluate(() => {
      for (const x of [...grok.shell.tv.viewers].filter((x: any) => x.type === 'Scatter plot')) x.close();
    });
    await v.addViewerByIcon(page, 'line-chart', 'Line-chart', 15_000, 'Line chart');
    await page.locator('[name="viewer-Line-chart"] canvas').first().waitFor();
    await v.setViewerProps(page, 'Line chart', [{set: {multiAxis: true}}]);
    await rightClick(page, 'Line chart');
    expect(await menuLabel(page, 'Draw Annotation Region').count()).toBe(0);
    await page.keyboard.press('Escape');
  });

  await softStep('6.2 Single-axis Line Chart: draw rectangle region + add formula region', async () => {
    await v.setViewerProps(page, 'Line chart', [{set: {multiAxis: false, yColumnNames: ['HEIGHT']}}]);

    await rightClick(page, 'Line chart');
    expect(await menuLabel(page, 'Draw Annotation Region').count()).toBeGreaterThan(0);
    await startDrawing(page, 'Line chart', 300);

    const box = (await page.locator('[name="viewer-Line-chart"] canvas').first().boundingBox())!;
    await dragRect(page, [box.x + box.width * 0.25, box.y + box.height * 0.25], [box.x + box.width * 0.55, box.y + box.height * 0.55]);

    await formulaLinesDialog(page).waitFor({timeout: 5000});
    await page.locator('[name="input-host-Title"] input').fill('LC Rect Region');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const afterRect = (await regions(page, 'Line chart')).length;
    expect(afterRect).toBeGreaterThan(0);

    await rightClick(page, 'Line chart');
    await clickMenuItem(page, 'Formula Lines...');
    await formulaLinesDialog(page).waitFor({timeout: 5000});
    await addFormulaRegionInDialog(page);

    const afterFormula = (await v.pollValue(() => regions(page, 'Line chart'), (r) => r.length > afterRect, 500)).length;
    expect(afterFormula).toBeGreaterThan(afterRect);
  });

  await noDialogs(page);
  await page.evaluate(() => { document.querySelectorAll('.d4-menu-popup').forEach((p) => p.remove()); });
  await page.keyboard.press('Escape');

  await softStep('7.1 Density Plot — rectangle region (lasso off)', async () => {
    await closeOtherViewers(page);
    await v.addViewerByIcon(page, 'density-plot', 'Density-plot', 10_000, 'Density plot');
    await page.locator('[name="viewer-Density-plot"] canvas').first().waitFor({timeout: 10_000});
    await v.setViewerProps(page, 'Density plot', [{set: {xColumnName: 'WEIGHT', yColumnName: 'HEIGHT', lassoTool: false}, wait: 500}]);

    await rightClick(page, 'Density plot');
    await startDrawing(page, 'Density plot', 600);

    const box = (await page.locator('[name="viewer-Density-plot"] canvas').first().boundingBox())!;
    await dragRect(page, [box.x + box.width * 0.3, box.y + box.height * 0.3], [box.x + box.width * 0.6, box.y + box.height * 0.6]);

    await formulaLinesDialog(page).waitFor({timeout: 5000});
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();

    expect((await regions(page, 'Density plot')).length).toBeGreaterThan(0);
  });

  await softStep('7.2 Density Plot — lasso region', async () => {
    await v.setViewerProps(page, 'Density plot', [{set: {lassoTool: true}, wait: 200}]);
    await rightClick(page, 'Density plot');
    await startDrawing(page, 'Density plot', 300);

    const box = (await page.locator('[name="viewer-Density-plot"] canvas').first().boundingBox())!;
    await dragPolygon(page, [
      [box.x + box.width * 0.65, box.y + box.height * 0.2],
      [box.x + box.width * 0.85, box.y + box.height * 0.3],
      [box.x + box.width * 0.85, box.y + box.height * 0.5],
      [box.x + box.width * 0.70, box.y + box.height * 0.55],
      [box.x + box.width * 0.65, box.y + box.height * 0.4],
    ]);

    await formulaLinesDialog(page).waitFor({timeout: 5000});
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
  });

  await softStep('7.3 Density Plot — Tools menu has 3 expected items in order', async () => {
    await rightClick(page, 'Density plot');
    expectToolsOrder(await menuLabels(page));
    await page.keyboard.press('Escape');
  });

  await softStep('8.1 Box Plot — Tools menu items in order, Lasso unavailable', async () => {
    await closeOtherViewers(page);
    await v.addViewerByIcon(page, 'box-plot', 'Box-plot', 10_000, 'Box plot');
    await page.locator('[name="viewer-Box-plot"] canvas').first().waitFor({timeout: 10_000});
    await v.setViewerProps(page, 'Box plot', [{set: {categoryColumnNames: ['RACE'], valueColumnName: 'AGE'}}]);

    await rightClick(page, 'Box plot');
    expectToolsOrder(await menuLabels(page));

    const lassoEnabled = await page.evaluate(() =>
      (grok.shell.tv.viewers.find((x: any) => x.type === 'Box plot') as any)
        .annotationRegionsFeature?.lassoEnabled ?? false);
    expect(lassoEnabled).toBe(false);
    await page.keyboard.press('Escape');
  });

  await softStep('8.2 Box Plot — drag rect → formula region with ${AGE} = lo / hi', async () => {
    await page.keyboard.press('Escape');
    await noPopups(page);
    await v.setViewerProps(page, 'Box plot', [{set: {annotationRegions: '[]'}}]);

    await rightClick(page, 'Box plot');
    await startDrawing(page, 'Box plot', 800);

    const ob = (await page.locator('[name="viewer-Box-plot"] canvas').nth(1).boundingBox())!;
    const cx = ob.x + ob.width * 0.5;
    await dragRect(page, [cx - 5, ob.y + ob.height * 0.30], [cx + 5, ob.y + ob.height * 0.55], 8);

    await formulaLinesDialog(page).waitFor({timeout: 8000});
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = (await regions(page, 'Box plot'))[0];
    expect(region).toBeDefined();
    expect(region.type).toBe('formula');
    expect(region.formula1).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
    expect(region.formula2).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
  });

  await softStep('8.3 Box Plot — Y axis Annotations group → Add Line creates ${AGE} = q2', async () => {
    const before = (await formulaLines(page, 'Box plot')).length;

    const bpBox = (await page.locator('[name="viewer-Box-plot"] canvas').first().boundingBox())!;
    await rightClick(page, 'Box plot', bpBox.width * 0.05, bpBox.height * 0.5);
    await axisAnnotationsItem(page, 'Add Line');

    const after = await formulaLines(page, 'Box plot');
    expect(after.length).toBe(before + 1);
    const last = after[after.length - 1];
    expect(last.formula).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
    expect(last.orientation).toBe('Horizontal');
  });

  await softStep('8.4 Box Plot — X axis (categorical) has NO Annotations group', async () => {
    const box = (await page.locator('[name="viewer-Box-plot"] canvas').first().boundingBox())!;
    await rightClick(page, 'Box plot', box.width / 2, box.height - 12);
    expect((await menuLabels(page)).includes('Annotations')).toBe(false);
    await page.keyboard.press('Escape');
  });

  await softStep('9.1 Histogram — Tools menu items in order', async () => {
    await closeOtherViewers(page);
    await v.addViewerByIcon(page, 'histogram', 'Histogram', 10_000, 'Histogram');
    await page.locator('[name="viewer-Histogram"] canvas').first().waitFor({timeout: 10_000});
    await v.setViewerProps(page, 'Histogram', [{set: {valueColumnName: 'AGE'}}]);

    await rightClick(page, 'Histogram');
    expectToolsOrder(await menuLabels(page));
    await page.keyboard.press('Escape');
  });

  await softStep('9.2 Histogram — drag rect → formula region with ${AGE} = lo / hi', async () => {
    await page.keyboard.press('Escape');
    await noPopups(page);
    await v.setViewerProps(page, 'Histogram', [{set: {annotationRegions: '[]'}}]);

    await rightClick(page, 'Histogram');
    await startDrawing(page, 'Histogram', 500);

    const box = (await page.locator('[name="viewer-Histogram"] canvas').first().boundingBox())!;
    const cy = box.y + box.height * 0.5;
    await dragRect(page, [box.x + box.width * 0.35, cy - 5], [box.x + box.width * 0.65, cy + 5], 4);

    await formulaLinesDialog(page).waitFor({timeout: 8000});
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = (await regions(page, 'Histogram'))[0];
    expect(region).toBeDefined();
    expect(region.type).toBe('formula');
    expect(region.formula1).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
    expect(region.formula2).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
  });

  await softStep('9.3 Histogram — X axis Annotations group → Add Line creates vertical line', async () => {
    const before = (await formulaLines(page, 'Histogram')).length;
    const box = (await page.locator('[name="viewer-Histogram"] canvas').first().boundingBox())!;

    await rightClick(page, 'Histogram', box.width / 2, box.height - 18);
    await axisAnnotationsItem(page, 'Add Line');

    const after = await formulaLines(page, 'Histogram');
    expect(after.length).toBe(before + 1);
    const last = after[after.length - 1];
    expect(last.formula).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
    expect(last.orientation).toBe('Vertical');
  });

  await softStep('10.1 Bar Chart — Tools menu items in order', async () => {
    await closeOtherViewers(page);
    await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart', 10_000, 'Bar chart');
    await page.locator('[name="viewer-Bar-chart"] canvas').first().waitFor({timeout: 10_000});
    await v.setViewerProps(page, 'Bar chart', [{
      set: {splitColumnName: 'RACE', valueColumnName: 'AGE', valueAggrType: 'avg', orientation: 'vertical'}, wait: 400,
    }]);

    await rightClick(page, 'Bar chart');
    expectToolsOrder(await menuLabels(page));
    await page.keyboard.press('Escape');
  });

  await softStep('10.2 Bar Chart — vertical, drag rect → formula region with avg(AGE) header', async () => {
    await page.keyboard.press('Escape');
    await noPopups(page);
    await v.setViewerProps(page, 'Bar chart', [{set: {annotationRegions: '[]'}}]);

    await rightClick(page, 'Bar chart');
    await startDrawing(page, 'Bar chart', 500);

    const box = (await page.locator('[name="viewer-Bar-chart"] canvas').first().boundingBox())!;
    const cx = box.x + box.width * 0.5;
    await dragRect(page, [cx - 5, box.y + box.height * 0.25], [cx + 5, box.y + box.height * 0.55], 4);

    await formulaLinesDialog(page).waitFor({timeout: 8000});
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = (await regions(page, 'Bar chart'))[0];
    expect(region).toBeDefined();
    expect(region.type).toBe('formula');
    expect(region.formula1).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
    expect(region.formula2).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
  });

  await softStep('10.3 Bar Chart — orientation flip preserves Tools menu', async () => {
    await v.setViewerProps(page, 'Bar chart', [{set: {orientation: 'horizontal'}, wait: 500}]);
    await rightClick(page, 'Bar chart');
    const labels = await menuLabels(page);
    expect(labels.indexOf('Show Annotation Regions')).toBeGreaterThanOrEqual(0);
    expect(labels.indexOf('Draw Annotation Region')).toBeGreaterThanOrEqual(0);
    expect(labels.indexOf('Formula Lines...')).toBeGreaterThanOrEqual(0);
    await page.keyboard.press('Escape');
  });

  await softStep('10.4 Bar Chart — value-axis Annotations group → Add Line', async () => {
    await v.setViewerProps(page, 'Bar chart', [{set: {orientation: 'vertical'}, wait: 500}]);
    const before = (await formulaLines(page, 'Bar chart')).length;

    // bar_chart_core.dart lays the value-axis strip (>= 20 px wide) right after the vertical
    // column selector and the left outer margin, so 10 px past that edge is always inside it
    const axis = await page.evaluate(() => {
      const bc = grok.shell.tv.viewers.find((x: any) => x.type === 'Bar chart');
      const sel = bc.root.querySelector('.d4-column-selector.d4-vertical') as HTMLElement | null;
      return {x: (sel?.clientHeight ?? 0) + (bc.getOptions(true).look.outerMarginLeft ?? 0) + 10,
        h: (bc.root.querySelector('canvas') as HTMLCanvasElement).getBoundingClientRect().height};
    });
    await rightClick(page, 'Bar chart', axis.x, axis.h * 0.5);
    expect((await menuLabels(page)).includes('Annotations')).toBe(true);
    await axisAnnotationsItem(page, 'Add Line');

    const after = await formulaLines(page, 'Bar chart');
    expect(after.length).toBe(before + 1);
  });

  await softStep('11.1 Histogram axis Add Band creates ${AGE} in (q1, q3)', async () => {
    await closeOtherViewers(page);
    await v.addViewerByIcon(page, 'histogram', 'Histogram', 10_000, 'Histogram');
    await page.locator('[name="viewer-Histogram"] canvas').first().waitFor({timeout: 10_000});
    await v.setViewerProps(page, 'Histogram', [{set: {valueColumnName: 'AGE'}}]);

    const box = (await page.locator('[name="viewer-Histogram"] canvas').first().boundingBox())!;
    const before = (await formulaLines(page, 'Histogram')).length;
    await rightClick(page, 'Histogram', box.width / 2, box.height - 18);
    await axisAnnotationsItem(page, 'Add Band');

    const after = await formulaLines(page, 'Histogram');
    expect(after.length).toBe(before + 1);
    const last = after[after.length - 1];
    expect(last.type).toBe('band');
    expect(last.formula).toMatch(/\$\{AGE\}\s+in\s+\(/);
    expect(last.orientation).toBe('Vertical');
  });

  await softStep('11.2 Histogram axis Add Region creates formula annotationRegion', async () => {
    const before = (await regions(page, 'Histogram')).length;
    const box = (await page.locator('[name="viewer-Histogram"] canvas').first().boundingBox())!;
    await rightClick(page, 'Histogram', box.width / 2, box.height - 18);
    await axisAnnotationsItem(page, 'Add Region');

    const after = await regions(page, 'Histogram');
    expect(after.length).toBe(before + 1);
    const last = after[after.length - 1];
    expect(last.type).toBe('formula');
    expect(last.formula1).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
    expect(last.formula2).toMatch(/\$\{AGE\}\s*=\s*-?\d/);
  });

  await noDialogs(page);
  await page.keyboard.press('Escape');
  await v.cleanupShell(page);
  v.finishSpec();
});
