/* ---
realizes: []
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const VIEWER_NAME = '3d-scatter-plot';
const VIEWER = `[name="viewer-${VIEWER_NAME}"]`;
const VIEWER_TYPE = '3d scatter plot';
const datasetPath = 'System:DemoFiles/demog.csv';
const HOME_CAMERA = [0, 0, 4];

// The plot is a three.js WebGL canvas without preserveDrawingBuffer, so its pixels are readable
// only in the task that painted them: the signature forces a paint and reads the buffer in one
// evaluate. An element screenshot cost ~1s and every step paid it at least twice.
async function installSignature(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    const proto = w.JsScatterPlot3dCore.prototype;
    if (!proto.__sp3dHooked) {
      const render = proto.render;
      proto.render = function(this: any) { w.__sp3d = this; return render.apply(this, arguments as any); };
      proto.__sp3dHooked = true;
    }
    // The whole drawing buffer, not a sample of it: sampling three full-width bands hashed the
    // axes shown and the axes hidden identically, and "Show Axes repaints" stopped detecting.
    w.__sp3dSig = () => {
      const plot = w.__sp3d;
      if (!plot) return null;
      plot.render();
      const gl = plot.renderer.getContext();
      const cv = plot.renderer.domElement;
      const buf = new Uint8Array(cv.width * cv.height * 4);
      gl.readPixels(0, 0, cv.width, cv.height, gl.RGBA, gl.UNSIGNED_BYTE, buf);
      let h = 0;
      for (let i = 0; i < buf.length; i += 16) h = (h * 31 + buf[i] + buf[i + 1] * 7 + buf[i + 2] * 13) % 2147483647;
      return h;
    };
    // The poll runs in the page: every read of the signature was a round trip queued behind the
    // scene's own render loop, and the read itself is a GPU sync worth spacing out.
    w.__sp3dRepaint = async (before: number, capMs: number) => {
      const deadline = Date.now() + capMs;
      let h = w.__sp3dSig();
      while (h === before && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 75));
        h = w.__sp3dSig();
      }
      return h;
    };
    w.__sp3dCamera = () => {
      const p = w.__sp3d?.camera?.position;
      return p ? [p.x, p.y, p.z].map((x: number) => Math.round(x * 1000) / 1000) : null;
    };
    const canvas = document.querySelector('[name="viewer-3d-scatter-plot"] canvas');
    canvas?.dispatchEvent(new MouseEvent('mousemove', {bubbles: true}));
  });
  await v.pollValue(() => page.evaluate(() => !!(window as any).__sp3d), (ok) => ok, 5000, 50);
}

const signature = (page: Page) => page.evaluate(() => (window as any).__sp3dSig() as number);
const camera = (page: Page) => page.evaluate(() => (window as any).__sp3dCamera() as number[]);
const repaints = (page: Page, before: number) =>
  page.evaluate((b) => (window as any).__sp3dRepaint(b, 3000) as Promise<number>, before);
const shownValue = (page: Page, prop: string) => v.propertyGridValue(page, prop);
// The probe read, the category expand and the settle are one in-page call: every read here was
// a Playwright actionability round trip queued behind the scene's render loop.
const rowVisible = (page: Page, probe: string) => page.evaluate((p) => {
  const row = document.querySelector(`.property-grid tr[name="prop-${p}"]`);
  return !!row && row.getClientRects().length > 0;
}, probe);

async function category(page: Page, cat: string, probe: string): Promise<void> {
  for (let attempt = 0; attempt < 3; attempt++) {
    const state = await page.evaluate(async ({c, p}) => {
      const w = window as any;
      const row = () => document.querySelector(`.property-grid tr[name="prop-${p}"]`);
      const shown = () => !!row() && row()!.getClientRects().length > 0;
      if (shown()) return 'shown';
      const header = document.querySelector(`[name="prop-category-${c}"]`) as HTMLElement | null;
      if (!header) return 'no-header';
      header.click();
      await w.__poll(shown, (ok: boolean) => ok, 1500, 25);
      return shown() ? 'shown' : 'hidden';
    }, {c: cat, p: probe});
    if (state === 'shown') break;
    if (state === 'no-header') {
      await v.clickViewerTitlebarIcon(page, VIEWER_NAME, 'icon-font-icon-settings').catch(() => {});
      await page.locator(`[name="prop-category-${cat}"]`).first().waitFor({timeout: 3000}).catch(() => {});
    }
  }
  expect(await rowVisible(page, probe)).toBe(true);
}

// Restoring a column the step has already asserted on is plumbing, not coverage: a trusted
// selector drive is ~1.3s of typing and commit polling, and the four picks the assertions read
// still go through the on-viewer selector.
async function setColumns(page: Page, props: Record<string, string>): Promise<void> {
  await page.evaluate(({vt, p}) => {
    const norm = (s: string) => s.replace(/[\s-]+/g, ' ').toLowerCase();
    const view = Array.from((window as any).grok.shell.tv.viewers)
      .find((x: any) => norm(x.type) === norm(vt)) as any;
    for (const k of Object.keys(p)) view.props[k] = p[k];
  }, {vt: VIEWER_TYPE, p: props});
}

async function selectorText(page: Page, role: string): Promise<string> {
  return (await page.locator(`${VIEWER} [name="div-column-combobox-${role}"]`).first().innerText())
    .replace(/\s+/g, ' ').trim();
}

async function plotCentre(page: Page): Promise<{x: number; y: number; box: any}> {
  const box = (await page.locator(VIEWER).boundingBox())!;
  return {x: box.x + box.width / 2, y: box.y + box.height / 2, box};
}

// The point to hover on the bar chart is read from its pixels: the bars are the only saturated
// colour on that canvas, so their centroid lands on a bar whatever the split column draws.
const barCentroid = (page: Page) => page.evaluate(() => {
  const cv = document.querySelector('[name="viewer-Bar-chart"] canvas') as HTMLCanvasElement | null;
  const ctx = cv?.getContext('2d');
  if (!cv || !ctx) return null;
  const d = ctx.getImageData(0, 0, cv.width, cv.height).data;
  let sx = 0, sy = 0, n = 0;
  for (let i = 0; i < d.length; i += 4) {
    const r = d[i], g = d[i + 1], b = d[i + 2];
    if (Math.max(r, g, b) - Math.min(r, g, b) < 60) continue;
    const px = (i / 4) % cv.width;
    sx += px; sy += (i / 4 - px) / cv.width; n++;
  }
  if (n === 0) return null;
  const rect = cv.getBoundingClientRect();
  return {x: rect.x + (sx / n) * rect.width / cv.width, y: rect.y + (sy / n) * rect.height / cv.height};
});

const clickMenuItem = (page: Page, label: string) => page.evaluate((text) => {
  const item = Array.from(document.querySelectorAll('.d4-menu-item'))
    .find((el) => el.querySelector('.d4-menu-item-label')?.textContent?.trim() === text) as HTMLElement | undefined;
  item?.click();
}, label);

const legendBox = (page: Page) => page.evaluate(() => {
  const el = document.querySelector('[name="viewer-3d-scatter-plot"] .d4-legend') as HTMLElement | null;
  if (!el) return null;
  const r = el.getBoundingClientRect();
  return {x: Math.round(r.x), y: Math.round(r.y), w: Math.round(r.width), h: Math.round(r.height), cls: el.className};
});

test('3D scatter plot', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await softStep('Add 3D scatter plot from the Viewers toolbox', async () => {
    await page.locator('[name="icon-3d-scatter-plot"]').first().click();
    await page.locator(VIEWER).first().waitFor({timeout: 30_000});
    await page.locator(`${VIEWER} canvas`).first().waitFor({timeout: 30_000});
    await installSignature(page);

    await expect.poll(() => selectorText(page, 'x'), {timeout: 30_000}).toBe('X: AGE');
    expect(await selectorText(page, 'y')).toBe('Y: HEIGHT');
    expect(await selectorText(page, 'z')).toBe('Z: WEIGHT');
  });

  await softStep('Reassign X and Z with the on-viewer selectors', async () => {
    const before = await signature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'x', columnName: 'WEIGHT', viewerType: VIEWER_TYPE, propName: 'xColumnName',
    });
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'z', columnName: 'AGE', viewerType: VIEWER_TYPE, propName: 'zColumnName',
    });
    expect(await selectorText(page, 'x')).toBe('X: WEIGHT');
    expect(await selectorText(page, 'z')).toBe('Z: AGE');
    expect(await repaints(page, before)).not.toBe(before);

    await setColumns(page, {xColumnName: 'AGE', zColumnName: 'WEIGHT'});
  });

  await softStep('Color by SEX shows a categorical legend', async () => {
    const before = await signature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'SEX', viewerType: VIEWER_TYPE, propName: 'colorColumnName',
    });

    await expect.poll(async () => (await v.readLegend(page, VIEWER_TYPE)).labels.sort(),
      {timeout: 10_000}).toEqual(['F', 'M']);
    expect((await v.readLegend(page, VIEWER_TYPE)).legendRendered).toBe(true);
    expect(await repaints(page, before)).not.toBe(before);
  });

  await softStep('Color by AGE switches the legend to a gradient', async () => {
    const before = await signature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'AGE', viewerType: VIEWER_TYPE, propName: 'colorColumnName',
    });

    await expect.poll(async () => (await v.readLegend(page, VIEWER_TYPE)).labels,
      {timeout: 10_000}).not.toEqual(['F', 'M']);
    expect(await repaints(page, before)).not.toBe(before);
  });

  await softStep('Marker type redraws the markers', async () => {
    await v.openViewerProperties(page, VIEWER_NAME);
    await category(page, 'marker', 'marker-type');

    const shapes: Record<string, number> = {};
    let current = await signature(page);
    for (const shape of ['box', 'sphere', 'cylinder']) {
      await v.selectPropertyGridChoice(page, 'marker-type', shape);
      expect(await shownValue(page, 'marker-type')).toBe(shape);
      current = await repaints(page, current);
      shapes[shape] = current;
    }
    expect(shapes['box']).not.toBe(shapes['sphere']);
    expect(shapes['cylinder']).not.toBe(shapes['sphere']);
  });

  await softStep('Marker opacity redraws the markers', async () => {
    await category(page, 'marker', 'marker-opacity');
    const before = await signature(page);
    await v.setPropertyGridValue(page, 'marker-opacity', '25');
    expect(await shownValue(page, 'marker-opacity')).toBe('25');
    const faded = await repaints(page, before);
    expect(faded).not.toBe(before);
    await v.setPropertyGridValue(page, 'marker-opacity', '100');
    await repaints(page, faded);
  });

  await softStep('Show Axes hides and restores the axes', async () => {
    await category(page, 'misc', 'show-axes');
    const before = await signature(page);
    expect(await v.togglePropertyGridCheckbox(page, 'show-axes')).toBe(false);
    const withoutAxes = await repaints(page, before);
    expect(withoutAxes).not.toBe(before);

    expect(await v.togglePropertyGridCheckbox(page, 'show-axes')).toBe(true);
    expect(await repaints(page, withoutAxes)).not.toBe(withoutAxes);
  });

  await softStep('X axis type switches to logarithmic', async () => {
    await category(page, 'x-axis', 'x-axis-type');
    const before = await signature(page);
    await v.selectPropertyGridChoice(page, 'x-axis-type', 'logarithmic');
    expect(await shownValue(page, 'x-axis-type')).toBe('logarithmic');
    const log = await repaints(page, before);
    expect(log).not.toBe(before);
    await v.selectPropertyGridChoice(page, 'x-axis-type', 'linear');
    expect(await shownValue(page, 'x-axis-type')).toBe('linear');
    await repaints(page, log);
  });

  await softStep('Drag rotates the scene and Reset View restores it', async () => {
    const {x, y, box} = await plotCentre(page);
    const before = await signature(page);
    const cameraBefore = await camera(page);

    await page.mouse.move(x, y);
    await page.mouse.down();
    // every mouse step costs ~0.9s of hit-testing in the plot, and the rotation is the same
    await page.mouse.move(box.x + box.width * 0.75, box.y + box.height * 0.3, {steps: 2});
    await page.mouse.up();
    const rotated = await repaints(page, before);
    expect(rotated).not.toBe(before);
    expect(await camera(page)).not.toEqual(cameraBefore);

    await page.mouse.click(x, y, {button: 'right'});
    await page.locator('.d4-menu-popup').first().waitFor({timeout: 5000});
    await clickMenuItem(page, 'Reset View');
    await expect(page.locator('.d4-menu-popup')).toHaveCount(0);
    expect(await repaints(page, rotated)).not.toBe(rotated);
    expect(await camera(page)).toEqual(HOME_CAMERA);
  });

  await softStep('Mouse wheel zooms the scene', async () => {
    const {x, y} = await plotCentre(page);
    await page.mouse.move(x, y);
    const before = await signature(page);
    const distanceBefore = (await camera(page))[2];

    await page.mouse.wheel(0, -600);
    const zoomedIn = await repaints(page, before);
    expect(zoomedIn).not.toBe(before);
    expect((await camera(page))[2]).toBeLessThan(distanceBefore);

    await page.mouse.wheel(0, 600);
    expect(await repaints(page, zoomedIn)).not.toBe(zoomedIn);
  });

  await softStep('Click makes a row current, Shift+click selects it', async () => {
    const {x, y, box} = await plotCentre(page);
    const startRow = await page.evaluate(() => (window as any).grok.shell.t.currentRowIdx);

    const currentRow = () =>
      page.evaluate(() => (window as any).grok.shell.t.currentRowIdx as number);
    let current = startRow;
    for (const [dx, dy] of [[0.5, 0.5], [0.45, 0.55], [0.55, 0.45], [0.5, 0.6]]) {
      await page.mouse.click(box.x + box.width * dx, box.y + box.height * dy);

      current = await v.pollValue(currentRow, (n) => n !== startRow && n >= 0, 1500);
      if (current !== startRow && current >= 0) break;
    }
    expect(current).toBeGreaterThanOrEqual(0);
    expect(current).not.toBe(startRow);

    const selected = () =>
      page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);
    const selectedBefore = await selected();
    await page.keyboard.down('Shift');
    await page.mouse.click(x, y);
    await page.keyboard.up('Shift');
    await expect.poll(selected, {timeout: 8000}).toBeGreaterThan(selectedBefore);

    await page.evaluate(() => (window as any).grok.shell.t.selection.setAll(false));
  });

  await softStep('Show Filtered Out Points repaints the filtered-away rows', async () => {
    const {filteredCount} = await v.applyCategoricalFilter(page, 'SEX', ['F']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeLessThan(total);

    await category(page, 'data', 'show-filtered-out-points');
    const filteredOnly = await signature(page);
    expect(await v.togglePropertyGridCheckbox(page, 'show-filtered-out-points')).toBe(true);
    expect(await repaints(page, filteredOnly)).not.toBe(filteredOnly);

    expect(await v.togglePropertyGridCheckbox(page, 'show-filtered-out-points')).toBe(false);
    await v.resetFilters(page);
  });

  await softStep('Hovering a bar chart bin highlights the matching 3D points', async () => {
    await page.locator('[name="icon-bar-chart"]').first().click();
    await page.locator('[name="viewer-Bar-chart"]').first().waitFor({timeout: 30_000});
    await expect.poll(async () =>
      (await v.countCanvasPixels(page, 'Bar chart')).total, {timeout: 30_000}).toBeGreaterThan(1000);

    const bar = (await page.locator('[name="viewer-Bar-chart"]').boundingBox())!;
    const target = (await barCentroid(page))!;
    expect(target).not.toBeNull();
    const idle = await signature(page);
    await page.mouse.move(target.x, target.y);
    expect(await repaints(page, idle)).not.toBe(idle);

    await page.mouse.move(bar.x + bar.width * 0.9, bar.y + bar.height * 0.05);
    await category(page, 'selection', 'show-mouse-over-row-group');
    expect(await v.togglePropertyGridCheckbox(page, 'show-mouse-over-row-group')).toBe(false);
    expect(await shownValue(page, 'show-mouse-over-row-group')).toBe('false');

    await category(page, 'selection', 'show-mouse-over-row-group');
    expect(await v.togglePropertyGridCheckbox(page, 'show-mouse-over-row-group')).toBe(true);
    await v.clickViewerTitlebarIcon(page, 'Bar-chart', 'Close');
    await expect(page.locator('[name="viewer-Bar-chart"]')).toHaveCount(0);
  });

  await softStep('Legend position moves the legend', async () => {
    // an earlier step left Color = AGE (numeric), which renders no legend element at all,
    // so colour by a categorical column first to have a legend to reposition
    await setColumns(page, {colorColumnName: 'SEX'});
    await expect.poll(async () => (await v.readLegend(page, VIEWER_TYPE)).legendRendered,
      {timeout: 10_000}).toBe(true);

    await category(page, 'legend', 'legend-position');
    await v.selectPropertyGridChoice(page, 'legend-visibility', 'Always');
    const before = await signature(page);
    const legendBefore = await legendBox(page);
    await v.selectPropertyGridChoice(page, 'legend-position', 'Left');
    expect(await shownValue(page, 'legend-position')).toBe('Left');
    expect(await repaints(page, before)).not.toBe(before);
    const legendLeft = await v.pollValue(() => legendBox(page), (b) => !!b && b.cls.includes('d4-legend-left'), 3000, 50);
    expect(legendLeft).not.toEqual(legendBefore);
    expect(legendLeft!.cls).toContain('d4-legend-left');
    await v.selectPropertyGridChoice(page, 'legend-position', 'Auto');
  });

  await page.evaluate(() => {
    const w = window as any;
    delete w.__sp3d; delete w.__sp3dSig; delete w.__sp3dCamera; delete w.__sp3dRepaint;
  });
  await v.cleanupShell(page);

  v.finishSpec();
});
