import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/geo/earthquakes.csv';

// Layers-panel rows are drawn on the mini-grid canvas, so a layer's visibility
// checkbox has no DOM node — it is hit by geometry: 24px rows, checkbox column
// at the panel's left edge. Row order matches the OpenLayers layer stack.
const LAYER_ROW = {'Bing sat': 0, 'BaseLayer': 1, 'Heatmap': 2, 'Markers GL': 3, 'Markers GL Selection': 4};

/** Visibility of every OpenLayers layer, keyed by the name shown in the layers panel. */
async function layerVisibility(page: Page): Promise<Record<string, boolean>> {
  return page.evaluate(() => {
    const mv = Array.from((window as any).grok.shell.tv.viewers).find((x: any) => x.type === 'Map') as any;
    const out: Record<string, boolean> = {};
    for (const l of mv.ol.olMap.getLayers().getArray())
      out[l.get('layerName') ?? l.get('name')] = l.getVisible();
    return out;
  });
}

async function openLayersPanel(page: Page): Promise<void> {
  await page.locator('[name="button-Map-layers"]').click();
  await page.locator('.panel-layers').waitFor({state: 'visible', timeout: 5000});
}

async function toggleLayer(page: Page, layer: keyof typeof LAYER_ROW): Promise<void> {
  const panel = await page.locator('.panel-layers').boundingBox();
  if (!panel) throw new Error('layers panel is not on screen');
  await page.mouse.click(panel.x + 11, panel.y + 11 + LAYER_ROW[layer] * 24);
}

const mapSignature = (page: Page) => v.viewerSignature(page, 'Map');
const mapRepaints = (page: Page, before: string) => v.waitForViewerRepaint(page, 'Map', before);
const shownValue = (page: Page, prop: string) => v.propertyGridValue(page, prop);
const ensureCategory = (page: Page, category: string, probeProp: string) =>
  v.ensurePropertyCategory(page, 'Map', category, probeProp);
const editProperty = (page: Page, prop: string, value: string) =>
  v.setPropertyGridValue(page, prop, value);

test('Map viewer', async ({page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(e.message));
  page.on('console', (m) => {
    if (m.type() === 'error') pageErrors.push(m.text());
  });

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 5000});

  // #### Add the viewer and auto-detect the geo columns
  await softStep('Add Map viewer from the Viewers toolbox', async () => {
    await page.locator('[name="icon-Map"]').first().click();
    await page.locator('[name="viewer-Map"]').waitFor({timeout: 30_000});

    // The geo columns are detected asynchronously once the map is up.
    await v.waitForPropertyValue(page, 'latitude', 'Latitude', undefined, 30_000);
    expect(await shownValue(page, 'longitude')).toBe('Longitude');

    const layers = await layerVisibility(page);
    expect(layers['Markers GL']).toBe(true);
    expect(layers['Heatmap']).toBe(false);
  });

  // #### Color and size coding through the Context Panel column selectors
  await softStep('Color by Magnitude', async () => {
    const before = await mapSignature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'color', columnName: 'Magnitude', viewerType: 'Map',
      propName: 'colorColumnName', scopeSelector: '.property-grid',
    });
    expect(await shownValue(page, 'color')).toBe('Magnitude');
    await mapRepaints(page, before);
  });

  await softStep('Size by Depth', async () => {
    const before = await mapSignature(page);
    await v.pickColumnViaSelectorTrusted(page, {
      role: 'size', columnName: 'Depth', viewerType: 'Map',
      propName: 'sizeColumnName', scopeSelector: '.property-grid',
    });
    expect(await shownValue(page, 'size')).toBe('Depth');
    await mapRepaints(page, before);
  });

  // #### Marker sizing through the property-grid editor
  await softStep('Marker Min Size redraws the markers', async () => {
    await ensureCategory(page, 'markers', 'marker-min-size');
    const before = await mapSignature(page);
    await editProperty(page, 'marker-min-size', '12');
    expect(await shownValue(page, 'marker-min-size')).toBe('12');
    await mapRepaints(page, before);
    await editProperty(page, 'marker-min-size', '2');
  });

  // #### Layers panel — the checkbox drives the OpenLayers stack
  await softStep('Layers panel toggles Heatmap and Markers GL', async () => {
    await openLayersPanel(page);

    const beforeHeatmap = await mapSignature(page);
    await toggleLayer(page, 'Heatmap');
    expect((await layerVisibility(page))['Heatmap']).toBe(true);
    await mapRepaints(page, beforeHeatmap);

    // Toggling Markers GL used to freeze the page (exponential re-subscription
    // of onCurrentCellChanged). The click must return and the next UI
    // interaction must still be served.
    await toggleLayer(page, 'Markers GL');
    expect((await layerVisibility(page))['Markers GL']).toBe(false);
    await page.locator('[name="button-Map-layers"]').click({timeout: 5000});
    await expect(page.locator('.panel-layers')).toBeHidden();

    await openLayersPanel(page);
    await toggleLayer(page, 'Markers GL');
    await toggleLayer(page, 'Heatmap');
    const restored = await layerVisibility(page);
    expect(restored['Markers GL']).toBe(true);
    expect(restored['Heatmap']).toBe(false);
    await page.locator('[name="button-Map-layers"]').click();
  });

  // #### Render type
  await softStep('Render Type cycles markers / heatmap / both', async () => {
    await ensureCategory(page, 'misc', 'render-type');
    const signatures: Record<string, string> = {};
    for (const mode of ['heatmap', 'both', 'markers']) {
      const before = await mapSignature(page);
      await page.locator('.property-grid tr[name="prop-render-type"] td').last().click();
      await page.locator('.property-grid select').selectOption(mode);
      await v.waitForPropertyValue(page, 'render-type', mode);
      signatures[mode] = await mapRepaints(page, before);
    }
    expect(signatures['heatmap']).not.toBe(signatures['markers']);
    expect(signatures['both']).not.toBe(signatures['markers']);
  });

  // #### Zoom controls
  await softStep('Zoom in and out with the map buttons', async () => {
    const zoom = () => page.evaluate(() => {
      const mv = Array.from((window as any).grok.shell.tv.viewers).find((x: any) => x.type === 'Map') as any;
      return mv.ol.olMap.getView().getZoom();
    });
    const start = await zoom();
    const before = await mapSignature(page);

    await page.locator('button.ol-zoom-in').click();
    await expect.poll(zoom, {timeout: 10_000}).toBeCloseTo(start + 1, 2);
    await mapRepaints(page, before);

    await page.locator('button.ol-zoom-out').click();
    await expect.poll(zoom, {timeout: 10_000}).toBeCloseTo(start, 2);
  });

  // #### Selection with Ctrl + drag, cleared with Escape
  await softStep('Ctrl+drag selects points, Escape clears', async () => {
    const box = (await page.locator('[name="viewer-Map"]').boundingBox())!;
    const before = await mapSignature(page);

    await page.keyboard.down('Control');
    await page.mouse.move(box.x + box.width * 0.25, box.y + box.height * 0.25);
    await page.mouse.down();
    await page.mouse.move(box.x + box.width * 0.75, box.y + box.height * 0.75, {steps: 25});
    await page.mouse.up();
    await page.keyboard.up('Control');

    const selected = () =>
      page.evaluate(() => (window as any).grok.shell.t.selection.trueCount as number);
    await expect.poll(selected, {timeout: 10_000}).toBeGreaterThan(0);
    await mapRepaints(page, before);

    await page.keyboard.press('Escape');
    await expect.poll(selected, {timeout: 10_000}).toBe(0);
  });

  // #### Only filtered rows stay on the map
  await softStep('Filtering the table narrows the map', async () => {
    const before = await mapSignature(page);
    const {filteredCount} = await v.applyCategoricalFilter(page, 'MagType', ['Mw']);
    const total = await page.evaluate(() => (window as any).grok.shell.t.rowCount);
    expect(filteredCount).toBeGreaterThan(0);
    expect(filteredCount).toBeLessThan(total);
    await mapRepaints(page, before);

    await v.resetFilters(page);
  });

  // #### Tooltip
  await softStep('Show Tooltip reveals a tooltip over a point', async () => {
    await ensureCategory(page, 'misc', 'show-tooltip');
    await page.locator('[name="prop-view-show-tooltip"]').click();

    // The tooltip only fires when the pointer is over a rendered marker, so the
    // hover targets come from the features themselves — screen coordinates of
    // markers that currently fall inside the viewport.
    const targets = await page.evaluate(() => {
      const mv = Array.from((window as any).grok.shell.tv.viewers).find((x: any) => x.type === 'Map') as any;
      const r = (mv.root as HTMLElement).getBoundingClientRect();
      const out: {x: number; y: number}[] = [];
      for (const f of (mv.features ?? []).slice(0, 400)) {
        let px: number[] | null = null;
        try { px = mv.ol.olMap.getPixelFromCoordinate(f.getGeometry().getCoordinates()); }
        catch { continue; }
        if (!px || px[0] < 5 || px[1] < 5 || px[0] > r.width - 5 || px[1] > r.height - 5) continue;
        out.push({x: r.x + px[0], y: r.y + px[1]});
        if (out.length >= 5) break;
      }
      return out;
    });
    expect(targets.length).toBeGreaterThan(0);

    let tooltip = {visible: false, text: ''};
    for (const t of targets) {
      await page.mouse.move(t.x - 15, t.y - 15);
      await page.mouse.move(t.x, t.y, {steps: 6});
      await page.waitForTimeout(400);
      tooltip = await page.evaluate(() => {
        const el = document.querySelector('.d4-tooltip') as HTMLElement | null;
        return {visible: !!el && getComputedStyle(el).display !== 'none', text: el?.innerText ?? ''};
      });
      if (tooltip.visible) break;
    }
    expect(tooltip.visible).toBe(true);
    expect(tooltip.text).toContain('Latitude');
    expect(tooltip.text).toContain('Magnitude');
  });

  // #### Closing the viewer must not raise
  await softStep('Close the viewer after moving the pointer over it', async () => {
    const box = (await page.locator('[name="viewer-Map"]').boundingBox())!;
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
    pageErrors.length = 0;

    await v.clickViewerTitlebarIcon(page, 'Map', 'Close');
    await expect(page.locator('[name="viewer-Map"]')).toHaveCount(0);

    // Pointer events arriving after disposal used to crash the viewer.
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2 + 40);
    await page.waitForTimeout(700);
    expect(pageErrors.join('\n')).toBe('');
  });

  await v.cleanupShell(page);

  v.finishSpec();
});
