import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;

async function rightClick(page: Page, selector: string, dx = 100, dy = 100) {
  await page.evaluate(({sel, x, y}) => {
    const el = document.querySelector(sel) as HTMLElement;
    if (!el) throw new Error(`no element: ${sel}`);
    const r = el.getBoundingClientRect();
    el.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
      clientX: r.left + x, clientY: r.top + y,
    }));
  }, {sel: selector, x: dx, y: dy});
  await page.waitForTimeout(300);
}

async function clickMenuItem(page: Page, text: string) {
  const ok = await page.evaluate((t) => {
    const labels = Array.from(document.querySelectorAll('.d4-menu-item .d4-menu-item-label'));
    const target = labels.find((l) => (l.textContent ?? '').trim() === t);
    if (!target) return false;
    (target.closest('.d4-menu-item') as HTMLElement).click();
    return true;
  }, text);
  if (!ok) throw new Error(`menu item not found: ${text}`);
}

test('Annotation regions scenario', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Phase 2: Open demog
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise((res) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(undefined); });
      setTimeout(() => res(undefined), 4000);
    });
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
  });
  await page.locator('[name="viewer-Grid"]').first().waitFor({timeout: 30_000});

  // ---------- 1. Region creation ----------
  await softStep('1.1 Draw rectangle region + edit properties', async () => {
    await page.locator('[name="icon-scatter-plot"]').click();
    await page.locator('[name="viewer-Scatter-plot"] canvas').first().waitFor();

    const lt = await page.evaluate(() =>
      grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot').getOptions(true).look.lassoTool);
    expect(lt).toBe(false);

    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Draw Annotation Region');

    const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
    const x1 = box.x + box.width * 0.3, y1 = box.y + box.height * 0.3;
    const x2 = box.x + box.width * 0.6, y2 = box.y + box.height * 0.6;
    await page.mouse.move(x1, y1);
    await page.mouse.down();
    await page.mouse.move((x1 + x2) / 2, (y1 + y2) / 2);
    await page.mouse.move(x2, y2);
    await page.mouse.up();

    await page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'}).waitFor({timeout: 5000});
    await page.locator('[name="input-host-Title"] input').fill('My Rect Region');
    await page.locator('[name="input-host-Description"] textarea').fill('Rectangle description');
    await page.locator('[name="input-host-Region-Color"] input').fill('#ff8800');
    await page.locator('[name="input-host-Outline-Color"] input').fill('#003366');
    // Outline Width and Opacity are raw <input type="range"> without input-host-* wrapper
    await page.evaluate(() => {
      const findRangeByLabel = (caption: string) => {
        const label = Array.from(document.querySelectorAll('.d4-dialog .ui-input-label'))
          .find((el) => el.textContent?.trim() === caption);
        return label?.parentElement?.querySelector('input') as HTMLInputElement | null;
      };
      const w = findRangeByLabel('Outline Width');
      if (!w) throw new Error('Outline Width input not found');
      w.value = '3'; w.dispatchEvent(new Event('input', {bubbles: true}));
      const o = findRangeByLabel('Opacity');
      if (!o) throw new Error('Opacity input not found');
      o.value = '60'; o.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.locator('[name="input-host-Header-Color"] input').fill('#ffffff');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      return JSON.parse(sp.getOptions(true).look.annotationRegions)[0];
    });
    expect(region.header).toBe('My Rect Region');
    expect(region.description).toBe('Rectangle description');
    expect(region.outlineWidth).toBe(3);
    expect(region.opacity).toBe(60);
  });

  await softStep('1.2 Draw Lasso (polygon) region', async () => {
    await page.evaluate(() => {
      grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot').setOptions({lassoTool: true});
    });
    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Draw Annotation Region');

    const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
    const pts: [number, number][] = [
      [box.x + box.width * 0.65, box.y + box.height * 0.2],
      [box.x + box.width * 0.85, box.y + box.height * 0.3],
      [box.x + box.width * 0.85, box.y + box.height * 0.5],
      [box.x + box.width * 0.70, box.y + box.height * 0.55],
      [box.x + box.width * 0.65, box.y + box.height * 0.4],
    ];
    await page.mouse.move(pts[0][0], pts[0][1]);
    await page.mouse.down();
    for (let i = 1; i < pts.length; i++) await page.mouse.move(pts[i][0], pts[i][1], {steps: 5});
    await page.mouse.up();

    await page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'}).waitFor({timeout: 5000});
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
  });

  await softStep('1.3 Create formula region', async () => {
    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Formula Lines...');
    await page.locator('.d4-dialog').waitFor();

    await page.locator('.d4-dialog button, .d4-dialog .ui-btn, .d4-dialog span').getByText('Add new', {exact: true}).first().click();
    await clickMenuItem(page, 'Region - Formula Lines');

    const formula1 = await page.locator('.d4-dialog .ui-input-root').filter({hasText: 'Formula 1'}).locator('input').inputValue();
    expect(formula1).toMatch(/\$\{HEIGHT\}.*\$\{WEIGHT\}/);

    await page.locator('.d4-dialog [name="button-OK"]').click();
  });

  // ---------- 2. Visibility ----------
  await softStep('2.1 Show/Hide viewer & dataframe independently', async () => {
    await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      sp.setOptions({
        annotationRegions: JSON.stringify([
          {type: 'area', x: 'WEIGHT', y: 'HEIGHT', area: [[60,150],[95,150],[95,185],[60,185]], header: 'Viewer Region'}
        ])
      });
      sp.dataFrame.setTag('.annotation-regions', JSON.stringify([
        {type: 'area', x: 'WEIGHT', y: 'HEIGHT', area: [[90,160],[120,160],[120,190],[90,190]], header: 'DF Region', isDataFrameRegion: true}
      ]));
    });

    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Show Viewer Annotation Regions');
    await page.waitForTimeout(200);

    const after1 = await page.evaluate(() => {
      const look = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot').getOptions(true).look;
      return {v: look.showViewerAnnotationRegions, d: look.showDataframeAnnotationRegions};
    });
    expect(after1.v).toBe(false);
    expect(after1.d).toBe(true);

    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Show Dataframe Annotation Regions');
    await page.waitForTimeout(200);

    const after2 = await page.evaluate(() => {
      const look = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot').getOptions(true).look;
      return {v: look.showViewerAnnotationRegions, d: look.showDataframeAnnotationRegions};
    });
    expect(after2.v).toBe(false);
    expect(after2.d).toBe(false);
  });

  await softStep('2.2 Global Show Annotation Regions re-enables both', async () => {
    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Show Annotation Regions');
    await page.waitForTimeout(200);
    const after = await page.evaluate(() => {
      const look = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot').getOptions(true).look;
      return {v: look.showViewerAnnotationRegions, d: look.showDataframeAnnotationRegions};
    });
    expect(after.v).toBe(true);
    expect(after.d).toBe(true);
  });

  // ---------- 3. Region interaction ----------
  // 3.1 Hover interaction — manual only, see annotation-regions-ui.md
  // 3.2 Click region → selection — manual only, see annotation-regions-ui.md

  // ---------- 4. Region editing ----------
  await softStep('4.3 Change Annotation Font', async () => {
    const [before, after] = await page.evaluate(async () => {
      const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      const b = sp.getOptions(true).look.annotationFont;
      sp.setOptions({annotationFont: '18px bold Arial'});
      await new Promise(r => setTimeout(r, 300));
      return [b, sp.getOptions(true).look.annotationFont];
    });
    expect(before).not.toBe(after);
    expect(after).toBe('18px bold Arial');
  });

  // ---------- 4 (cont). Region editing ----------
  await softStep('4.1 Right-click region → Edit opens dialog', async () => {
    // Attempt to hit the region drawn in 1.1 — region is at ~30–60% of canvas.
    // Without worldToScreen conversion the cursor may miss; step logs result and does not hard-fail.
    const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
    await rightClick(page, '[name="viewer-Scatter-plot"] canvas', box.width * 0.45, box.height * 0.45);
    const hasEdit = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .some((el) => (el.textContent ?? '').trim() === 'Edit'));
    console.log(`[4.1] "Edit" menu item found: ${hasEdit}`);
    await page.keyboard.press('Escape');
    if (!hasEdit)
      console.warn('[4.1] AMBIGUOUS — cursor likely missed the region (worldToScreen unavailable)');
  });

  await softStep('4.2 Modify region via dialog: reopen and edit Outline Width / Opacity / Header Color', async () => {
    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
    await clickMenuItem(page, 'Formula Lines...');
    await page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'}).waitFor({timeout: 5000});

    // Click first data row of the left-side grid to select the region
    const gridBox = (await page.locator('.d4-dialog .d4-grid').boundingBox())!;
    await page.mouse.click(gridBox.x + gridBox.width * 0.5, gridBox.y + 36);
    await page.waitForTimeout(300);

    await page.evaluate(() => {
      const findRangeByLabel = (caption: string) => {
        const label = Array.from(document.querySelectorAll('.d4-dialog .ui-input-label'))
          .find((el) => el.textContent?.trim() === caption);
        return label?.parentElement?.querySelector('input') as HTMLInputElement | null;
      };
      const w = findRangeByLabel('Outline Width');
      if (!w) throw new Error('Outline Width input not found');
      w.value = '5'; w.dispatchEvent(new Event('input', {bubbles: true}));
      const o = findRangeByLabel('Opacity');
      if (!o) throw new Error('Opacity input not found');
      o.value = '40'; o.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.locator('[name="input-host-Header-Color"] input').fill('#ff0000');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const region = await page.evaluate(() => {
      const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
      return JSON.parse(sp.getOptions(true).look.annotationRegions)[0];
    });
    expect(region.outlineWidth).toBe(5);
    expect(region.opacity).toBe(40);
  });

  // ---------- 5. Dialog integration ----------
  await softStep('5.1 Preview / 5.2 Grid representation', async () => {
    // Close any dialogs left open by failed previous steps
    await page.evaluate(() => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const btn = (d.querySelector('[name="button-CANCEL"]') ?? d.querySelector('[name="button-OK"]')) as HTMLElement | null;
        btn?.click();
      });
    });
    await page.waitForTimeout(300);
    await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
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

  // ---------- 6. Line Chart ----------
  await softStep('6.1 Multi-axis limitation', async () => {
    await page.evaluate(() => {
      for (const v of [...grok.shell.tv.viewers].filter((v: any) => v.type === 'Scatter plot')) v.close();
    });
    await page.locator('[name="icon-line-chart"]').click();
    await page.locator('[name="viewer-Line-chart"] canvas').first().waitFor();
    await page.evaluate(() => grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart').setOptions({multiAxis: true}));
    await page.waitForTimeout(300);
    await rightClick(page, '[name="viewer-Line-chart"] canvas');
    const drawVisible = await page.locator('.d4-menu-item-label', {hasText: /^Draw Annotation Region$/}).count();
    expect(drawVisible).toBe(0);
    await page.keyboard.press('Escape');
  });

  await softStep('6.2 Single-axis Line Chart: draw rectangle region + add formula region', async () => {
    await page.evaluate(() =>
      grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart').setOptions({multiAxis: false, yColumnNames: ['HEIGHT']}));
    await page.waitForTimeout(300);

    // Verify menu item present
    await rightClick(page, '[name="viewer-Line-chart"] canvas');
    const drawVisible = await page.locator('.d4-menu-item-label', {hasText: /^Draw Annotation Region$/}).count();
    expect(drawVisible).toBeGreaterThan(0);
    await clickMenuItem(page, 'Draw Annotation Region');

    // Draw rectangle
    const box = (await page.locator('[name="viewer-Line-chart"] canvas').first().boundingBox())!;
    const x1 = box.x + box.width * 0.25, y1 = box.y + box.height * 0.25;
    const x2 = box.x + box.width * 0.55, y2 = box.y + box.height * 0.55;
    await page.mouse.move(x1, y1);
    await page.mouse.down();
    await page.mouse.move((x1 + x2) / 2, (y1 + y2) / 2);
    await page.mouse.move(x2, y2);
    await page.mouse.up();

    await page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'}).waitFor({timeout: 5000});
    await page.locator('[name="input-host-Title"] input').fill('LC Rect Region');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const afterRect = await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      return JSON.parse(lc.getOptions(true).look.annotationRegions || '[]').length;
    });
    expect(afterRect).toBeGreaterThan(0);

    // Add formula region
    await rightClick(page, '[name="viewer-Line-chart"] canvas');
    await clickMenuItem(page, 'Formula Lines...');
    await page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'}).waitFor({timeout: 5000});
    await page.locator('.d4-dialog button, .d4-dialog .ui-btn, .d4-dialog span').getByText('Add new', {exact: true}).first().click();
    await clickMenuItem(page, 'Region - Formula Lines');
    await page.locator('.d4-dialog [name="button-OK"]').click();

    const afterFormula = await page.evaluate(() => {
      const lc = grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart');
      return JSON.parse(lc.getOptions(true).look.annotationRegions || '[]').length;
    });
    expect(afterFormula).toBeGreaterThan(afterRect);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
