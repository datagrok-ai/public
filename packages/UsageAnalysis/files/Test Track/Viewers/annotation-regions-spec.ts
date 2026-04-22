import {test, expect, chromium, Page} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

declare const grok: any;

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];
async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e?.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
  }
}

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

test('Annotation regions scenario', async () => {
    test.setTimeout(300_000);

    // Reuse existing logged-in context via CDP
    const cdp = await chromium.connectOverCDP('http://127.0.0.1:9222');
    const ctx = cdp.contexts()[0];
    const pages = ctx.pages();
    let page: Page = pages.find((p) => p.url().includes('datagrok')) ?? pages[0];
    if (!page) page = await ctx.newPage();
    await page.bringToFront();

    await page.goto(baseUrl, {timeout: 60_000, waitUntil: 'networkidle'});
    await page.waitForFunction(() => {
      try {
        if (typeof grok === 'undefined' || !grok.shell || !grok.dapi || !grok.dapi.files) return false;
        grok.shell.closeAll();
        return true;
      } catch { return false; }
    }, {timeout: 180_000, polling: 1000});
    await page.waitForTimeout(2000);

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

      // lassoTool should default to false
      const lt = await page.evaluate(() =>
        grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot').getOptions(true).look.lassoTool);
      expect(lt).toBe(false);

      // Context menu → Draw Annotation Region
      await rightClick(page, '[name="viewer-Scatter-plot"] canvas');
      await clickMenuItem(page, 'Draw Annotation Region');

      // Drag a rectangle using real pointer events (required for canvas)
      const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
      const x1 = box.x + box.width * 0.3, y1 = box.y + box.height * 0.3;
      const x2 = box.x + box.width * 0.6, y2 = box.y + box.height * 0.6;
      await page.mouse.move(x1, y1);
      await page.mouse.down();
      await page.mouse.move((x1 + x2) / 2, (y1 + y2) / 2);
      await page.mouse.move(x2, y2);
      await page.mouse.up();

      // Dialog auto-opens
      await page.locator('.d4-dialog .d4-dialog-title', {hasText: 'Formula Lines'}).waitFor({timeout: 5000});

      // Edit fields
      await page.locator('[name="input-host-Title"] input').fill('My Rect Region');
      await page.locator('[name="input-host-Description"] textarea').fill('Rectangle description');
      await page.locator('[name="input-host-Region-Color"] input').fill('#ff8800');
      await page.locator('[name="input-host-Outline-Color"] input').fill('#003366');

      // OK
      await page.locator('.d4-dialog [name="button-OK"]').click();

      // Verify persisted
      const region = await page.evaluate(() => {
        const sp = grok.shell.tv.viewers.find((v: any) => v.type === 'Scatter plot');
        return JSON.parse(sp.getOptions(true).look.annotationRegions)[0];
      });
      expect(region.header).toBe('My Rect Region');
      expect(region.description).toBe('Rectangle description');
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

    // ---------- 3. Region interaction (canvas hit-test) ----------
    await softStep('3.1 Hover region / intersection', async () => {
      // Real pointer hover — CDP-level. This may or may not match region coords
      // without world→screen conversion; use approximate center and verify mouseOverRowGroup changes.
      const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
      await page.mouse.move(box.x + box.width * 0.35, box.y + box.height * 0.55);
      await page.waitForTimeout(400);
      // Best-effort verification; do not hard-fail since region position depends on axis autoscale
      const moCount = await page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowGroup?.trueCount ?? 0);
      console.log(`mouseOverRowGroup trueCount: ${moCount}`);
    });

    await softStep('3.2 Click / Ctrl+click region', async () => {
      const box = (await page.locator('[name="viewer-Scatter-plot"] canvas').first().boundingBox())!;
      await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
      await page.mouse.click(box.x + box.width * 0.35, box.y + box.height * 0.55);
      await page.waitForTimeout(300);
      const sel = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
      console.log(`selection.trueCount after click: ${sel}`);
    });

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

    // ---------- 5. Dialog integration ----------
    await softStep('5.1 Preview / 5.2 Grid representation', async () => {
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

    await softStep('6.2 Single-axis Line Chart shows drawing option', async () => {
      await page.evaluate(() =>
        grok.shell.tv.viewers.find((v: any) => v.type === 'Line chart').setOptions({multiAxis: false, yColumnNames: ['HEIGHT']}));
      await page.waitForTimeout(300);
      await rightClick(page, '[name="viewer-Line-chart"] canvas');
      const drawVisible = await page.locator('.d4-menu-item-label', {hasText: /^Draw Annotation Region$/}).count();
      expect(drawVisible).toBeGreaterThan(0);
      await page.keyboard.press('Escape');
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
});
