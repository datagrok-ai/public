import {test, expect, chromium} from '@playwright/test';

const baseUrl = 'https://dev.datagrok.ai/';

test('Radar viewer (Charts package)', async () => {
  test.setTimeout(600_000);

  const browser = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find((p) => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
  }

  await page.waitForFunction(() => typeof (globalThis as any).grok !== 'undefined' && (globalThis as any).grok.shell, {timeout: 15000});

  const stepErrors: {step: string; error: string}[] = [];
  async function softStep(name: string, fn: () => Promise<void>) {
    try {
      await test.step(name, fn);
    } catch (e: any) {
      stepErrors.push({step: name, error: e?.message ?? String(e)});
      console.error(`[STEP FAILED] ${name}: ${e?.message ?? e}`);
    }
  }

  await softStep('Setup environment', async () => {
    await page!.evaluate(async () => {
      document.body.classList.add('selenium');
      (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
      (window as any).grok.shell.windows.simpleMode = false;
      (window as any).grok.shell.closeAll();
    });
  });

  await softStep('Step 1 — Open earthquakes.csv and add Radar viewer', async () => {
    await page!.evaluate(async () => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/geo/earthquakes.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
    });
    await page!.locator('[name="viewer-Grid"]').first().waitFor({timeout: 30000});

    await page!.evaluate(() => {
      const el = document.querySelector('i.svg-add-viewer') as HTMLElement;
      el.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      el.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      el.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page!.locator('.d4-dialog').waitFor({timeout: 5000});
    await page!.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog')!;
      const cards = dialog.querySelectorAll('.d4-item-card.viewer-gallery');
      for (const card of Array.from(cards)) {
        const lbl = card.querySelector('.card-label');
        if (lbl && lbl.textContent?.trim() === 'Radar') {
          (card as HTMLElement).dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
          (card as HTMLElement).dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
          (card as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
          (card as HTMLElement).dispatchEvent(new MouseEvent('dblclick', {bubbles: true}));
          break;
        }
      }
    });
    await page!.locator('[name="viewer-Radar"]').waitFor({timeout: 15000});
  });

  await softStep('Step 2 — Open demog.csv and add Radar viewer', async () => {
    await page!.evaluate(async () => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
    });
    await page!.evaluate(() => {
      const el = document.querySelector('i.svg-add-viewer') as HTMLElement;
      el.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      el.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      el.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page!.locator('.d4-dialog').waitFor({timeout: 5000});
    await page!.evaluate(() => {
      const dialog = document.querySelector('.d4-dialog')!;
      const cards = dialog.querySelectorAll('.d4-item-card.viewer-gallery');
      for (const card of Array.from(cards)) {
        const lbl = card.querySelector('.card-label');
        if (lbl && lbl.textContent?.trim() === 'Radar') {
          (card as HTMLElement).dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
          (card as HTMLElement).dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
          (card as HTMLElement).dispatchEvent(new MouseEvent('click', {bubbles: true}));
          (card as HTMLElement).dispatchEvent(new MouseEvent('dblclick', {bubbles: true}));
          break;
        }
      }
    });
    await page!.locator('[name="viewer-Radar"]').waitFor({timeout: 15000});
    const viewerTypes = await page!.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      return types;
    });
    expect(viewerTypes).toContain('Radar');
  });

  await softStep('Step 3 — Open Radar settings (gear icon)', async () => {
    await page!.evaluate(() => {
      const radar = document.querySelector('[name="viewer-Radar"]')!;
      const panelBase = radar.closest('.panel-base')!;
      const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      gear.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      gear.dispatchEvent(new MouseEvent('click', {bubbles: true}));
    });
    await page!.waitForTimeout(800);
    const categories = await page!.evaluate(() => {
      const pg = document.querySelector('.property-grid-base');
      if (!pg) return [];
      return Array.from(pg.querySelectorAll('.property-grid-category')).map((c) => c.textContent?.trim());
    });
    expect(categories).toEqual(expect.arrayContaining(['Data', 'Selection', 'Color', 'Value', 'Style']));
  });

  await softStep('Step 3a — Switch tables (Table dropdown)', async () => {
    const result = await page!.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      radar.setOptions({table: 'Table'});
      await new Promise((r) => setTimeout(r, 1500));
      const switched = radar.dataFrame?.name;
      radar.setOptions({table: 'Table (2)'});
      await new Promise((r) => setTimeout(r, 1500));
      const restored = radar.dataFrame?.name;
      return {switched, restored};
    });
    expect(result.switched).toBe('Table');
    expect(result.restored).toBe('Table (2)');
  });

  await softStep('Step 3b — Selection check-boxes', async () => {
    const result = await page!.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      radar.setOptions({showCurrentRow: false, showMouseOverRow: false, showMouseOverRowGroup: false});
      await new Promise((r) => setTimeout(r, 500));
      const off = radar.getOptions().look;
      radar.setOptions({showCurrentRow: true, showMouseOverRow: true, showMouseOverRowGroup: true});
      await new Promise((r) => setTimeout(r, 500));
      const on = radar.getOptions().look;
      return {off, on};
    });
    expect(result.off.showCurrentRow).toBe(false);
    expect(result.on.showCurrentRow).toBe(true);
  });

  await softStep('Step 3c — Increase/decrease Values', async () => {
    const result = await page!.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      radar.setOptions({valuesColumnNames: ['AGE', 'HEIGHT']});
      await new Promise((r) => setTimeout(r, 800));
      const small = radar.getOptions().look.valuesColumnNames;
      radar.setOptions({valuesColumnNames: ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']});
      await new Promise((r) => setTimeout(r, 800));
      const big = radar.getOptions().look.valuesColumnNames;
      return {small, big};
    });
    expect(result.small.length).toBe(2);
    expect(result.big.length).toBe(4);
  });

  await softStep('Step 3d — Style (color) changes', async () => {
    const result = await page!.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      let radar: any = null;
      for (const v of tv.viewers) if (v.type === 'Radar') { radar = v; break; }
      radar.setOptions({backgroundMinColor: 0xFFFFE0E0, backgroundMaxColor: 0xFFFF6666, lineColor: 0xFF0000FF});
      await new Promise((r) => setTimeout(r, 800));
      const after = radar.getOptions().look;
      radar.setOptions({backgroundMinColor: 4290479197, backgroundMaxColor: 4293381581, lineColor: 11393254});
      return after;
    });
    expect(result.lineColor).toBe(0xFF0000FF);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
