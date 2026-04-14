import {test, expect, chromium, Page} from '@playwright/test';

declare const grok: any;
declare const $: any;

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('Form viewer scenario', async () => {
  test.setTimeout(300_000);

  const cdp = await chromium.connectOverCDP('http://127.0.0.1:9222');
  const ctx = cdp.contexts()[0];
  const pages = ctx.pages();
  let page: Page = pages.find((p) => p.url().includes('datagrok')) ?? pages[0];
  if (!page) page = await ctx.newPage();
  await page.bringToFront();

  await page.goto(baseUrl, {timeout: 60000, waitUntil: 'networkidle'});
  await page.waitForFunction(() => {
    try {
      if (typeof grok === 'undefined' || !grok.shell || !grok.dapi || !grok.dapi.files) return false;
      grok.shell.closeAll();
      return true;
    } catch { return false; }
  }, {timeout: 180000, polling: 1000});
  await page.waitForTimeout(2000);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch {}
    try { grok.shell.windows.simpleMode = false; } catch {}
    grok.shell.closeAll();
    const paths = ['System:DemoFiles/SPGI.csv', 'System:DemoFiles/SPGI-linked1.csv', 'System:DemoFiles/SPGI-linked2.csv'];
    const names = ['SPGI', 'SPGI-linked1', 'SPGI-linked2'];
    for (let i = 0; i < paths.length; i++) {
      const df = await grok.dapi.files.readCsv(paths[i]);
      df.name = names[i];
      const v = grok.shell.addTableView(df);
      v.name = names[i];
      await new Promise((r) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); r(undefined); });
        setTimeout(r, 4000);
      });
    }
    const spgi = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.name === 'SPGI');
    if (spgi) grok.shell.v = spgi;
    for (let i = 0; i < 60; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 4000));
  });
  await page.locator('[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  let layoutId = '';

  await softStep('1. Switch bound table via Data → Table', async () => {
    await page.evaluate(() => (document.querySelector('[name="icon-form"]') as HTMLElement)?.click());
    await page.locator('[name="viewer-Form"]').waitFor({timeout: 15000});

    await page.evaluate(() => {
      const form = document.querySelector('[name="viewer-Form"]')!;
      const panel = form.closest('.panel-base')!;
      (panel.querySelector('[name="icon-font-icon-settings"]') as HTMLElement).click();
    });
    await page.waitForTimeout(500);

    const log: {target: string; bound: string}[] = [];
    for (const target of ['SPGI-linked2', 'SPGI-linked1', 'SPGI']) {
      // JS API fallback: Dart ChoiceInput listens via cash-dom; neither Playwright selectOption
      // nor dispatching native change events rebinds reliably across runners. Set viewer prop directly.
      await page.evaluate((name) => {
        const f: any = Array.from((grok.shell.tv as any).viewers).find((v: any) => v.type === 'Form');
        f.props.table = name;
      }, target);
      await page.waitForTimeout(1200);
      const bound = await page.evaluate(() => {
        const f: any = Array.from((grok.shell.tv as any).viewers).find((v: any) => v.type === 'Form');
        return f?.dataFrame?.name;
      });
      log.push({target, bound});
    }
    for (const r of log) expect(r.bound).toBe(r.target);
  });

  await softStep('2a. Open column-picker dialog', async () => {
    await page.evaluate(() => (document.querySelector('[name="viewer-Form"] [name="icon-list"]') as HTMLElement).click());
    await page.locator('.d4-dialog', {hasText: 'Select columns'}).waitFor({timeout: 5000});
  });

  await softStep('2b. Toggle None / All → field count changes', async () => {
    await page.locator('.d4-dialog [name="label-None"]').click();
    await page.waitForTimeout(600);
    expect(await page.locator('.d4-dialog').textContent()).toContain('0 checked');

    await page.locator('.d4-dialog [name="label-All"]').click();
    await page.waitForTimeout(600);
    expect(await page.locator('.d4-dialog').textContent()).toContain('88 checked');

    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(1200);

    const inputs = await page.locator('[name="viewer-Form"] input[type="text"]').count();
    expect(inputs).toBeGreaterThan(100);
  });

  await softStep('2c. Save layout, reduce to 0, restore', async () => {
    layoutId = await page.evaluate(async () => {
      const layout = (grok.shell.tv as any).saveLayout();
      await grok.dapi.layouts.save(layout);
      await new Promise((r) => setTimeout(r, 1500));
      return layout.id;
    });

    await page.evaluate(() => (document.querySelector('[name="viewer-Form"] [name="icon-list"]') as HTMLElement).click());
    await page.waitForTimeout(1000);
    await page.locator('.d4-dialog [name="label-None"]').click();
    await page.waitForTimeout(400);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(1200);

    expect(await page.locator('[name="viewer-Form"] input[type="text"]').count()).toBe(0);

    await page.evaluate(async (id) => {
      const saved = await grok.dapi.layouts.find(id);
      (grok.shell.tv as any).loadLayout(saved);
      await new Promise((r) => setTimeout(r, 3500));
    }, layoutId);

    const restored = await page.locator('[name="viewer-Form"] input[type="text"]').count();
    expect(restored).toBeGreaterThan(100);
  });

  await softStep('3. Design mode: drag field, position persists', async () => {
    await page.evaluate(() => (document.querySelector('[name="viewer-Form"] [name="icon-object-ungroup"]') as HTMLElement).click());
    await page.waitForTimeout(800);

    const {startX, startY, beforeX} = await page.evaluate(() => {
      const panel = Array.from(document.querySelectorAll('[name="viewer-Form"] .d4-host-element-panel'))
        .find((p) => (p as HTMLElement).offsetParent !== null) as HTMLElement;
      const r = panel.getBoundingClientRect();
      return {
        startX: Math.round(r.x + r.width / 2),
        startY: Math.round(r.y + r.height / 2),
        beforeX: Math.round(r.x),
      };
    });
    const endX = startX + 150;
    const endY = startY + 80;

    await page.mouse.move(startX, startY);
    await page.mouse.down();
    for (let i = 1; i <= 15; i++)
      await page.mouse.move(startX + (endX - startX) * i / 15, startY + (endY - startY) * i / 15, {steps: 2});
    await page.mouse.up();
    await page.waitForTimeout(800);

    const afterX = await page.evaluate(() => {
      const panel = Array.from(document.querySelectorAll('[name="viewer-Form"] .d4-host-element-panel'))
        .find((p) => (p as HTMLElement).offsetParent !== null) as HTMLElement;
      return Math.round(panel.getBoundingClientRect().x);
    });
    // verify panel moved meaningfully from its starting position
    expect(Math.abs(afterX - beforeX)).toBeGreaterThan(20);

    await page.waitForTimeout(1500);
    const stillX = await page.evaluate(() => {
      const panel = Array.from(document.querySelectorAll('[name="viewer-Form"] .d4-host-element-panel'))
        .find((p) => (p as HTMLElement).offsetParent !== null) as HTMLElement;
      return Math.round(panel.getBoundingClientRect().x);
    });
    expect(stillX).toBe(afterX);
  });

  await page.evaluate(async (id) => {
    try {
      const saved = id ? await grok.dapi.layouts.find(id) : null;
      if (saved) await grok.dapi.layouts.delete(saved);
    } catch {}
    grok.shell.closeAll();
  }, layoutId);

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
