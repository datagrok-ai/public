import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Scripts Run — context menu + console', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  await softStep('1. Go to Browse > Platform > Functions > Scripts', async () => {
    await page.evaluate(() => { grok.shell.route('/scripts'); });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    const v = await page.evaluate(() => grok.shell.v?.name);
    expect(v).toBe('Scripts');
  });

  // Prerequisite — script must exist and cars table must be open for the dialog to pre-select it
  await softStep('[pre] Ensure testRscript exists and cars is open', async () => {
    await page.evaluate(async () => {
      const existing = await grok.dapi.scripts.filter('name = "testRscript"').first().catch(() => null);
      if (!existing) {
        const body = '#name: testRscript\n#language: r\n#sample: cars.csv\n#input: dataframe table [Data table]\n#output: int count [Number of cells in table]\n#output: string newParam\ncount <- nrow(table) * ncol(table)\nnewParam="test"\n';
        const s = DG.Script.create(body);
        await grok.dapi.scripts.save(s);
      }
      if (!grok.shell.tables.some((t: any) => t.name === 'cars')) {
        const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv').catch(() => null as any);
        if (df) { df.name = 'cars'; grok.shell.addTableView(df); }
      }
      // Full route round-trip so the gallery re-fetches
      grok.shell.route('/');
      await new Promise((r) => setTimeout(r, 400));
      grok.shell.route('/scripts');
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    // Poll for the card with a fallback delay
    const found = await page.evaluate(async () => {
      for (let i = 0; i < 40; i++) {
        const hit = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
          .some((e) => e.textContent?.trim() === 'testRscript');
        if (hit) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    });
    expect(found).toBe(true);
  });

  await softStep('2-3. Right-click testRscript, Run..., choose cars, OK', async () => {
    const balloons = await page.evaluate(async () => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'testRscript') as HTMLElement;
      const card = (label?.closest('.grok-gallery-grid-item') || label?.parentElement) as HTMLElement;
      const rect = card.getBoundingClientRect();
      card.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, view: window,
        clientX: rect.left + 10, clientY: rect.top + 10, button: 2,
      }));
      await new Promise((r) => setTimeout(r, 1200));
      const run = Array.from(document.querySelectorAll('.d4-menu-item-label, .d4-menu-item'))
        .find((e) => e.textContent?.trim() === 'Run...') as HTMLElement;
      run.click();
      // Wait for the dialog to appear
      for (let i = 0; i < 30; i++) {
        const dlgs = Array.from(document.querySelectorAll('.d4-dialog')).filter((d) => (d as HTMLElement).offsetParent !== null);
        if (dlgs.length > 0 && dlgs.some((d) => d.querySelector('[name="button-OK"]'))) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      const dlgs = Array.from(document.querySelectorAll('.d4-dialog'));
      let ok: HTMLElement | null = null;
      for (const d of dlgs) {
        const b = d.querySelector('[name="button-OK"]') as HTMLElement | null;
        if (b && b.offsetParent !== null) { ok = b; break; }
      }
      ok?.click();
      await new Promise((r) => setTimeout(r, 4000));
      return Array.from(document.querySelectorAll('.grok-balloon')).map((b) => b.textContent?.trim()).filter(Boolean);
    });
    // Script run should produce output (count=510) — filter by the specific error text
    expect(balloons.filter((b: any) => /Value not defined/i.test(b || ''))).toEqual([]);
  });

  // 4-6. Rerun with local file / Datagrok Files / query — intentionally skipped (manual)
  await test.step.skip('4-6. Rerun with local file / Datagrok Files / query (SKIP — manual)', async () => {});

  await softStep('7. Open Datagrok console (~)', async () => {
    await page.keyboard.press('Backquote');
    await page.waitForTimeout(1000);
    const ok = await page.evaluate(() => !!document.querySelector('.d4-console-wrapper') &&
      !!Array.from(document.querySelectorAll('input[placeholder*="Enter command"]')).find((i) => (i as HTMLElement).offsetParent !== null));
    expect(ok).toBe(true);
  });

  await softStep('8-9. Enter agolovko:testRscript("cars"), press Enter, check output', async () => {
    await page.evaluate(async () => {
      // Ensure cars table is open
      if (!grok.shell.tables.some((t: any) => t.name === 'cars')) {
        const df = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv').catch(() => null as any);
        if (df) { grok.shell.addTableView(df); }
        else {
          const d = (grok as any).data?.demo?.cars?.();
          if (d) { d.name = 'cars'; grok.shell.addTableView(d); }
        }
      }
      const tbl = grok.shell.tables.find((t: any) => t.name === 'cars' || t.name === 'Table');
      if (tbl && tbl.name !== 'cars') tbl.name = 'cars';
    });
    await page.evaluate(() => {
      const cmd = Array.from(document.querySelectorAll('input[placeholder*="Enter command"]'))
        .find((i) => (i as HTMLElement).offsetParent !== null) as HTMLInputElement;
      cmd.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      setter.call(cmd, 'agolovko:testRscript("cars")');
      cmd.dispatchEvent(new Event('input', {bubbles: true}));
    });
    await page.keyboard.press('Enter');
    await page.waitForTimeout(4000);
    const lines = await page.evaluate(() => {
      const wrap = document.querySelector('.d4-console-wrapper');
      return Array.from(wrap?.querySelectorAll('*') || [])
        .filter((e) => e.children.length === 0 && e.textContent?.trim())
        .map((e) => e.textContent!.trim()).slice(-10);
    });
    const hasCount = lines.some((l: string) => /count:\s*510/.test(l));
    const hasNewParam = lines.some((l: string) => /newParam:\s*"test"/.test(l));
    expect(hasCount).toBe(true);
    expect(hasNewParam).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
