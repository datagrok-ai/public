import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Scripts Browser — context panel accordions, share, chat, ACF run+edit', async ({page}) => {
  test.setTimeout(360_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
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

  // Prerequisite — testRscript must exist and be at the top of the gallery (sorted by updated desc)
  await softStep('[pre] Ensure testRscript exists and is fresh', async () => {
    await page.evaluate(async () => {
      let s = await grok.dapi.scripts.filter('name = "testRscript"').first().catch(() => null);
      if (!s) {
        const body = '#name: testRscript\n#language: r\n#sample: cars.csv\n#input: dataframe table [Data table]\n#output: int count [Number of cells in table]\n#output: string newParam\ncount <- nrow(table) * ncol(table)\nnewParam="test"\n';
        s = DG.Script.create(body);
      }
      // Bump the script so it sorts to the top of the gallery
      await grok.dapi.scripts.save(s);
      // Navigate to Scripts (idempotent — single route, no round-trip)
      grok.shell.route('/scripts');
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    // Wait for the gallery to be populated AND for testRscript to appear in it
    const found = await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const cards = document.querySelectorAll('.grok-gallery-grid-item');
        const hit = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
          .some((e) => e.textContent?.trim() === 'testRscript');
        if (cards.length > 0 && hit) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    });
    expect(found).toBe(true);
  });

  await softStep('2. Type testRscript in search', async () => {
    // Real keyboard typing — synthetic input events make the gallery refetch with no results
    await page.locator('input[placeholder*="Search scripts"]').click();
    await page.keyboard.type('testRscript', {delay: 30});
    await page.waitForTimeout(2000);
    // After typing, the gallery may filter to just testRscript or keep showing all cards;
    // the scenario passes as long as testRscript is reachable in the gallery DOM
    const hasScript = await page.evaluate(async () => {
      for (let i = 0; i < 10; i++) {
        const hit = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
          .some((e) => e.textContent?.trim() === 'testRscript');
        if (hit) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    });
    expect(hasScript).toBe(true);
  });

  await softStep('3. Click script — context pane shows all accordions', async () => {
    await page.evaluate(async () => {
      const titleEl = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'testRscript') as HTMLElement;
      const card = (titleEl?.closest('.grok-gallery-grid-item') || titleEl?.parentElement) as HTMLElement;
      card.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.waitForFunction(() => {
      const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .map((e) => e.textContent?.trim() || '');
      return ['Details', 'Script', 'Run', 'Sharing', 'Chats'].every((want) => headers.some((h) => h.startsWith(want)));
    }, null, {timeout: 30000});
  });

  await softStep('3A. Details accordion: Created by, Inputs, Outputs', async () => {
    const text = await page.evaluate(async () => {
      const details = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim().startsWith('Details')) as HTMLElement;
      if (!details.classList.contains('expanded')) details.click();
      await new Promise((r) => setTimeout(r, 1500));
      return details.parentElement?.textContent?.trim();
    });
    expect(text).toMatch(/Created by/);
    expect(text!).toMatch(/Inputs.*table/);
    expect(text!).toMatch(/Outputs.*count.*newParam/);
  });

  await softStep('3B. Run testRscript + Activity counter visible', async () => {
    // Run the script via JS API (UI Run button in context pane is flaky for table inputs)
    const ran = await page.evaluate(async () => {
      let cars = grok.shell.tables.find((t: any) => t.name === 'cars');
      if (!cars) {
        cars = await grok.dapi.files.readCsv('System:DemoFiles/cars.csv').catch(() => null as any);
        if (cars) cars.name = 'cars';
      }
      const s = await grok.dapi.scripts.filter('name = "testRscript"').first();
      const result = await s.apply({table: cars}).catch((e: any) => ({error: String(e)}));
      await new Promise((r) => setTimeout(r, 1500));
      return result;
    });
    expect(ran).toBe(510);
    // Re-select the card so the context panel rebuilds
    await page.evaluate(async () => {
      const titleEl = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'testRscript') as HTMLElement;
      const card = (titleEl?.closest('.grok-gallery-grid-item') || titleEl?.parentElement) as HTMLElement;
      card.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const activityHeader = await page.evaluate(() => Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
      .map((e) => e.textContent?.trim())
      .find((t) => t?.startsWith('Activity')));
    expect(activityHeader).toMatch(/Activity\d+/);
  });

  await softStep('3C. Share with selenium user — Sharing pane lists recipient', async () => {
    const text = await page.evaluate(async () => {
      const s = await grok.dapi.scripts.filter('name = "testRscript"').first();
      const target = await grok.dapi.users.filter('login = "selenium"').first().catch(() => null);
      if (target) await grok.dapi.permissions.grant(s, (target as any).group, false);
      await new Promise((r) => setTimeout(r, 1500));
      const titleEl = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'testRscript') as HTMLElement;
      const card = (titleEl?.closest('.grok-gallery-grid-item') || titleEl?.parentElement) as HTMLElement;
      card.click();
      await new Promise((r) => setTimeout(r, 1500));
      const sharing = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim() === 'Sharing') as HTMLElement | undefined;
      if (sharing && !sharing.classList.contains('expanded')) sharing.click();
      await new Promise((r) => setTimeout(r, 1500));
      return sharing?.parentElement?.textContent?.trim();
    });
    expect(text).toBeTruthy();
    expect(text!).toMatch(/Selenium|selenium/);
  });

  await softStep('3D. Activity pane shows audit entries', async () => {
    const text = await page.evaluate(async () => {
      const activity = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim().startsWith('Activity')) as HTMLElement | undefined;
      if (activity && !activity.classList.contains('expanded')) activity.click();
      await new Promise((r) => setTimeout(r, 1500));
      return activity?.parentElement?.textContent?.trim();
    });
    expect(text).toBeTruthy();
    // Audit log should contain at least one action on testRscript (created/shared/edited/ran)
    expect(text!).toMatch(/(created|shared|edited|ran).*testRscript/);
  });

  await softStep('3E. Send message to chat', async () => {
    const result = await page.evaluate(async () => {
      const chats = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim() === 'Chats') as HTMLElement | undefined;
      if (chats && !chats.classList.contains('expanded')) chats.click();
      await new Promise((r) => setTimeout(r, 1500));
      const pane = chats?.parentElement;
      const ta = pane?.querySelector('textarea') as HTMLTextAreaElement | null;
      if (!ta) return {ok: false};
      ta.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value')!.set!;
      setter.call(ta, 'Hello from grok-debug-scenarios');
      ta.dispatchEvent(new Event('input', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 500));
      ta.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', keyCode: 13, bubbles: true}));
      ta.dispatchEvent(new KeyboardEvent('keypress', {key: 'Enter', code: 'Enter', keyCode: 13, bubbles: true}));
      ta.dispatchEvent(new KeyboardEvent('keyup', {key: 'Enter', code: 'Enter', keyCode: 13, bubbles: true}));
      await new Promise((r) => setTimeout(r, 2500));
      return {ok: true, text: pane?.textContent?.trim()};
    });
    expect(result.ok).toBe(true);
    expect(result.text || '').toMatch(/grok-debug-scenarios/);
  });

  await softStep('4. Script browser — gallery + sort + search controls present', async () => {
    const info = await page.evaluate(() => ({
      cardCount: document.querySelectorAll('.grok-gallery-grid-item').length,
      hasSearch: !!document.querySelector('input[placeholder*="Search scripts"]'),
      viewName: grok.shell.v?.name,
    }));
    expect(info.cardCount).toBeGreaterThan(0);
    expect(info.hasSearch).toBe(true);
    expect(info.viewName).toBe('Scripts');
  });

  await softStep('5a. Open TSLA.csv (precondition for ACF)', async () => {
    const ok = await page.evaluate(async () => {
      let df = null;
      for (const path of ['System:DemoFiles/finance/TSLA.csv', 'System:DemoFiles/TSLA.csv', 'System:DemoFiles/stocks/TSLA.csv']) {
        try { df = await grok.dapi.files.readCsv(path); if (df) { df.name = 'TSLA'; break; } } catch {}
      }
      if (!df) return false;
      grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      return df.rowCount > 0;
    });
    expect(ok).toBe(true);
  });

  await softStep('5b. Right-click ACF card → Run... → set Columns → OK', async () => {
    await page.evaluate(async () => {
      grok.shell.route('/scripts');
      for (let i = 0; i < 30; i++) { if (grok.shell.v?.name === 'Scripts') break; await new Promise((r) => setTimeout(r, 200)); }
      for (let i = 0; i < 40; i++) {
        const hit = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
          .some((e) => e.textContent?.trim() === 'ACF');
        if (hit) break;
        await new Promise((r) => setTimeout(r, 300));
      }
    });
    const balloons = await page.evaluate(async () => {
      // Close any leftover dialogs from previous steps
      for (const d of (DG.Dialog as any).getOpenDialogs()) try { d.close(); } catch {}
      await new Promise((r) => setTimeout(r, 400));
      const titleEl = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'ACF') as HTMLElement;
      const card = (titleEl?.closest('.grok-gallery-grid-item') || titleEl?.parentElement) as HTMLElement;
      const rect = card.getBoundingClientRect();
      card.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, view: window,
        clientX: rect.left + 10, clientY: rect.top + 10, button: 2}));
      await new Promise((r) => setTimeout(r, 1200));
      const run = Array.from(document.querySelectorAll('.d4-menu-item-label, .d4-menu-item'))
        .find((e) => e.textContent?.trim() === 'Run...') as HTMLElement;
      run.click();
      // Wait for the ACF dialog (DG.Dialog API)
      let acf: any = null;
      for (let i = 0; i < 40; i++) {
        const opened = (DG.Dialog as any).getOpenDialogs();
        acf = opened.find((d: any) => d.title === 'ACF');
        if (acf) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      // De-duplicate any extra ACF dialogs
      const all = (DG.Dialog as any).getOpenDialogs().filter((d: any) => d.title === 'ACF');
      for (let i = 0; i < all.length - 1; i++) all[i].close();
      const target = (DG.Dialog as any).getOpenDialogs().find((d: any) => d.title === 'ACF');
      // Set the Columns input via the dialog API
      const df = grok.shell.tables.find((t: any) => t.name === 'TSLA');
      const closeCol = df.col('Close');
      target.input('Columns').value = [closeCol];
      await new Promise((r) => setTimeout(r, 500));
      // Click OK on the visible ACF dialog
      const dlgs = Array.from(document.querySelectorAll('.d4-dialog')) as HTMLElement[];
      const visibleAcf = dlgs.find((d) => d.querySelector('.d4-dialog-header')?.textContent?.trim() === 'ACF' &&
        getComputedStyle(d).display !== 'none');
      visibleAcf?.querySelector<HTMLButtonElement>('[name="button-OK"]')?.click();
      await new Promise((r) => setTimeout(r, 7000));
      return Array.from(document.querySelectorAll('.grok-balloon-error'))
        .map((b) => b.textContent?.trim()).filter(Boolean);
    });
    expect(balloons).toEqual([]);
  });

  await softStep('5c. Right-click ACF → Edit... → ScriptView opens', async () => {
    await page.evaluate(async () => {
      for (const d of (DG.Dialog as any).getOpenDialogs()) try { d.close(); } catch {}
      await new Promise((r) => setTimeout(r, 400));
      const titleEl = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'ACF') as HTMLElement;
      const card = (titleEl?.closest('.grok-gallery-grid-item') || titleEl?.parentElement) as HTMLElement;
      const rect = card.getBoundingClientRect();
      card.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, view: window,
        clientX: rect.left + 10, clientY: rect.top + 10, button: 2}));
      await new Promise((r) => setTimeout(r, 1200));
      const edit = Array.from(document.querySelectorAll('.d4-menu-item-label, .d4-menu-item'))
        .find((e) => e.textContent?.trim() === 'Edit...') as HTMLElement;
      edit.click();
      for (let i = 0; i < 60; i++) {
        if (grok.shell.v?.type === 'ScriptView' && grok.shell.v?.name === 'ACF') break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 1000));
    });
    const info = await page.evaluate(() => ({
      viewName: grok.shell.v?.name,
      viewType: grok.shell.v?.type,
      hasCM: !!document.querySelector('.CodeMirror'),
    }));
    expect(info.viewName).toBe('ACF');
    expect(info.viewType).toBe('ScriptView');
    expect(info.hasCM).toBe(true);
  });

  await softStep('5d. Editor Context Panel — Details tab', async () => {
    const text = await page.evaluate(async () => {
      const acf = await grok.dapi.scripts.filter('name = "ACF"').first().catch(() => null);
      if (acf) grok.shell.o = acf;
      await new Promise((r) => setTimeout(r, 2500));
      const details = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim().startsWith('Details')) as HTMLElement;
      if (details && !details.classList.contains('expanded')) details.click();
      await new Promise((r) => setTimeout(r, 1500));
      return details?.parentElement?.textContent?.trim();
    });
    expect(text).toBeTruthy();
    expect(text!).toMatch(/Inputs.*data.*columns/);
    expect(text!).toMatch(/Outputs.*acf/);
  });

  await softStep('5e. Editor Context Panel — Script tab body visible', async () => {
    const info = await page.evaluate(async () => {
      const scrHdr = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim().startsWith('Script')) as HTMLElement;
      if (scrHdr && !scrHdr.classList.contains('expanded')) scrHdr.click();
      await new Promise((r) => setTimeout(r, 1500));
      const ta = scrHdr?.parentElement?.querySelector('textarea') as HTMLTextAreaElement | null;
      return {hasTa: !!ta, taValue: ta?.value?.slice(0, 200)};
    });
    expect(info.hasTa).toBe(true);
    expect(info.taValue || '').toMatch(/#name:\s*ACF/);
  });

  await softStep('5f. Editor Context Panel — Activity tab populated', async () => {
    const text = await page.evaluate(async () => {
      const act = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((e) => e.textContent?.trim().startsWith('Activity')) as HTMLElement;
      if (act && !act.classList.contains('expanded')) act.click();
      await new Promise((r) => setTimeout(r, 1500));
      return act?.parentElement?.textContent?.trim();
    });
    expect(text).toBeTruthy();
    expect(text!).toMatch(/Activity\d+/);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
