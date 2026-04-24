import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — adding a new SQL query', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: Tabs mode, selenium class, closeAll, show Browse.
  // Pre-clean: delete any stale test_query left from previous runs.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;
    try {
      const existing = await (window as any).grok.dapi.queries
        .filter('name in ("test_query", "TestQuery")').list();
      for (const q of existing) {
        try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
      }
    } catch (e) {}
  });

  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await softStep('Go to Browse → Databases → Postgres', async () => {
    // Expand Databases group, then Postgres, wait for NorthwindTest to appear.
    const expanded = await page.evaluate(async () => {
      const db = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'Databases') as HTMLElement | undefined;
      if (!db) return {ok: false, stage: 'Databases'};
      db.click();
      // Wait for Postgres
      for (let i = 0; i < 30; i++) {
        const pg = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
          .find((el) => el.textContent?.trim() === 'Postgres') as HTMLElement | undefined;
        if (pg) {
          pg.click();
          pg.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
          break;
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      for (let i = 0; i < 30; i++) {
        const nw = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
          .find((el) => el.textContent?.trim() === 'NorthwindTest');
        if (nw) return {ok: true, stage: 'NorthwindTest-visible'};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'NorthwindTest'};
    });
    expect(expanded.ok).toBe(true);
  });

  await softStep('Right-click NorthwindTest → New Query...', async () => {
    const opened = await page.evaluate(async () => {
      const nw = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
        .find((el) => el.textContent?.trim() === 'NorthwindTest') as HTMLElement | undefined;
      if (!nw) return {ok: false, stage: 'NorthwindTest-missing'};
      const node = nw.closest('.d4-tree-view-node');
      node!.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
      for (let i = 0; i < 20; i++) {
        const item = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'New Query...') as HTMLElement | undefined;
        if (item) {
          item.click();
          return {ok: true};
        }
        await new Promise((r) => setTimeout(r, 200));
      }
      return {ok: false, stage: 'menu-item-missing'};
    });
    expect(opened.ok).toBe(true);
    // Wait for the editor (CodeMirror) to mount.
    await page.locator('.CodeMirror').first().waitFor({timeout: 15_000});
  });

  await softStep('Enter test_query into the Name field', async () => {
    const nameInput = page.locator('input[type="text"]:visible').filter({hasText: ''}).first();
    // Find the name input — it's the visible input that is not the console / spotlight.
    const nameHandle = await page.evaluateHandle(() => {
      const inputs = Array.from(document.querySelectorAll('input[type="text"]')) as HTMLInputElement[];
      return inputs.find((i) =>
        i.offsetParent !== null &&
        i.placeholder !== 'Alt+Q to search commands' &&
        i.placeholder !== '> Enter command');
    });
    await nameHandle.asElement()!.focus();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('test_query');
    const value = await page.evaluate(() => {
      const inputs = Array.from(document.querySelectorAll('input[type="text"]')) as HTMLInputElement[];
      const input = inputs.find((i) =>
        i.offsetParent !== null &&
        i.placeholder !== 'Alt+Q to search commands' &&
        i.placeholder !== '> Enter command');
      return input?.value ?? '';
    });
    expect(value).toBe('test_query');
    // Silence unused-locator lint: nameInput declared for future UI-first upgrade.
    void nameInput;
  });

  await softStep('Enter query body: select * from products', async () => {
    const body = await page.evaluate(() => {
      const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
      if (!cm) return null;
      cm.setValue('select * from products');
      return cm.getValue();
    });
    expect(body).toBe('select * from products');
  });

  await softStep('Run via Ribbon Play button → inline preview grid', async () => {
    const play = page.locator('[name="icon-play"]').first();
    await play.click();
    // Expect at least one grid canvas to render inside the current view.
    await page.waitForFunction(
      () => document.querySelectorAll('[name="viewer-Grid"] canvas').length >= 1,
      null, {timeout: 30_000});
    const count = await page.evaluate(
      () => document.querySelectorAll('[name="viewer-Grid"] canvas').length);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Run via Toolbox → Actions → Run query... → new view', async () => {
    const opened = await page.evaluate(async () => {
      const viewsBefore = Array.from((window as any).grok.shell.views).length;
      const label = Array.from(document.querySelectorAll('label, .d4-toolbox-item, a, .d4-link-label'))
        .find((el) => (el.textContent?.trim() ?? '').startsWith('Run query')) as HTMLElement | undefined;
      if (!label) return {ok: false, stage: 'run-query-missing'};
      label.click();
      for (let i = 0; i < 80; i++) {
        const views = Array.from((window as any).grok.shell.views) as any[];
        if (views.length > viewsBefore || views.some((v) => v.type === 'TableView' && v.name === 'test_query'))
          return {ok: true, viewsBefore, viewsAfter: views.length};
        await new Promise((r) => setTimeout(r, 250));
      }
      return {ok: false, stage: 'new-view-missing'};
    });
    expect(opened.ok).toBe(true);
  });

  await softStep('Save the query', async () => {
    // Switch back to the query editor view, then click SAVE.
    await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.views) as any[];
      const qv = views.find((v) => v.type === 'DataQueryView');
      if (qv) (window as any).grok.shell.v = qv;
    });
    await page.waitForTimeout(500);
    const saveBtn = page.locator('[name="button-SAVE"], [name="button-Save"]').first();
    await saveBtn.click();
    // Verify the query is persisted (connection name on dev is "PostgresTest", the
    // NorthwindTest node's underlying connection; accept either name).
    const saved = await page.evaluate(async () => {
      for (let i = 0; i < 40; i++) {
        const q = await (window as any).grok.dapi.queries
          .filter('name in ("test_query", "TestQuery")').first().catch(() => null);
        if (q && q.connection) return {ok: true, connName: q.connection.name, id: q.id};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false};
    });
    expect(saved.ok).toBe(true);
  });

  // Cleanup
  await page.evaluate(async () => {
    try {
      const all = await (window as any).grok.dapi.queries
        .filter('name in ("test_query", "TestQuery")').list();
      for (const q of all) {
        try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
      }
    } catch (e) {}
    (window as any).grok.shell.closeAll();
  });

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
