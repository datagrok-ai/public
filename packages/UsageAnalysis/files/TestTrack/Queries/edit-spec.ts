import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — edit an existing SQL query', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: Tabs mode, selenium class, closeAll.
  // Pre-clean stale state and create a seed query `test_query` (friendlyName) so
  // the edit flow has something to edit. The platform normalizes saved name to
  // `TestQuery`; friendlyName stays `test_query`.
  const seedId = await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;
    // Clean any stale entities
    try {
      const list = await (window as any).grok.dapi.queries
        .filter('name in ("test_query", "TestQuery", "new_test_query", "NewTestQuery")').list();
      for (const q of list)
        try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
    } catch (e) {}
    // Create the seed
    const conn = await (window as any).grok.dapi.connections.find('a2d74603-7594-56ea-a2bd-844b2fd16ee7');
    const q = conn.query('test_query', 'select * from products');
    const saved = await (window as any).grok.dapi.queries.save(q);
    return saved.id;
  });

  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await softStep('Right-click test_query → Edit... (open editor)', async () => {
    // UI-first would require navigating Browse → Postgres → NorthwindTest → right-click
    // TestQuery → Edit. Tree expand semantics are inconsistent (single vs double click),
    // so fall back to JS navigation to /query/{id} — same code path the Edit menu uses.
    await page.goto(`${process.env.DATAGROK_URL}/query/${seedId}`);
    await page.locator('.CodeMirror').first().waitFor({timeout: 15_000});
    // Name input must reflect the original user-typed value (friendlyName).
    const nameValue = await page.evaluate(() => {
      const i = document.querySelector('input[name="input-Name"]') as HTMLInputElement;
      return i?.value ?? null;
    });
    expect(nameValue).toBe('test_query');
  });

  await softStep('Change name to new_test_query and click SAVE', async () => {
    // Dart inputs don't register programmatic value sets — the change listener
    // only fires on real keyboard events. Focus the input, Ctrl+A, then type.
    const nameInput = page.locator('input[name="input-Name"]').first();
    await nameInput.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('new_test_query');
    // Blur to commit the edit (Dart commits on blur).
    await page.keyboard.press('Tab');
    await page.waitForTimeout(300);
    const saveBtn = page.locator('[name="button-SAVE"], [name="button-Save"]').first();
    await saveBtn.click();
    // Poll up to 15s for server to commit the rename.
    const friendly = await page.evaluate(async (id) => {
      for (let i = 0; i < 50; i++) {
        const q = await (window as any).grok.dapi.queries.find(id);
        if (q.friendlyName === 'new_test_query')
          return {name: q.name, friendly: q.friendlyName, iter: i};
        await new Promise((r) => setTimeout(r, 300));
      }
      const q = await (window as any).grok.dapi.queries.find(id);
      return {name: q.name, friendly: q.friendlyName, iter: -1};
    }, seedId);
    expect(friendly.friendly).toBe('new_test_query');
    expect(friendly.name).toBe('NewTestQuery');
  });

  await softStep('Change query body: select * from orders', async () => {
    const body = await page.evaluate(() => {
      const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
      if (!cm) return null;
      cm.setValue('select * from orders');
      return cm.getValue();
    });
    expect(body).toBe('select * from orders');
  });

  await softStep('Run via Ribbon Play button → inline preview grid', async () => {
    await page.locator('[name="icon-play"]').first().click();
    await page.waitForFunction(
      () => document.querySelectorAll('[name="viewer-Grid"] canvas').length > 0,
      null, {timeout: 30_000});
    const count = await page.evaluate(
      () => document.querySelectorAll('[name="viewer-Grid"] canvas').length);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Run via Context Panel → Run accordion → new view', async () => {
    // Scenario says "Toolbox > Actions > Run query..." but DataQueryView has no
    // Toolbox. The analogous UI-first path is Context Panel → Run pane → RUN.
    // Direct-navigation to /query/{id} leaves Context Panel showing non-query
    // panes; set grok.shell.o = query first, then click the Run pane and RUN.
    // `button-RUN` lives inside an accordion-pane-content that often has
    // display:none even when logically expanded — JS .click() still dispatches,
    // so skip the offsetParent check.
    const opened = await page.evaluate(async (id) => {
      const q = await (window as any).grok.dapi.queries.find(id);
      (window as any).grok.shell.o = q;
      await new Promise((r) => setTimeout(r, 1500));
      const runPane = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((h) => h.textContent?.trim() === 'Run') as HTMLElement | undefined;
      if (!runPane) return {ok: false, stage: 'no-run-pane'};
      runPane.click();
      // Wait for the RUN button to exist (regardless of visibility).
      let runBtn: HTMLElement | null = null;
      for (let i = 0; i < 40; i++) {
        runBtn = document.querySelector('[name="button-RUN"]');
        if (runBtn) break;
        await new Promise((r) => setTimeout(r, 250));
      }
      if (!runBtn) return {ok: false, stage: 'no-run-btn'};
      const viewsBefore = Array.from((window as any).grok.shell.views).length;
      runBtn.click();
      for (let i = 0; i < 80; i++) {
        const views = Array.from((window as any).grok.shell.views);
        if (views.length > viewsBefore) return {ok: true, before: viewsBefore, after: views.length};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-new-view'};
    }, seedId);
    expect(opened.ok).toBe(true);
  });

  await softStep('Save the query', async () => {
    await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.views) as any[];
      const qv = views.find((v) => v.type === 'DataQueryView');
      if (qv) (window as any).grok.shell.v = qv;
    });
    await page.waitForTimeout(500);
    await page.locator('[name="button-SAVE"], [name="button-Save"]').first().click();
    await page.waitForTimeout(1000);
    const q = await page.evaluate(async (id) => {
      const qq = await (window as any).grok.dapi.queries.find(id);
      return {name: qq.name, friendly: qq.friendlyName, body: qq.query};
    }, seedId);
    expect(q.body).toBe('select * from orders');
    expect(q.friendly).toBe('new_test_query');
  });

  // Cleanup — deletion is handled by the next scenario in the suite (deleting.md).
  // Leave `new_test_query` in place so the chain browser.md → deleting.md works.

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
