import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — edit an existing SQL query', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup: Tabs mode, selenium class, closeAll. Pre-clean stale entities and
  // seed `test_query` on PostgresTest (Browse label `NorthwindTest`) so the
  // edit flow has something to edit. Server normalizes `name` to `TestQuery`;
  // `friendlyName` stays `test_query`.
  const seedId = await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;
    try {
      const list = await (window as any).grok.dapi.queries
        .filter('name in ("test_query", "TestQuery", "new_test_query", "NewTestQuery")').list();
      for (const q of list)
        try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
    } catch (e) {}
    const conn = await (window as any).grok.dapi.connections
      .find('a2d74603-7594-56ea-a2bd-844b2fd16ee7');
    const q = conn.query('test_query', 'select * from products');
    const saved = await (window as any).grok.dapi.queries.save(q);
    return saved.id;
  });

  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await softStep('Right-click test_query → Edit... (open editor)', async () => {
    // UI-first: open the connection's queries gallery via URL (the same page
    // shown when you double-click NorthwindTest in Browse), then right-click
    // the test_query card and pick Edit...
    await page.goto(`${process.env.DATAGROK_URL ?? 'http://localhost:8888'}/queries/Dbtests.PostgresTest?browse=db`);
    await page.locator('div[name="div-TestQuery"]').waitFor({timeout: 30_000});

    const opened = await page.evaluate(async () => {
      const card = document.querySelector('div[name="div-TestQuery"]') as HTMLElement | null;
      if (!card) return {ok: false, stage: 'no-card'};
      // contextmenu must include viewport coords for the menu items to render
      const rect = card.getBoundingClientRect();
      card.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2, buttons: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      let editItem: HTMLElement | undefined;
      for (let i = 0; i < 30; i++) {
        editItem = Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
        if (editItem) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      if (!editItem) return {ok: false, stage: 'no-edit-item'};
      editItem.click();
      return {ok: true};
    });
    expect(opened.ok).toBe(true);

    await page.locator('.CodeMirror').first().waitFor({timeout: 15_000});
    await page.locator('input[name="input-Name"]').waitFor({timeout: 15_000});
    const nameValue = await page.evaluate(() =>
      (document.querySelector('input[name="input-Name"]') as HTMLInputElement)?.value ?? null);
    expect(nameValue).toBe('test_query');
  });

  await softStep('Change name to new_test_query and click SAVE', async () => {
    // Dart inputs commit only on real keyboard events.
    const nameInput = page.locator('input[name="input-Name"]').first();
    await nameInput.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('new_test_query');
    await page.keyboard.press('Tab');
    await page.waitForTimeout(300);

    await page.locator('[name="button-Save"]').first().click();

    const persisted = await page.evaluate(async (id) => {
      for (let i = 0; i < 50; i++) {
        const q = await (window as any).grok.dapi.queries.find(id).catch(() => null);
        if (q && q.friendlyName === 'new_test_query')
          return {name: q.name, friendly: q.friendlyName};
        await new Promise((r) => setTimeout(r, 300));
      }
      const q = await (window as any).grok.dapi.queries.find(id);
      return {name: q.name, friendly: q.friendlyName};
    }, seedId);
    expect(persisted.friendly).toBe('new_test_query');
    expect(persisted.name).toBe('NewTestQuery');
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

  await softStep('Run via Toolbox → Actions → Run query... → new view', async () => {
    // DataQueryView's Toolbox includes `Actions` > `Run query...` — rendered
    // as `label.d4-link-action`. Clicking it opens a new TableView with the
    // saved query result and detaches the editor's CodeMirror from the DOM.
    // We must wait for the detach (`cmCount === 0`) before switching back —
    // otherwise the editor's pending body isn't yet cached and the switch-
    // back will reload server state, losing the unsaved `orders` edit.
    const opened = await page.evaluate(async () => {
      const viewsBefore = Array.from((window as any).grok.shell.views).length;
      const label = Array.from(document.querySelectorAll('label.d4-link-action'))
        .find((el) => el.textContent?.trim() === 'Run query...') as HTMLElement | undefined;
      if (!label) return {ok: false, stage: 'no-run-query'};
      label.click();
      let viewsAfter = viewsBefore;
      for (let i = 0; i < 80; i++) {
        viewsAfter = Array.from((window as any).grok.shell.views).length;
        if (viewsAfter > viewsBefore) break;
        await new Promise((r) => setTimeout(r, 250));
      }
      if (viewsAfter <= viewsBefore) return {ok: false, stage: 'no-new-view'};
      // Wait for the editor's CodeMirror to detach — this is the signal that
      // Dart has finished caching the DataQueryView's pending editor state.
      for (let i = 0; i < 40; i++) {
        if (document.querySelectorAll('.CodeMirror').length === 0) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      return {ok: true, viewsBefore, viewsAfter};
    });
    expect(opened.ok).toBe(true);
  });

  await softStep('Save the query', async () => {
    // Switch back to the DataQueryView, wait for the editor to remount with
    // the cached body (`select * from orders`), then click SAVE. The single
    // evaluate matches the MCP timing and avoids Playwright's actionability
    // delay between view switch and click.
    const q = await page.evaluate(async (id) => {
      const views = Array.from((window as any).grok.shell.views) as any[];
      const qv = views.find((v) => v.type === 'DataQueryView');
      if (qv) (window as any).grok.shell.v = qv;
      // Wait for the editor's CodeMirror to remount with the cached body.
      let cm: any = null;
      for (let i = 0; i < 50; i++) {
        cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
        if (cm && cm.getValue() === 'select * from orders') break;
        await new Promise((r) => setTimeout(r, 200));
      }
      if (!cm) return {ok: false, stage: 'no-cm-after-switch'};
      // Defensive: if remount restored the server-side body, re-apply ours.
      if (cm.getValue() !== 'select * from orders') {
        cm.setValue('select * from orders');
        await new Promise((r) => setTimeout(r, 300));
      }
      const saveBtn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      if (!saveBtn) return {ok: false, stage: 'no-save-btn'};
      saveBtn.click();
      for (let i = 0; i < 50; i++) {
        const qq = await (window as any).grok.dapi.queries.find(id);
        if (qq.query === 'select * from orders' && qq.friendlyName === 'new_test_query')
          return {ok: true, name: qq.name, friendly: qq.friendlyName, body: qq.query};
        await new Promise((r) => setTimeout(r, 300));
      }
      const qq = await (window as any).grok.dapi.queries.find(id);
      return {ok: false, name: qq.name, friendly: qq.friendlyName, body: qq.query};
    }, seedId);
    expect(q.body).toBe('select * from orders');
    expect(q.friendly).toBe('new_test_query');
  });

  // Cleanup — leave `new_test_query` in place so the chain → deleting.md
  // continues to work.

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
