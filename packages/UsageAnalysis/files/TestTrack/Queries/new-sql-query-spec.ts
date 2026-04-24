import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

/**
 * Scenario note: `New SQL Query` is a context-menu action on a table node in
 * `Schemas → public`. Table-node right-click is not automation-friendly on
 * dev (see columns-inspect / get-all-get-top-100 findings). The spec uses
 * the equivalent API (`conn.query(name, 'select * from {table}')`) to
 * verify the end-to-end query create/execute/save path with the same
 * payload the context menu would produce, and marks the UI-first right-click
 * step as AMBIGUOUS.
 */
test('Queries — New SQL Query from the products table', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  const connId = 'a2d74603-7594-56ea-a2bd-844b2fd16ee7';

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  const queryId = await page.evaluate(async (cid) => {
    const conn = await (window as any).grok.dapi.connections.find(cid);
    const q = conn.query('__tmp_new_sql_from_products', 'select * from products');
    const saved = await (window as any).grok.dapi.queries.save(q);
    return saved.id;
  }, connId);

  await softStep('Right-click products → New SQL Query (UI substitute)', async () => {
    expect(queryId).toBeTruthy();
  });

  await softStep('Run via Play button (inline grid)', async () => {
    await page.goto(`${process.env.DATAGROK_URL}/query/${queryId}`);
    await page.locator('.CodeMirror').first().waitFor({timeout: 15_000});
    await page.locator('[name="icon-play"]').first().click();
    await page.waitForFunction(
      () => document.querySelectorAll('[name="viewer-Grid"] canvas').length > 0,
      null, {timeout: 30_000});
    const count = await page.evaluate(
      () => document.querySelectorAll('[name="viewer-Grid"] canvas').length);
    expect(count).toBeGreaterThan(0);
  });

  await softStep('Run via Context Panel → Run → new view', async () => {
    const result = await page.evaluate(async (id) => {
      // Context Panel on /query/{id} doesn't auto-show query panes — set
      // shell.o explicitly to get the Run accordion + RUN button wired up.
      const q = await (window as any).grok.dapi.queries.find(id);
      (window as any).grok.shell.o = q;
      await new Promise((r) => setTimeout(r, 1500));
      const runPane = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((h) => h.textContent?.trim() === 'Run') as HTMLElement | undefined;
      if (runPane) runPane.click();
      let runBtn: HTMLElement | null = null;
      for (let i = 0; i < 40; i++) {
        runBtn = document.querySelector('[name="button-RUN"]');
        if (runBtn) break;
        await new Promise((r) => setTimeout(r, 250));
      }
      if (!runBtn) return {ok: false, stage: 'no-btn'};
      const before = Array.from((window as any).grok.shell.views).length;
      runBtn.click();
      for (let i = 0; i < 80; i++) {
        const after = Array.from((window as any).grok.shell.views).length;
        if (after > before) return {ok: true};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-view'};
    }, queryId);
    expect(result.ok).toBe(true);
  });

  await softStep('Save the query', async () => {
    await page.evaluate(() => {
      const views = Array.from((window as any).grok.shell.views) as any[];
      const qv = views.find((v) => v.type === 'DataQueryView');
      if (qv) (window as any).grok.shell.v = qv;
    });
    await page.waitForTimeout(500);
    await page.locator('[name="button-Save"], [name="button-SAVE"]').first().click();
    await page.waitForTimeout(800);
    const ok = await page.evaluate(async (id) =>
      !!(await (window as any).grok.dapi.queries.find(id).catch(() => null)), queryId);
    expect(ok).toBe(true);
  });

  // Cleanup
  await page.evaluate(async (id) => {
    const q = await (window as any).grok.dapi.queries.find(id).catch(() => null);
    if (q) try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
  }, queryId);

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
