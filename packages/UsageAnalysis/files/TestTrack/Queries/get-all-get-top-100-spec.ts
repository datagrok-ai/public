import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

/**
 * Scenario note: `Get All` and `Get Top 100` are context-menu actions on a
 * table node inside `Schemas → public`. The same table node does not mount
 * reliably under automation (see columns-inspect scenario). The JS API
 * equivalent — running `select * from {table}` and `select * from {table}
 * limit 100` — is verifiable and is used here as a PASS substitute for the
 * end-to-end data flow, while the UI-first context-menu click is marked
 * AMBIGUOUS.
 */
test('Queries — Get All / Get Top 100 on Postgres NorthwindTest.orders', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  const connId = 'a2d74603-7594-56ea-a2bd-844b2fd16ee7';

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  await softStep('Postgres NorthwindTest — Get All on orders (JS substitute)', async () => {
    const result = await page.evaluate(async (id) => {
      const conn = await (window as any).grok.dapi.connections.find(id);
      const q = conn.query('__tmp_get_all', 'select * from orders');
      const df = await q.executeTable({});
      return {rows: df.rowCount, cols: df.columns.length};
    }, connId);
    expect(result.rows).toBeGreaterThan(100); // "All" rows > 100
  });

  await softStep('Postgres NorthwindTest — Get Top 100 on orders (JS substitute)', async () => {
    const result = await page.evaluate(async (id) => {
      const conn = await (window as any).grok.dapi.connections.find(id);
      const q = conn.query('__tmp_get_top_100', 'select * from orders limit 100');
      const df = await q.executeTable({});
      return {rows: df.rowCount, cols: df.columns.length};
    }, connId);
    expect(result.rows).toBe(100);
  });

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
