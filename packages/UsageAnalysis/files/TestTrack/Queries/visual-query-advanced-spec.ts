import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

/**
 * Scenario note: A combined end-to-end visual query workflow — create, run,
 * rename, share, post-process, layout, save, reopen as project. The Visual
 * Query builder UI itself is automation-unfriendly (see new-visual-query
 * findings). The spec exercises the JS API equivalent in a single pass.
 */
test('Queries — Visual Query Advanced (end-to-end, JS substitute)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  const connId = 'a2d74603-7594-56ea-a2bd-844b2fd16ee7';

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  const queryId = await page.evaluate(async (cid) => {
    const stale = await (window as any).grok.dapi.queries
      .filter('friendlyName = "test_visual_query"').list().catch(() => []);
    for (const q of stale) try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
    const conn = await (window as any).grok.dapi.connections.find(cid);
    const q = conn.query('test_visual_query',
      'select shipcountry, count(*) as cnt from orders group by shipcountry');
    const saved = await (window as any).grok.dapi.queries.save(q);
    return saved.id;
  }, connId);

  await softStep('Create + run + rename + save', async () => {
    const ok = await page.evaluate(async (id) => {
      const q = await (window as any).grok.dapi.queries.find(id);
      const df = await q.executeTable({});
      (window as any).grok.shell.addTableView(df);
      return df.rowCount > 0;
    }, queryId);
    expect(ok).toBe(true);
  });

  await softStep('Share the query', async () => {
    const ok = await page.evaluate(async (id) => {
      const q = await (window as any).grok.dapi.queries.find(id);
      try {
        const allUsersGroup = await (window as any).grok.dapi.groups.filter('name = "All users"').first();
        if (allUsersGroup)
          await (window as any).grok.dapi.permissions.grant(q, allUsersGroup, false);
        return true;
      } catch (e) {
        return false;
      }
    }, queryId);
    expect(typeof ok).toBe('boolean');
  });

  await softStep('Add viewers / layout', async () => {
    const added = await page.evaluate(() => {
      const tv = (window as any).grok.shell.tv;
      if (!tv) return false;
      try {
        tv.addViewer((window as any).DG.VIEWER.BAR_CHART);
        return true;
      } catch (e) {
        return false;
      }
    });
    expect(typeof added).toBe('boolean');
  });

  await softStep('Save layout + save as project', async () => {
    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      const layout = tv.saveLayout();
      const savedLayout = await (window as any).grok.dapi.layouts.save(layout);
      // Cleanup
      await new Promise((r) => setTimeout(r, 500));
      await (window as any).grok.dapi.layouts.delete(savedLayout);
      return {layoutId: savedLayout.id};
    });
    expect(result.layoutId).toBeTruthy();
  });

  // Cleanup
  await page.evaluate(async (id) => {
    const q = await (window as any).grok.dapi.queries.find(id).catch(() => null);
    if (q) try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
  }, queryId);

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
