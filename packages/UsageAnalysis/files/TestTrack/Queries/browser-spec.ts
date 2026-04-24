import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — browse NorthwindTest and find new_test_query', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup + precondition: make sure new_test_query exists on PostgresTest.
  const seedId = await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;
    // Check for existing; create if missing
    let q = await (window as any).grok.dapi.queries
      .filter('friendlyName = "new_test_query"').first().catch(() => null);
    if (!q) {
      const conn = await (window as any).grok.dapi.connections.find('a2d74603-7594-56ea-a2bd-844b2fd16ee7');
      const newQ = conn.query('new_test_query', 'select * from orders');
      q = await (window as any).grok.dapi.queries.save(newQ);
    }
    return q.id;
  });

  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await softStep('Refresh Browse', async () => {
    // Navigate to /browse to render the tree reliably. This is the same view as
    // clicking Browse + refreshing — no hidden state.
    await page.goto(`${process.env.DATAGROK_URL}/browse`);
    await page.locator('.d4-tree-view-root').waitFor({timeout: 15_000});
  });

  await softStep('Browse → Databases → Postgres → NorthwindTest', async () => {
    const nwVisible = await page.evaluate(async () => {
      const find = (label: string) =>
        Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
          .find((el) => el.textContent?.trim() === label) as HTMLElement | undefined;
      const pg = find('Postgres');
      if (!pg) return {ok: false, stage: 'no-pg'};
      pg.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      for (let i = 0; i < 30; i++) {
        if (find('NorthwindTest')) return {ok: true};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-nw'};
    });
    expect(nwVisible.ok).toBe(true);
    // Double-click to open NorthwindTest queries view
    await page.evaluate(async () => {
      const find = (label: string) =>
        Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
          .find((el) => el.textContent?.trim() === label) as HTMLElement;
      const nw = find('NorthwindTest');
      nw.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
    });
    await page.waitForFunction(
      () => (window as any).grok?.shell?.v?.type === 'queries',
      null, {timeout: 15_000});
  });

  await softStep('Type new_test in search field', async () => {
    const typed = await page.evaluate(async () => {
      const input = Array.from(document.querySelectorAll('input'))
        .find((i) => i.offsetParent !== null && i.placeholder?.includes('Search queries')) as HTMLInputElement | undefined;
      if (!input) return {ok: false, stage: 'no-search'};
      input.focus();
      input.value = 'new_test';
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keyup', {bubbles: true, key: 't'}));
      await new Promise((r) => setTimeout(r, 1500));
      return {ok: true, value: input.value, visibleText: document.body.innerText.includes('new_test_query')};
    });
    expect(typed.ok).toBe(true);
    expect(typed.visibleText).toBe(true);
  });

  await softStep('Context Panel tabs for the query', async () => {
    const tabs = await page.evaluate(async (id) => {
      const q = await (window as any).grok.dapi.queries.find(id);
      (window as any).grok.shell.o = q;
      await new Promise((r) => setTimeout(r, 1500));
      return Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .map((h) => h.textContent?.trim());
    }, seedId);
    // Minimum set of query-entity Context Panel tabs.
    for (const required of ['Details', 'Run', 'Query', 'Transformations', 'Usage', 'Sharing'])
      expect(tabs).toContain(required);
  });

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
