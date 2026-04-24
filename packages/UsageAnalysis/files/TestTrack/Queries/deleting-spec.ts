import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — delete new_test_query via context menu', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup + precondition: ensure new_test_query exists.
  const seedId = await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
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

  await softStep('Go to Browse → Databases → Postgres → NorthwindTest', async () => {
    await page.goto(`${process.env.DATAGROK_URL}/browse`);
    await page.locator('.d4-tree-view-root').waitFor({timeout: 15_000});
    const opened = await page.evaluate(async () => {
      const find = (label: string) =>
        Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label'))
          .find((el) => el.textContent?.trim() === label) as HTMLElement | undefined;
      const pg = find('Postgres');
      if (!pg) return {ok: false, stage: 'no-pg'};
      pg.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      for (let i = 0; i < 30; i++) {
        if (find('NorthwindTest')) break;
        await new Promise((r) => setTimeout(r, 300));
      }
      const nw = find('NorthwindTest');
      if (!nw) return {ok: false, stage: 'no-nw'};
      nw.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      for (let i = 0; i < 30; i++) {
        if ((window as any).grok.shell.v?.type === 'queries') return {ok: true};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-queries-view'};
    });
    expect(opened.ok).toBe(true);
  });

  await softStep('Right-click new_test_query → Delete → confirm DELETE', async () => {
    await page.waitForFunction(() => {
      return Array.from(document.querySelectorAll('*'))
        .some((el) => el.children.length === 0 && (el as HTMLElement).textContent?.trim() === 'new_test_query'
          && (el as HTMLElement).offsetParent !== null);
    }, null, {timeout: 15_000});
    // Dispatch contextmenu on the card containing the label.
    await page.evaluate(async () => {
      const label = Array.from(document.querySelectorAll('*'))
        .find((el) => el.children.length === 0 && (el as HTMLElement).textContent?.trim() === 'new_test_query'
          && (el as HTMLElement).offsetParent !== null) as HTMLElement;
      const card = label.closest('.grok-gallery-grid-item, .d4-gallery-item, .grok-entity, .grok-item, .d4-tile, .d4-card')
        ?? label.parentElement!;
      card.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
    });
    // Click Delete menu item.
    const deleteItem = page.locator('.d4-menu-item-label', {hasText: /^Delete$/}).first();
    await deleteItem.waitFor({timeout: 10_000});
    await deleteItem.click();
    // Confirmation dialog — click DELETE button.
    const deleteBtn = page.locator('.d4-dialog button', {hasText: 'DELETE'}).first();
    await deleteBtn.waitFor({timeout: 10_000});
    await deleteBtn.click();
    // Verify deletion.
    const gone = await page.evaluate(async (id) => {
      for (let i = 0; i < 40; i++) {
        const q = await (window as any).grok.dapi.queries.find(id).catch(() => null);
        if (!q) return true;
        await new Promise((r) => setTimeout(r, 300));
      }
      return false;
    }, seedId);
    expect(gone).toBe(true);
  });

  await softStep('Refresh Browse — verify query is no longer present', async () => {
    const visible = await page.evaluate(async () => {
      await new Promise((r) => setTimeout(r, 1000));
      return Array.from(document.querySelectorAll('*'))
        .some((el) => el.children.length === 0 && (el as HTMLElement).textContent?.trim() === 'new_test_query'
          && (el as HTMLElement).offsetParent !== null);
    });
    expect(visible).toBe(false);
  });

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
