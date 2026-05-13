import {test, expect, Page} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Self-contained: creates `new_test_postgres` and `test_postgres_2` via JS API if missing,
// then deletes them through the UI to exercise the right-click → Delete… → confirm flow.

const ensureConnection = async (page: Page, friendlyName: string, slug: string): Promise<string> => {
  return await page.evaluate(async ({friendlyName, slug}) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const list = await grok.dapi.connections.filter(`name = "${slug}"`).list();
    if (list.length > 0) {
      const conn = list[0];
      if (conn.friendlyName !== friendlyName) {
        conn.friendlyName = friendlyName;
        await grok.dapi.connections.save(conn);
      }
      return conn.nqName;
    }
    const fresh = DG.DataConnection.create(slug, {
      dataSource: 'Postgres', server: 'db.datagrok.ai', port: 54322, db: 'northwind', ssl: false,
    });
    fresh.name = slug;
    fresh.friendlyName = friendlyName;
    const saved = await grok.dapi.connections.save(fresh);
    return saved.nqName;
  }, {friendlyName, slug});
};

const deleteConnectionByDataLink = async (page: Page, dataLink: string) => {
  await page.evaluate(async (link) => {
    const target = document.querySelector(`[data-link="${link}"]`) as HTMLElement | null;
    if (!target) return;
    target.scrollIntoView({behavior: 'instant' as ScrollBehavior, block: 'center'});
    const r = target.getBoundingClientRect();
    target.dispatchEvent(new MouseEvent('contextmenu', {
      bubbles: true, cancelable: true, button: 2,
      clientX: r.left + 10, clientY: r.top + 5,
    }));
    await new Promise((res) => setTimeout(res, 500));
    const del = Array.from(document.querySelectorAll('.d4-menu-item'))
      .find((el) => el.textContent?.trim() === 'Delete...') as HTMLElement | undefined;
    del?.click();
    await new Promise((res) => setTimeout(res, 1000));
    const dialog = document.querySelector('.d4-dialog');
    const confirmBtn = (dialog?.querySelector('[name="button-DELETE"]')
      ?? Array.from(dialog?.querySelectorAll('button, .ui-btn') ?? [])
          .find((b) => b.textContent?.trim() === 'DELETE')) as HTMLElement | undefined;
    confirmBtn?.click();
    await new Promise((res) => setTimeout(res, 2500));
    (document.querySelector('[name="icon-sync"]') as HTMLElement | null)?.click();
    await new Promise((res) => setTimeout(res, 1500));
  }, dataLink);
};

test('Connections / Delete', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Setup: ensure both targets exist
  const ntpNq = await ensureConnection(page, 'new_test_postgres', 'NewTestPostgres');
  const tp2Nq = await ensureConnection(page, 'test_postgres_2', 'test_postgres_2');
  const ntpLink = '/db/' + ntpNq.replace(/^([^:]+):/, '$1.');
  const tp2Link = '/db/' + tp2Nq.replace(/^([^:]+):/, '$1.');

  await page.goto(`${baseUrl}/connect?browse=connections`, {waitUntil: 'networkidle', timeout: 30000});
  await page.waitForTimeout(3000);
  await page.evaluate(async () => {
    const postgresNode = Array.from(document.querySelectorAll('.d4-tree-view-group-label'))
      .find((el) => el.textContent?.trim() === 'Postgres');
    (postgresNode as HTMLElement | undefined)?.click();
    await new Promise((r) => setTimeout(r, 2000));
    (document.querySelector('[name="icon-sync"]') as HTMLElement | null)?.click();
    await new Promise((r) => setTimeout(r, 1500));
  });

  await softStep('Steps 1-2: Delete new_test_postgres', async () => {
    await expect(page.locator(`[data-link="${ntpLink}"]`)).toBeVisible({timeout: 10000});
    await deleteConnectionByDataLink(page, ntpLink);
    const present = await page.evaluate((link) => !!document.querySelector(`[data-link="${link}"]`), ntpLink);
    expect(present).toBe(false);
  });

  await softStep('Steps 3-4: Delete test_postgres_2', async () => {
    await expect(page.locator(`[data-link="${tp2Link}"]`)).toBeVisible({timeout: 10000});
    await deleteConnectionByDataLink(page, tp2Link);
    const present = await page.evaluate((link) => !!document.querySelector(`[data-link="${link}"]`), tp2Link);
    expect(present).toBe(false);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
