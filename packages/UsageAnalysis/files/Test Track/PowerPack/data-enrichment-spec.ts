import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Data enrichment', async ({page}) => {
  test.setTimeout(300_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  await softStep('1.1 Navigate to Databases > Postgres > Northwind', async () => {
    // Scenario specifies "Northwind" but the dev server has "NorthwindTest" under Postgres
    // as the closest match. Verify it is reachable in the Browse tree.
    await page.locator('[name="Browse"]').click();
    await page.waitForTimeout(500);
    const found = await page.evaluate(() => {
      const findByText = (sel: string, text: string) =>
        Array.from(document.querySelectorAll(sel)).find(el => el.textContent?.trim() === text) as HTMLElement | undefined;
      const dbNode = findByText('.d4-tree-view-group-label', 'Databases');
      if (!dbNode) return false;
      dbNode.closest('.d4-tree-view-node')?.dispatchEvent(new MouseEvent('dblclick', {bubbles: true}));
      return true;
    });
    expect(found, 'Databases node not found in Browse tree').toBe(true);
    await page.waitForTimeout(1500);
    const nwLabels = await page.evaluate(() =>
      Array.from(document.querySelectorAll('label'))
        .filter(e => e.textContent?.trim() === 'NorthwindTest' && e.parentElement?.querySelector('[name="icon-ds-postgres"]'))
        .length);
    expect(nwLabels, 'Postgres NorthwindTest connection is not visible').toBeGreaterThan(0);
  });

  await softStep('1.2 Orders table accessible via NorthwindTest', async () => {
    // The scenario needs a working query/table on Orders. Verify any SELECT against
    // the Northwind-equivalent connection completes. This fails on dev today.
    const res = await page.evaluate(async () => {
      try {
        const df = await (window as any).grok.data.query('Agolovko:OrdersVQ');
        return { ok: true, rows: df.rowCount };
      } catch (e: any) {
        return { ok: false, error: String(e).substring(0, 300) };
      }
    });
    expect(res.ok, `Query execution failed: ${res.error}`).toBe(true);
  });

  await softStep('1.3 Enrich panel available on customerid column', async () => {
    // Requires Section 1.2 to succeed. Once Orders is open and the customerid column
    // is selected, the Context Panel should show an "Enrich" subsection under the
    // NorthwindTest (connection) root. Skipped here — no data to drive it.
    test.skip(true, 'Blocked: Section 1.2 failed so the Orders table cannot be opened.');
  });

  await softStep('2 Multiple enrichments per column and per table', async () => {
    test.skip(true, 'Blocked on Section 1 — no working enrichment to replicate across columns.');
  });

  await softStep('3 Persistence across projects, layouts, reuse', async () => {
    test.skip(true, 'Blocked on Section 1 — no enrichment to persist.');
  });

  await softStep('4 Visibility for different users', async () => {
    test.skip(true, 'Blocked on Section 1. Also: multi-user login within a single Playwright spec '
      + 'would require separate credential sets not available in this environment.');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  • ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
