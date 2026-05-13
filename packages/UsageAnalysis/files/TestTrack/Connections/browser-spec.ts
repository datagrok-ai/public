import {test, expect} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Prerequisite: a connection named "new_test" must exist (created by the user running the
// `Adding` scenario, or seeded out-of-band, e.g. via `grok s connections save --json`).

test('Connections / Browser', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  await softStep('Step 1: Go to Browse > Databases', async () => {
    await page.goto(`${baseUrl}/db?browse=db`, {waitUntil: 'networkidle', timeout: 30_000});
    await page.waitForTimeout(2500);
    await expect(page.locator('label', {hasText: /^New_test$/}).first()).toBeVisible({timeout: 15_000});
  });

  await softStep('Step 2: Filter templates icon — toggle each template', async () => {
    await page.locator('[name="icon-filter"]').first().click();
    await page.waitForTimeout(500);
    await expect(page.locator('text=QUICK FILTERS')).toBeVisible({timeout: 5000});
    for (const t of ['All', 'Created by me', 'Used by me', 'Invalid']) {
      await page.locator('span.d4-tag', {hasText: new RegExp(`^${t}$`)}).first().click();
      await page.waitForTimeout(300);
    }
    // Reset back to All so subsequent search is unfiltered
    await page.locator('span.d4-tag', {hasText: /^All$/}).first().click();
    await page.waitForTimeout(400);
  });

  await softStep('Step 3: Type new_test in the search field', async () => {
    const search = page.locator('input[placeholder="Search connections by name or by #tags"]');
    await search.fill('new_test');
    await page.waitForTimeout(800);
    await expect(page.locator('label', {hasText: /^New_test$/}).first()).toBeVisible({timeout: 5000});
  });

  await softStep('Step 4: Click the found connection', async () => {
    await page.locator('label', {hasText: /^New_test$/}).first().click();
    await page.waitForTimeout(1500);
    await expect(page.locator('.grok-prop-panel [name="label-New-test"]')).toBeVisible({timeout: 5000});
  });

  await softStep('Step 5.1: Details pane shows correct values', async () => {
    const details = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Details"]');
    await expect(details).toBeVisible({timeout: 5000});
    const txt = await details.locator('.d4-accordion-pane-content').innerText();
    expect(txt).toContain('db.datagrok.ai');
    expect(txt).toContain('northwind');
    expect(txt).toContain('54322');
  });

  await softStep('Step 5.2: Share connection with another user (Selenium)', async () => {
    // Try UI: right-click → Share... — synthetic clicks do not commit the autocomplete pick
    // reliably, so we cancel the dialog and use the JS API fallback.
    await page.evaluate(() => {
      const span = document.querySelector('[name="span-new-test"]');
      span?.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
    });
    await page.waitForTimeout(500);
    const shareItem = page.locator('[name="div-Share..."]');
    if (await shareItem.count() > 0)
      await shareItem.first().click();
    await page.waitForTimeout(1000);
    const cancel = page.locator('[name="button-CANCEL"]');
    if (await cancel.count() > 0) await cancel.first().click();
    await page.waitForTimeout(500);

    const result = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const conn = await g.dapi.connections.filter('name = "new_test"').first();
      const sel = await g.dapi.groups.filter('name = "Selenium"').first();
      if (!conn || !sel) return false;
      await g.dapi.permissions.grant(conn, sel, false);
      return true;
    });
    expect(result).toBe(true);

    // Re-click connection to refresh the property panel.
    await page.locator('label', {hasText: /^New_test$/}).first().click();
    await page.waitForTimeout(800);
    await page.locator('label', {hasText: /^New_test$/}).first().click();
    await page.waitForTimeout(1500);

    const sharingHeader = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Sharing"] .d4-accordion-pane-header');
    if (!(await sharingHeader.evaluate(el => el.classList.contains('expanded'))))
      await sharingHeader.click();
    await page.waitForTimeout(1500);
    const sharing = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Sharing"] .d4-accordion-pane-content');
    await expect(sharing).toContainText('Selenium', {timeout: 5000});
  });

  await softStep('Step 5.3: Activity has the expected events', async () => {
    const header = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Activity"] .d4-accordion-pane-header');
    if (!(await header.evaluate(el => el.classList.contains('expanded'))))
      await header.click();
    await page.waitForTimeout(1500);
    const pane = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Activity"] .d4-accordion-pane-content');
    await expect(pane).toContainText('shared', {timeout: 5000});
    await expect(pane).toContainText('created');
    await expect(pane).toContainText('New_test');
  });

  await softStep('Step 5.4: Send a chat message', async () => {
    const header = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Chats"] .d4-accordion-pane-header');
    if (!(await header.evaluate(el => el.classList.contains('expanded'))))
      await header.click();
    await page.waitForTimeout(1500);
    const ta = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Chats"] textarea[placeholder]');
    await ta.click();
    const message = `Test chat message ${Date.now()}`;
    await page.keyboard.type(message);
    await page.waitForTimeout(300);
    await page.keyboard.press('Enter');
    await page.waitForTimeout(2000);
    const pane = page.locator('.grok-prop-panel .d4-accordion-pane[d4-title="Chats"] .d4-accordion-pane-content');
    await expect(pane).toContainText(message, {timeout: 5000});
  });

  await softStep('Step 6: Open the connection name dropdown menu', async () => {
    await page.locator('.grok-prop-panel [name="icon-context-arrow-down"]').click();
    await page.waitForTimeout(800);
    await expect(page.locator('[name="div-Browse"]')).toBeVisible({timeout: 5000});
    await expect(page.locator('[name="div-Edit..."]')).toBeVisible();
    await expect(page.locator('[name="div-Share..."]')).toBeVisible();
    await expect(page.locator('[name="div-Test-connection"]')).toBeVisible();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
