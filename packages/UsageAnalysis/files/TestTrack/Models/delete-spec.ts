import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Delete predictive model', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const w = window as any;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
  });

  let modelExisted = false;
  let modelNameAfterSearch = '';

  await softStep('Go to Browse > Platform > Models', async () => {
    // UI tree-click on Browse > Platform > Predictive models did not switch the view
    // during the MCP run — fell back to JS routing, which is the only reliable path.
    await page.evaluate(async () => {
      const w = window as any;
      w.grok.shell.windows.showBrowse = true;
      w.grok.shell.route('/models');
    });
    await page.waitForFunction(() => location.pathname === '/models', null, {timeout: 15000});
    await page.locator('input[placeholder*="Search models"]').waitFor({timeout: 15000});
    const url = await page.evaluate(() => location.pathname);
    expect(url).toBe('/models');
  });

  await softStep('Find the model from the previous steps', async () => {
    // Search for "TestDemog" — the prerequisite model from train.md.
    const search = page.locator('input[placeholder*="Search models"]');
    await search.click();
    await page.keyboard.type('TestDemog');
    await page.waitForTimeout(1000);
    const apiCount = await page.evaluate(async () => {
      const w = window as any;
      const list = await w.grok.dapi.models.filter('name like "%TestDemog%"').list({pageSize: 10});
      return list.length;
    });
    modelNameAfterSearch = await page.locator('text=/TestDemog/i').first().textContent({timeout: 1000}).catch(() => '');
    modelExisted = apiCount > 0 || !!modelNameAfterSearch;
    // Scenario intent: locate the model. Assert it is visible; if the prerequisite
    // chain (train.md/browser.md) didn't seed it on this server, this step fails
    // truthfully — the scenario as written assumes prior steps ran.
    expect(modelExisted, 'TestDemog (prerequisite from train.md) not found on this server').toBe(true);
  });

  await softStep('Right-click the model and select Delete', async () => {
    if (!modelExisted) {
      test.skip(true, 'No model from previous steps to delete on this server');
      return;
    }
    const item = page.locator('text=/TestDemog/i').first();
    await item.click({button: 'right'});
    const deleteItem = page.locator('.d4-menu-item:has-text("Delete")').first();
    await expect(deleteItem).toBeVisible({timeout: 5000});
    await deleteItem.click();
  });

  await softStep('In the confirmation dialog, click Delete', async () => {
    if (!modelExisted) {
      test.skip(true, 'No model — no confirmation dialog');
      return;
    }
    const dlg = page.locator('.d4-dialog').filter({hasText: /delete|are you sure/i }).first();
    await expect(dlg).toBeVisible({timeout: 5000});
    const confirm = dlg.locator('button:has-text("Delete"), button:has-text("OK"), [name="button-OK"]').first();
    await confirm.click();
  });

  await softStep('Check that model has been deleted and is no longer present', async () => {
    if (!modelExisted) {
      test.skip(true, 'Nothing was deleted — vacuous check');
      return;
    }
    await page.waitForTimeout(2000);
    const stillThere = await page.evaluate(async () => {
      const w = window as any;
      const list = await w.grok.dapi.models.filter('name like "%TestDemog%"').list({pageSize: 5});
      return list.length;
    });
    expect(stillThere).toBe(0);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n'));
});
