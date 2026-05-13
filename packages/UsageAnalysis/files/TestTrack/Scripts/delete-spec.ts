import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Scripts Delete — testRscript', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  await softStep('1. Go to Browse > Platform > Functions > Scripts', async () => {
    await page.evaluate(() => { grok.shell.route('/scripts'); });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    const v = await page.evaluate(() => grok.shell.v?.name);
    expect(v).toBe('Scripts');
  });

  await softStep('2. Find the script via search', async () => {
    // Seed testRscript if missing, then route round-trip to refresh the gallery.
    await page.evaluate(async () => {
      const existing = await grok.dapi.scripts.filter('name = "testRscript"').first().catch(() => null);
      if (!existing) {
        const s = DG.Script.create(
          '#name: testRscript\n#language: r\n#input: dataframe table\n#output: int count\ncount <- nrow(table) * ncol(table)\n',
        );
        await grok.dapi.scripts.save(s);
      }
      grok.shell.route('/');
      await new Promise((r) => setTimeout(r, 400));
      grok.shell.route('/scripts');
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    const found = await page.evaluate(async () => {
      for (let i = 0; i < 40; i++) {
        const hit = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
          .some((e) => e.textContent?.trim() === 'testRscript');
        if (hit) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    });
    expect(found).toBe(true);
  });

  await softStep('3. Right-click the script, select Delete', async () => {
    const cardLocator = page.locator('.grok-gallery-grid-item')
      .filter({has: page.locator('.grok-gallery-grid-item-title', {hasText: /^testRscript$/})})
      .first();
    await cardLocator.scrollIntoViewIfNeeded();
    await cardLocator.click({button: 'right', force: true});
    const delLocator = page.locator('.d4-menu-item', {hasText: /^Delete$/}).first();
    await delLocator.waitFor({timeout: 15000});
    await delLocator.click({force: true});
  });

  await softStep('4. Click YES in the confirmation dialog', async () => {
    // .d4-dialog is position:fixed → offsetParent is always null even when
    // visible. Use bounding rect + computed style instead.
    const dialogHandled = await page.evaluate(async () => {
      const isVisible = (el: Element) => {
        const cs = getComputedStyle(el as HTMLElement);
        if (cs.display === 'none' || cs.visibility === 'hidden') return false;
        const r = el.getBoundingClientRect();
        return r.width > 0 && r.height > 0;
      };
      for (let i = 0; i < 60; i++) {
        const dialog = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => isVisible(d) && (d.textContent || '').includes('Delete script'));
        if (dialog) {
          const yes = dialog.querySelector('[name="button-YES"]') as HTMLElement | null;
          if (yes) { yes.click(); return true; }
          const byText = Array.from(dialog.querySelectorAll('button, .ui-btn, [role="button"]'))
            .find((b) => /^YES$/i.test((b.textContent || '').trim())) as HTMLElement | undefined;
          if (byText) { byText.click(); return true; }
        }
        await new Promise((r) => setTimeout(r, 100));
      }
      return false;
    });
    expect(dialogHandled).toBe(true);
  });

  await softStep('5. Check script is no longer present', async () => {
    const remaining = await page.evaluate(async () => {
      for (let i = 0; i < 30; i++) {
        const found = await grok.dapi.scripts.filter('name = "testRscript"').list();
        if (found.length === 0) return 0;
        await new Promise((r) => setTimeout(r, 500));
      }
      const found = await grok.dapi.scripts.filter('name = "testRscript"').list();
      return found.length;
    });
    expect(remaining).toBe(0);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
