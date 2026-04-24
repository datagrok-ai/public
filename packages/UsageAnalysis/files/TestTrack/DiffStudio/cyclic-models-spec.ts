import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('DiffStudio Cyclic Models (PK-PD): Load, Multiaxis+Facet, Count clickers, tooltips', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);
  await page.waitForTimeout(2000);

  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });

  await softStep('Step 1: Open Diff Studio + PK-PD from Library', async () => {
    await page.evaluate(async () => {
      const f = DG.Func.find({name: 'runDiffStudio', package: 'DiffStudio'})[0];
      const call = f.prepare();
      await call.call();
      grok.shell.addView(call.getOutputParamValue());
    });
    await page.waitForFunction(() => grok.shell.v?.name === 'Diff Studio', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const card = page.locator('.diff-studio-hub-card', {hasText: 'PK-PD'}).first();
    await card.waitFor({timeout: 15000});
    await card.dblclick();
    await page.waitForFunction(() => grok.shell.v?.name === 'PK-PD', null, {timeout: 30000});
    await page.waitForTimeout(2000);
    const countFound = await page.locator('[name="input-host-count"]').count();
    expect(countFound).toBeGreaterThan(0);
  });

  await softStep('Step 2: Check Multiaxis and Facet tabs are updated', async () => {
    const tabs = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.tab-handle')).map(t => t.textContent?.trim()));
    expect(tabs).toEqual(expect.arrayContaining(['Multiaxis', 'Facet']));
    const canvases = await page.evaluate(() => document.querySelectorAll('.d4-viewer canvas').length);
    expect(canvases).toBeGreaterThan(2);
  });

  await softStep('Step 3: Modify Count via clickers; real-time URL update', async () => {
    const before = await page.locator('[name="input-host-count"] input.ui-input-editor').inputValue();
    // Click [name="icon-plus"] inside the count host three times
    for (let i = 0; i < 3; i++) {
      await page.locator('[name="input-host-count"] [name="icon-plus"]').click();
      await page.waitForTimeout(400);
    }
    await page.waitForTimeout(1500);
    const after = await page.locator('[name="input-host-count"] input.ui-input-editor').inputValue();
    expect(parseInt(after)).toBeGreaterThan(parseInt(before));
    expect(page.url()).toContain('count=');
  });

  await softStep('Step 4: Tooltips on Begin, End, Step inputs', async () => {
    const tooltips = await page.evaluate(async () => {
      const names = ['begin', 'end', 'step'];
      const results: Record<string, string> = {};
      for (const n of names) {
        const host = document.querySelector(`[name="input-host-${n}"]`);
        const label = host?.querySelector('label') as HTMLElement;
        const rect = label?.getBoundingClientRect();
        label?.dispatchEvent(new MouseEvent('mouseover', {bubbles: true, clientX: rect.left + 5, clientY: rect.top + 5}));
        label?.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
        await new Promise(r => setTimeout(r, 800));
        results[n] = document.querySelector('.d4-tooltip')?.textContent?.trim() ?? '';
        label?.dispatchEvent(new MouseEvent('mouseleave', {bubbles: true}));
        await new Promise(r => setTimeout(r, 300));
      }
      return results;
    });
    expect(tooltips.begin.length).toBeGreaterThan(0);
    expect(tooltips.end.length).toBeGreaterThan(0);
    expect(tooltips.step.length).toBeGreaterThan(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
