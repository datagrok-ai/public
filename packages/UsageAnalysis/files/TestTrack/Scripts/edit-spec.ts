import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Scripts Edit — testRscript', async ({page}) => {
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
    await page.waitForFunction(() => !!document.querySelector('[name="button-New"]'),
      null, {timeout: 15000});
    const v = await page.evaluate(() => grok.shell.v?.name);
    expect(v).toBe('Scripts');
  });

  // Prerequisite — if testRscript was deleted by a previous run, seed it
  await softStep('[pre] Ensure testRscript exists', async () => {
    await page.evaluate(async () => {
      const existing = await grok.dapi.scripts.filter('name = "testRscript"').first().catch(() => null);
      if (!existing) {
        const body = '#name: testRscript\n#language: r\n#sample: cars.csv\n#input: dataframe table [Data table]\n#output: int count [Number of cells in table]\n#output: string newParam\ncount <- nrow(table) * ncol(table)\n';
        const s = DG.Script.create(body);
        await grok.dapi.scripts.save(s);
      }
      // Full route round-trip so the gallery re-fetches
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

  await softStep('2. Find testRscript and double-click', async () => {
    await page.evaluate(async () => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'testRscript') as HTMLElement;
      const card = (label?.closest('.grok-gallery-grid-item, .grok-entity-card') || label) as HTMLElement;
      card.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, view: window}));
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.type === 'ScriptView',
      null, {timeout: 30000});
    await page.waitForFunction(() => !!(document.querySelector('.CodeMirror') as any)?.CodeMirror,
      null, {timeout: 15000});
    const view = await page.evaluate(() => ({name: grok.shell.v?.name, type: grok.shell.v?.type}));
    expect(view.name).toBe('testRscript');
    expect(view.type).toBe('ScriptView');
  });

  await softStep('3. Add newParam="test" to script body', async () => {
    await page.evaluate(() => {
      const cm = (document.querySelector('.CodeMirror') as any).CodeMirror;
      const current = cm.getValue();
      if (!current.includes('newParam="test"'))
        cm.setValue(current + '\nnewParam="test"\n');
    });
    const code = await page.evaluate(() => (document.querySelector('.CodeMirror') as any)?.CodeMirror?.getValue());
    expect(code).toContain('newParam="test"');
  });

  await softStep('4. Click Save', async () => {
    await page.evaluate(() => {
      (document.querySelector('[name="button-Save"]') as HTMLElement)?.click();
    });
    // Poll the server until the saved body reflects the edit
    const saved = await page.evaluate(async () => {
      for (let i = 0; i < 30; i++) {
        const s = await grok.dapi.scripts.filter('name = "testRscript"').first().catch(() => null);
        if (s?.script?.includes('newParam="test"')) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    });
    expect(saved).toBe(true);
  });

  await softStep('5. Close script view', async () => {
    await page.evaluate(() => { grok.shell.v?.close(); });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    const v = await page.evaluate(() => grok.shell.v?.name);
    expect(v).toBe('Scripts');
  });

  await softStep('6. Double-click testRscript again; verify newParam="test"', async () => {
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
      .some((e) => e.textContent?.trim() === 'testRscript'),
      null, {timeout: 30000});
    await page.evaluate(() => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'testRscript') as HTMLElement;
      const card = (label?.closest('.grok-gallery-grid-item, .grok-entity-card') || label) as HTMLElement;
      card?.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, view: window}));
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'testRscript',
      null, {timeout: 30000});
    await page.waitForFunction(() => !!(document.querySelector('.CodeMirror') as any)?.CodeMirror,
      null, {timeout: 15000});
    const info = await page.evaluate(() => {
      const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
      return {view: grok.shell.v?.name, code: cm?.getValue()};
    });
    expect(info.view).toBe('testRscript');
    expect(info.code).toContain('newParam="test"');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
