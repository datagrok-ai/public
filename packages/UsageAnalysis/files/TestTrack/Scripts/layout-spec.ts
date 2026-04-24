import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Scripts Layout — test_Layout with Layout tab', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  await softStep('1. Create test_Layout JavaScript script', async () => {
    await page.evaluate(async () => {
      const prior = await grok.dapi.scripts.filter('name = "test_Layout"').list();
      for (const s of prior) await grok.dapi.scripts.delete(s);
      const body = `//name: test_Layout
//language: javascript
//input: int idx=1
//output: dataframe df

const fileList = await grok.dapi.files.list('System:DemoFiles/chem', true, '');
const csvFiles = fileList.filter((fi) => fi.fileName.endsWith('.csv'));
if (csvFiles.length === 0)
  throw new Error('No CSV files found in System:DemoFiles/chem');

csvFiles.sort((a, b) => b.createdOn - a.createdOn);
const lastModifiedFile = csvFiles[idx];
const csv = await grok.dapi.files.readAsText(lastModifiedFile.fullPath);
df = DG.DataFrame.fromCsv(csv);
`;
      const s = DG.Script.create(body);
      await grok.dapi.scripts.save(s);
    });
    const count = await page.evaluate(async () =>
      (await grok.dapi.scripts.filter('name = "test_Layout"').list()).length);
    expect(count).toBe(1);
  });

  await softStep('2. Open script editor, go to Layout tab', async () => {
    await page.evaluate(() => { grok.shell.route('/scripts'); });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'Scripts',
      null, {timeout: 30000});
    // The gallery doesn't auto-subscribe to `scripts.save`, so poll until the card appears
    await page.waitForFunction(() => Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
      .some((e) => e.textContent?.trim() === 'test_Layout'),
      null, {timeout: 30000});
    await page.evaluate(() => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find((e) => e.textContent?.trim() === 'test_Layout') as HTMLElement;
      const card = label?.closest('.grok-gallery-grid-item') as HTMLElement;
      card?.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, view: window}));
    });
    await page.waitForFunction(() => (window as any).grok?.shell?.v?.name === 'test_Layout',
      null, {timeout: 30000});
    await page.evaluate(async () => {
      const layout = Array.from(document.querySelectorAll('.d4-tab-header'))
        .find((e) => e.textContent?.trim() === 'Layout') as HTMLElement;
      layout?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const vName = await page.evaluate(() => grok.shell.v?.name);
    expect(vName).toBe('test_Layout');
  });

  await softStep('3. Click Run script link; dialog opens, click OK; grid appears', async () => {
    await page.evaluate(async () => {
      const runLink = Array.from(document.querySelectorAll('a, .d4-link-label, span, button'))
        .filter((e) => (e as HTMLElement).offsetParent !== null)
        .find((e) => e.textContent?.trim() === 'Run script') as HTMLElement;
      runLink?.click();
      await new Promise((r) => setTimeout(r, 1800));
      const ok = Array.from(document.querySelectorAll('[name="button-OK"]'))
        .find((b) => (b as HTMLElement).offsetParent !== null) as HTMLElement;
      ok?.click();
      await new Promise((r) => setTimeout(r, 12000));
    });
    const hasGrid = await page.evaluate(() => !!document.querySelector('[name="viewer-Grid"]'));
    expect(hasGrid).toBe(true);
  });

  // 4. Add viewers / change styling — intentionally skipped (canvas Toolbox)
  await test.step.skip('4. Add viewers / change styling (SKIP — canvas toolbox)', async () => {});

  await softStep('5. Save script (and layout)', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="button-Save"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 2500));
    });
  });

  await softStep('6. Close All', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1000));
    });
    const view = await page.evaluate(() => grok.shell.v?.name);
    expect(view).toBe('Home');
  });

  await softStep('7. Run test_Layout — result opens', async () => {
    const info = await page.evaluate(async () => {
      const s = await grok.dapi.scripts.filter('name = "test_Layout"').first();
      const res = await s.apply({idx: 1});
      if (res && res.rowCount) grok.shell.addTableView(res);
      await new Promise((r) => setTimeout(r, 2000));
      return {rows: res?.rowCount, viewers: Array.from(document.querySelectorAll('[name^="viewer-"]')).map((e) => e.getAttribute('name'))};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.viewers).toContain('viewer-Grid');
  });

  // 8-11. Project save / reopen / refresh — intentionally skipped (cross-cutting)
  await test.step.skip('8-11. Project save / reopen / refresh (SKIP — cross-cutting)', async () => {});

  await softStep('cleanup: remove test_Layout', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const scripts = await grok.dapi.scripts.filter('name = "test_Layout"').list();
      for (const s of scripts) await grok.dapi.scripts.delete(s);
    });
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
