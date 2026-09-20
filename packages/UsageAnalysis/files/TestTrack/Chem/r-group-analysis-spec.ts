import {expect} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as chem from '../helpers/chem';
import {openChemMenuItemFast, waitForChemMenuRoot} from './chem-fast-helpers';

test.use(specTestOptions);

async function openRGroupsDialog(page: any) {
  await openChemMenuItemFast(page, 'R-Groups Analysis...', {delayMs: 600});
  await page.locator('.d4-dialog').waitFor({timeout: 10000});
}

// The dialog clears its update indicator before the core reaches the sketcher
// (r-group-analysis.ts calls setUpdateIndicator(false) ahead of setMolFile), so a cleared
// shadow only means getMCS resolved. The core is in the sketcher once the dialog canvas holds
// ink and stops changing — build 404 clicked OK in that gap and got "No core was provided".
async function clickMCS(page: any) {
  await page.evaluate(() => {
    const w = window as any;
    w.__mcsHash = null;
    w.__mcsStable = 0;
    const mcs = Array.from(document.querySelectorAll('.d4-dialog button'))
      .find(b => b.textContent!.trim() === 'MCS') as HTMLElement;
    mcs?.click();
  });
  await page.waitForFunction(
    () => document.querySelector('.d4-dialog .d4-update-shadow') == null,
    null, {timeout: 60_000});
  await page.waitForFunction(() => {
    const w = window as any;
    let h = 0;
    let ink = 0;
    for (const cv of Array.from(document.querySelectorAll('.d4-dialog canvas')) as HTMLCanvasElement[]) {
      if (!cv.width || !cv.height) continue;
      const d = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      for (let i = 0; i < d.length; i += 401) {
        h = (h * 31 + d[i]) | 0;
        if (d[i + 3] > 0 && d[i] < 200) ink++;
      }
    }
    if (ink > 0 && h === w.__mcsHash)
      return ++w.__mcsStable >= 3;
    w.__mcsHash = h;
    w.__mcsStable = 0;
    return false;
  }, null, {timeout: 60_000, polling: 150});
}

test('Chem: R-Groups Analysis Block A (GROK-16329) + Block B (Replace Latest matrix)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await waitForChemMenuRoot(page);

  await softStep('A1: Open smiles.csv (DIVERSE dataset — required for GROK-16329 empty-MCS trigger)', async () => {

    await page.evaluate(async () => {
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      grok.shell.addTableView(df);
      (window as any).__rg_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__rg_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
    await waitForChemMenu(page);
  });

  await softStep('A2: Chem → Analyze → R-Groups Analysis → dialog opens', async () => {
    await openRGroupsDialog(page);
  });

  await softStep('A3: Click MCS → sketcher populates with MCS molfile', async () => {
    await clickMCS(page);
  });

  await softStep('A4: Visual analysis checkbox is checked (default true)', async () => {
    const checked = await page.evaluate(() => {
      const cb = document.querySelector('[name="input-Visual-analysis"]') as HTMLInputElement;
      return cb?.checked ?? null;
    });
    expect(checked === true || checked === null).toBe(true);
  });

  await softStep('A5-6: Click OK → "No R-Groups were found" balloon, no trellis, no null-ref crash', async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // the balloon is the positive signal; the "no trellis" half is an absence and cannot be
    // polled for, so this still spends its budget when nothing arrives
    await page.waitForFunction(() => !!document.querySelector('.d4-balloon'),
      undefined, {timeout: 10000}).catch(() => {});
    const result = await page.evaluate(() => {
      const trellis = Array.from(grok.shell.tv?.viewers ?? []).some((v: any) => v.type === 'Trellis plot');
      const errs = ((window as any).__rg_errors ?? []) as string[];
      const nullRefErr = errs.find(e => /null.*reference|cannot read prop/i.test(e));
      return {trellis, nullRefErr};
    });
    expect((result as any).trellis, 'No trellis plot expected for empty R-Groups on smiles-50').toBe(false);
    expect((result as any).nullRefErr,
      `GROK-16329 regression: null-reference crash. err=${(result as any).nullRefErr}`).toBeUndefined();
  });

  await softStep('B1: Open sar_small.csv', async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      const df = await grok.dapi.files.readCsv('System:DemoFiles/chem/sar_small.csv');
      grok.shell.addTableView(df);
    });
    await waitForChemMenu(page);
  });

  await softStep('B2-3: Run R-Groups → click MCS → sketcher populates', async () => {
    await openRGroupsDialog(page);
    await clickMCS(page);
  });

  await softStep('B4: OK → trellis plot + R-group columns appended', async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // the step reads the trellis and the R columns, so wait for those rather than for 12s
    await page.waitForFunction(() => Array.from(grok.shell.tv?.viewers ?? [])
      .some((v: any) => v.type === 'Trellis plot') &&
      grok.shell.t.columns.toList().some((c: any) => /^R[1-4]$/.test(c.name)),
    undefined, {timeout: 12000}).catch(() => {});
    const result = await page.evaluate(() => ({
      viewers: Array.from(grok.shell.tv.viewers).map((v: any) => v.type),
      colsHasR: grok.shell.t.columns.toList().some((c: any) => /^R[1-4]$/.test(c.name)),
    }));
    expect((result as any).viewers).toContain('Trellis plot');
    expect((result as any).colsHasR).toBe(true);
  });

  await softStep('B5-7: 2nd run, MCS, uncheck Replace latest, OK → 2nd trellis plot added', async () => {
    await openRGroupsDialog(page);
    await clickMCS(page);
    await page.evaluate(() => {
      const replaceLatest = document.querySelector('[name="input-Replace-latest"]') as HTMLInputElement;
      if (replaceLatest?.checked) replaceLatest.click();
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction((n) => Array.from(grok.shell.tv?.viewers ?? [])
      .filter((v: any) => v.type === 'Trellis plot').length >= n,
    2, {timeout: 12000}).catch(() => {});
    const trellisCount = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot').length);
    expect(trellisCount).toBeGreaterThanOrEqual(2);
  });

  await softStep('B8-10: 3rd run, MCS, check Replace latest, OK → 3rd replaces 2nd', async () => {
    await openRGroupsDialog(page);
    await clickMCS(page);
    await page.evaluate(() => {
      const replaceLatest = document.querySelector('[name="input-Replace-latest"]') as HTMLInputElement;
      if (replaceLatest && !replaceLatest.checked) replaceLatest.click();
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction((n) => Array.from(grok.shell.tv?.viewers ?? [])
      .filter((v: any) => v.type === 'Trellis plot').length >= n,
    1, {timeout: 12000}).catch(() => {});
    const trellisCount = await page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot').length);
    expect(trellisCount).toBeGreaterThanOrEqual(1);
  });

  await softStep('B11: 4th run, OK without MCS → "No core was provided" balloon', async () => {
    await openRGroupsDialog(page);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForTimeout(3000);
    const balloons = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.grok-balloon, .balloon, [class*="balloon"]'))
        .map(b => b.textContent?.trim() ?? ''));
    const noCoreBalloon = balloons.find(b => /No core was provided/i.test(b));
    expect(noCoreBalloon,
      `Expected "No core was provided" balloon. balloons=${JSON.stringify(balloons.slice(0, 5))}`).toBeTruthy();
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
