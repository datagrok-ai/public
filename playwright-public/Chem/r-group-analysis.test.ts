/* ---
sub_features_covered: [chem.analyze.r-groups, chem.analyze.r-groups.decomposition, chem.analyze.r-groups.top-menu, chem.sketcher]
--- */
// Paired scenario: r-group-analysis.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import * as chem from '../helpers/chem';

test.use(specTestOptions);

async function openRGroupsDialog(page: any) {
  // Retry the menu-open: on a cold session the first Chem-menu click can be
  // swallowed (semType detection still settling, transient build balloons), so
  // re-dispatch until the dialog actually mounts instead of a single 10s wait.
  for (let attempt = 0; attempt < 4; attempt++) {
    await chem.openChemMenuItem(page, 'R-Groups Analysis...', {delayMs: 800});
    try {
      await page.locator('.d4-dialog').waitFor({timeout: 8000});
      return;
    } catch (e) {
      if (attempt === 3) throw e;
      await page.keyboard.press('Escape').catch(() => {});
      await page.waitForTimeout(1000);
    }
  }
}

// Click MCS, then poll until the sketcher's update-indicator (`.d4-update-shadow`,
// added synchronously on click, removed when getMCS resolves) clears — the real
// "scaffold computed into sketcher" gate, not a blind 8s wait.
async function clickMCS(page: any) {
  await page.evaluate(() => {
    const mcs = Array.from(document.querySelectorAll('.d4-dialog button'))
      .find(b => b.textContent!.trim() === 'MCS') as HTMLElement;
    mcs?.click();
  });
  await page.waitForFunction(
    () => document.querySelector('.d4-dialog .d4-update-shadow') == null,
    null, {timeout: 60_000});
}

test('Chem: R-Groups Analysis Block A (GROK-16329) + Block B (Replace Latest matrix)', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  // ===== Block A — smiles-50.csv empty-result balloon (GROK-16329) =====

  await softStep('A1: Open smiles.csv (DIVERSE dataset — required for GROK-16329 empty-MCS trigger)', async () => {
    // smiles.csv diversity triggers MCS-cannot-decompose → empty result; smiles-50.csv is too uniform.
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
      // Capture shell.error balloon text on the grok.shell singleton — balloons
      // self-expire before a DOM/warnings scan catches them (see run.md), so hook
      // the source instead of racing the fade-out.
      (window as any).__rg_balloons = [];
      const sh: any = grok.shell;
      const origErr = sh.error.bind(sh);
      sh.error = function(msg: any) {
        (window as any).__rg_balloons.push(typeof msg === 'string' ? msg : (msg?.message ?? String(msg)));
        return origErr(msg);
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
    await page.locator('[name="input-Visual-analysis"]').waitFor({timeout: 8000});
    const checked = await page.evaluate(() =>
      (document.querySelector('[name="input-Visual-analysis"]') as HTMLInputElement).checked);
    expect(checked, 'Visual analysis defaults to true').toBe(true);
  });

  await softStep('A5-6: Click OK → "No R-Groups were found" balloon, no trellis, no null-ref crash', async () => {
    await page.evaluate(() => { (window as any).__rg_balloons = []; });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await expect.poll(async () => page.evaluate(() =>
      ((window as any).__rg_balloons as string[]).length > 0 ||
      Array.from(grok.shell.tv?.viewers ?? []).some((v: any) => v.type === 'Trellis plot')),
    {timeout: 30_000, intervals: [250, 500, 1000]}).toBe(true);
    const result = await page.evaluate(() => ({
      trellis: Array.from(grok.shell.tv?.viewers ?? []).some((v: any) => v.type === 'Trellis plot'),
      balloons: (window as any).__rg_balloons as string[],
      nullRefErr: ((window as any).__rg_errors as string[]).find(e => /null.*reference|cannot read prop/i.test(e)),
    }));
    expect(result.trellis, 'No trellis plot expected for empty R-Groups on smiles.csv').toBe(false);
    expect(result.nullRefErr,
      `GROK-16329 regression: null-reference crash. err=${result.nullRefErr}`).toBeUndefined();
    expect(result.balloons.some(b => /No R-Groups were found/i.test(b)),
      `GROK-16329: expected "No R-Groups were found" balloon. balloons=${JSON.stringify(result.balloons)}`).toBe(true);
  });

  // ===== Block B — sar_small.csv Replace Latest matrix =====

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

  await softStep('B4: OK → trellis plot + Core/R1-R4 columns appended', async () => {
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await expect.poll(async () => page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).some((v: any) => v.type === 'Trellis plot') &&
      grok.shell.t.columns.toList().some((c: any) => /^R[1-4]$/.test(c.name))),
    {timeout: 30_000, intervals: [250, 500, 1000]}).toBe(true);
    const result = await page.evaluate(() => {
      const names = grok.shell.t.columns.toList().map((c: any) => c.name);
      return {
        trellis: Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot').length,
        hasCore: names.includes('Core'),
        rCols: names.filter((n: string) => /^R[1-4]$/.test(n)).sort(),
      };
    });
    expect(result.trellis, 'first run adds exactly one trellis plot').toBe(1);
    expect(result.hasCore, 'Core scaffold column appended').toBe(true);
    expect(result.rCols, 'R1-R4 columns appended').toEqual(['R1', 'R2', 'R3', 'R4']);
  });

  await softStep('B5-7: 2nd run, MCS, uncheck Replace latest, OK → 2nd trellis + appended R-column set', async () => {
    await openRGroupsDialog(page);
    await clickMCS(page);
    const preRCount = await page.evaluate(() => {
      const replaceLatest = document.querySelector('[name="input-Replace-latest"]') as HTMLInputElement;
      if (replaceLatest?.checked) replaceLatest.click();
      return grok.shell.t.columns.toList().filter((c: any) => /^R\d/.test(c.name)).length;
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await expect.poll(async () => page.evaluate(() =>
      Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot').length),
    {timeout: 30_000, intervals: [250, 500, 1000]}).toBe(2);
    const state = await page.evaluate(() => {
      const names = grok.shell.t.columns.toList().map((c: any) => c.name);
      return {
        rCount: names.filter((n: string) => /^R\d/.test(n)).length,
        hasR14: ['R1', 'R2', 'R3', 'R4'].every((n) => names.includes(n)),
      };
    });
    expect(state.hasR14, 'first run R1-R4 preserved after append').toBe(true);
    expect(state.rCount, `second R-group set appended (was ${preRCount})`).toBeGreaterThan(preRCount);
  });

  await softStep('B8-10: 3rd run, MCS, check Replace latest, OK → 3rd replaces 2nd (2 trellis, first preserved)', async () => {
    await openRGroupsDialog(page);
    await clickMCS(page);
    // Tag existing trellis roots so we can gate on the 3rd run producing a fresh
    // (untagged) trellis while the total stays 2 — polling a bare count of 2 would
    // pass instantly on the pre-existing append state without the replace running.
    const preRCount = await page.evaluate(() => {
      const replaceLatest = document.querySelector('[name="input-Replace-latest"]') as HTMLInputElement;
      if (replaceLatest && !replaceLatest.checked) replaceLatest.click();
      Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot')
        .forEach((v: any) => v.root.setAttribute('data-rg-pre', '1'));
      return grok.shell.t.columns.toList().filter((c: any) => /^R\d/.test(c.name)).length;
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await expect.poll(async () => page.evaluate(() => {
      const trellis = Array.from(grok.shell.tv.viewers).filter((v: any) => v.type === 'Trellis plot');
      return {total: trellis.length, fresh: trellis.filter((v: any) => !v.root.hasAttribute('data-rg-pre')).length};
    }), {timeout: 30_000, intervals: [250, 500, 1000]}).toEqual({total: 2, fresh: 1});
    const state = await page.evaluate(() => {
      const names = grok.shell.t.columns.toList().map((c: any) => c.name);
      return {
        rCount: names.filter((n: string) => /^R\d/.test(n)).length,
        hasR14: ['R1', 'R2', 'R3', 'R4'].every((n) => names.includes(n)),
      };
    });
    expect(state.hasR14, 'first run R1-R4 remain intact after replace').toBe(true);
    expect(state.rCount, 'replace swaps latest set (no net R-column growth)').toBe(preRCount);
  });

  await softStep('B11: 4th run, OK without MCS → "No core was provided" balloon', async () => {
    await openRGroupsDialog(page);
    await page.evaluate(() => { (window as any).__rg_balloons = []; });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await expect.poll(async () => page.evaluate(() =>
      ((window as any).__rg_balloons as string[]).slice()),
    {timeout: 15_000, intervals: [200, 400, 800]}).toEqual(
      expect.arrayContaining([expect.stringMatching(/No core was provided/i)]));
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
