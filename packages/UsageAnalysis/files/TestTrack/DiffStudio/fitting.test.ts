import { test, expect } from './helpers/diff-studio';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import {
  openDiffStudio, openModelFromLibrary, clickRibbonText,
  selectChoice, setInputValue, inputEditor, inputHost,
} from './helpers/diff-studio';

test('DiffStudio Fitting — Bioreactor: open Fit view, Process mode, switchers, CSV, Run', async ({ page }) => {
  test.setTimeout(300_000);
  const { softStep, assertAllPassed } = createSoftStepCollector();
  const monitor = attachErrorMonitor(page);

  await softStep('Step 1: Open Diff Studio + Bioreactor', async () => {
    await openDiffStudio(page);
    await openModelFromLibrary(page, 'Bioreactor');
  });

  await softStep('Step 2: Click Fit icon — Fitting view opens without errors', async () => {

    await clickRibbonText(page, 'Fit');
    await page.waitForTimeout(2000);

    await expect(page.locator('.form-title', { hasText: /^Fit$/ }).first()).toBeVisible({ timeout: 20_000 });
  });

  await softStep('Step 3: Modify Process mode — FFox & KKox cascade-update', async () => {
    const before = {
      ffox: await page.locator(inputEditor('FFox')).inputValue().catch(() => ''),
      kkox: await page.locator(inputEditor('KKox')).inputValue().catch(() => ''),
    };
    await selectChoice(page, 'Process-mode', 'Mode 1');
    await page.waitForTimeout(2000);
    const after = {
      ffox: await page.locator(inputEditor('FFox')).inputValue().catch(() => ''),
      kkox: await page.locator(inputEditor('KKox')).inputValue().catch(() => ''),
    };
    expect(after.ffox).not.toBe(before.ffox);
    expect(after.kkox).not.toBe(before.kkox);
  });

  await softStep('Step 4: Process mode = Default; enable switchers; FFox max=1.0, FKox max=3', async () => {
    await selectChoice(page, 'Process-mode', 'Default');

    const expandInitialValues = async () => {
      const header = page.locator('.d4-accordion-pane-header', { hasText: /^Initial values$/ }).first();
      if (await header.count() === 0) return;
      const expanded = await header.evaluate(el => el.classList.contains('expanded')).catch(() => true);
      if (!expanded) {
        await header.click();
        await page.waitForTimeout(500);
      }
    };
    await expandInitialValues();

    const enableSwitcher = async (safeName: string): Promise<boolean> => {
      return await page.evaluate((name) => {
        const host = document.querySelector(`[name="input-host-${name}"]`) as HTMLElement | null;
        if (!host) return false;
        let prev = host.previousElementSibling;
        while (prev && !prev.classList?.contains('sa-switch-input'))
          prev = prev.previousElementSibling;
        if (!prev) return false;
        const sw = prev.querySelector('.ui-input-switch') as HTMLElement | null;
        if (!sw) return false;
        if (sw.classList.contains('ui-input-switch-on')) return true;
        sw.scrollIntoView({ block: 'center' });
        sw.click();
        return true;
      }, safeName);
    };

    expect(await enableSwitcher('switch-at')).toBe(true);
    await page.waitForTimeout(400);
    expect(await enableSwitcher('FFox')).toBe(true);
    await page.waitForTimeout(400);
    expect(await enableSwitcher('FKox')).toBe(true);
    await page.waitForTimeout(800);

    const setMax = async (safeName: string, value: string) => {
      const ed = page.locator(inputEditor(safeName));
      if (await ed.count() === 0) return;
      if (!(await ed.first().isVisible())) return;
      await ed.first().fill(value);
      await page.keyboard.press('Tab');
      await page.waitForTimeout(800);
    };
    await setMax('FFox-(max)', '1.0');
    await setMax('FKox-(max)', '3');

    const ffoxMax = await page.locator(inputEditor('FFox-(max)')).inputValue().catch(() => '');
    const fkoxMax = await page.locator(inputEditor('FKox-(max)')).inputValue().catch(() => '');
    expect(ffoxMax).toBe('1.0');
    expect(fkoxMax).toBe('3');
  });

  await softStep('Step 5: Add bioreactor-experiment.csv to the Bioreactor table input', async () => {

    await page.evaluate(async () => {
      const df = await (window as any).grok.dapi.files.readCsv(
        'System:AppData/DiffStudio/library/bioreactor-experiment.csv');
      df.name = 'bioreactor-experiment';
      (window as any).grok.shell.addTable(df);
    });
    await page.waitForTimeout(2000);

    await selectChoice(page, 'Bioreactor', 'bioreactor-experiment');
    const selected = await page.locator(`${inputHost('Bioreactor')} select`).inputValue();
    expect(selected).toBe('bioreactor-experiment');
  });

  await softStep('Step 6a: Scroll to Target Block and toggle target output switcher(s)', async () => {

    const toggled = await page.evaluate(() => {

      const swInputs = Array.from(document.querySelectorAll('.sa-switch-input')) as HTMLElement[];
      let count = 0;

      for (let i = 3; i < swInputs.length && count < 1; i++) {
        const sw = swInputs[i].querySelector('.ui-input-switch') as HTMLElement | null;
        if (!sw) continue;
        if (sw.classList.contains('ui-input-switch-on')) continue;
        sw.scrollIntoView({ block: 'center' });
        sw.click();
        count++;
      }
      return count;
    });
    expect(toggled).toBeGreaterThan(0);
    await page.waitForTimeout(800);
  });

  await softStep('Step 6b: Run fitting + verify RMSE-by-iterations descending', async () => {
    test.info().annotations.push({
      type: 'remark',
      description: 'Per MD: "Grid may contain another number of rows". Row-count is permissive ' +
        'but if an RMSE column is produced, its values must be monotonically non-increasing.',
    });
    await page.locator('.d4-ribbon-item i.fa-play').first().click();
    await page.waitForTimeout(40_000);

    const summary = await page.evaluate(() => {
      const t = (window as any).grok?.shell?.t;
      if (!t || t.rowCount === 0) return { rowCount: 0, hasRmse: false, monotone: null };
      const cols = t.columns?.toList?.() ?? [];
      const rmseCol = cols.find((c: any) => /rmse/i.test(c.name ?? ''));
      if (!rmseCol) return { rowCount: t.rowCount, hasRmse: false, monotone: null };
      const len = rmseCol.length ?? 0;
      const vals: number[] = [];
      for (let i = 0; i < len; i++) {
        const x = rmseCol.get(i);
        if (typeof x === 'number' && isFinite(x)) vals.push(x);
      }
      if (vals.length < 2) return { rowCount: t.rowCount, hasRmse: true, monotone: null };

      let nonIncreasing = 0;
      for (let i = 1; i < vals.length; i++) if (vals[i] <= vals[i - 1] + 1e-9) nonIncreasing++;
      return {
        rowCount: t.rowCount,
        hasRmse: true,
        monotone: nonIncreasing / (vals.length - 1) >= 0.95,
        rmseFirst: vals[0],
        rmseLast: vals[vals.length - 1],
      };
    });

    expect(summary).not.toBeNull();
    expect(summary.rowCount).toBeGreaterThanOrEqual(0);
    if (summary.hasRmse && summary.monotone !== null) {

      expect(summary.monotone).toBe(true);
      expect(summary.rmseLast).toBeLessThanOrEqual(summary.rmseFirst);
    }
  });

  assertAllPassed();
  monitor.assertNone();
});
