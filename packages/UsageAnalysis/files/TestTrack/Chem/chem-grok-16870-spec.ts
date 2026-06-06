// GROK-16870: hovering a Box Plot point on a Molecule-column table must not crash the RDKit cell renderer (fixed 1.22.0).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Chem: GROK-16870 RDKit cell renderer does not crash in Box Plot tooltip context', async ({page}) => {
  test.setTimeout(180_000);

  await loginToDatagrok(page);

  await softStep('Setup: close all + selenium flags', async () => {
    await page.evaluate(() => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
    });
    await page.waitForTimeout(500);
  });

  await softStep('Read smiles-50.csv + addTableView + hook console.error', async () => {
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');

      grok.shell.addTableView(df);
      (window as any).__df = df;
      (window as any).__grok16870_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__grok16870_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    });
  });

  await softStep('Wait for Chem menu registration (Molecule detector + RDKit renderer ready)', async () => {
    await waitForChemMenu(page);
    await waitForMolecule(page);
  });

  await softStep('Verify Molecule semType + add Box Plot (Chem now warm)', async () => {
    const result = await page.evaluate(() => {
      const df = (window as any).__df;
      const hasMol = df.columns.toList().some((c: any) => c.semType === 'Molecule');
      const cols = df.columns.toList().map((c: any) => ({name: c.name, semType: c.semType}));
      if (!hasMol) return {ok: false, cols};
      grok.shell.tv.addViewer('Box plot');
      return {ok: true};
    });
    if (!result.ok)
      throw new Error(`Setup failed: no Molecule column detected after 30s settle. cols=${JSON.stringify(result.cols)}`);
    await page.waitForTimeout(3000);
  });

  await softStep('Hover over Box Plot at 5 positions and let tooltip render', async () => {
    const rect = await page.evaluate(() => {
      const bp = document.querySelector('[name="viewer-Box-plot"]') as HTMLElement | null;
      if (!bp) return null;
      const r = bp.getBoundingClientRect();
      return {x: r.x, y: r.y, w: r.width, h: r.height};
    });
    expect(rect, 'Box Plot viewer node not found at [name="viewer-Box-plot"]').not.toBeNull();
    const r = rect!;
    const positions: Array<[number, number]> = [
      [r.x + r.w * 0.5, r.y + r.h * 0.5],
      [r.x + r.w * 0.3, r.y + r.h * 0.5],
      [r.x + r.w * 0.7, r.y + r.h * 0.5],
      [r.x + r.w * 0.5, r.y + r.h * 0.3],
      [r.x + r.w * 0.5, r.y + r.h * 0.7],
    ];
    for (const [hx, hy] of positions) {
      await page.mouse.move(hx, hy, {steps: 5});
      await page.waitForTimeout(1200);
    }
    await page.waitForTimeout(2000);
  });

  await softStep('Assert no rdkit-cell-renderer errors during hover sweep', async () => {
    const result = await page.evaluate(() => {
      const errs = ((window as any).__grok16870_errors as string[] | undefined) ?? [];
      const rdkitErrors = errs.filter(e =>
        /rdkit[-_]cell[-_]renderer|method not found|gS|cellRenderer\.render|NullError/i.test(e));
      return {totalErrors: errs.length, rdkitErrors};
    });
    expect(
      result.rdkitErrors.length,
      `GROK-16870 regression: RDKit cell renderer errors fired during Box Plot tooltip render. rdkitErrors=${JSON.stringify(result.rdkitErrors.slice(0, 5))}`,
    ).toBe(0);
  });

  await softStep('Assert Box Plot remains responsive after hover sweep', async () => {
    const status = await page.evaluate(() => {
      const bp = document.querySelector('[name="viewer-Box-plot"]') as HTMLElement | null;
      if (!bp) return {present: false, rect: null, canvasCount: 0};
      const r = bp.getBoundingClientRect();
      const canvasCount = bp.querySelectorAll('canvas').length;
      return {present: true, rect: {w: r.width, h: r.height}, canvasCount};
    });
    expect(status.present, 'Box Plot viewer disappeared after hover sweep — possible renderer crash').toBe(true);
    expect(
      status.rect!.w * status.rect!.h,
      `Box Plot viewer has zero size after hover — possible layout poison from renderer crash. rect=${JSON.stringify(status.rect)}`,
    ).toBeGreaterThan(0);
    expect(
      status.canvasCount,
      'Box Plot canvases unmounted after hover sweep — render loop disrupted',
    ).toBeGreaterThanOrEqual(1);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
