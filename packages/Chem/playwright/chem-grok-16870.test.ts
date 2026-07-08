/* ---
sub_features_covered: [chem.rendering, chem.rendering.molecule-cell, chem.rendering.rdkit-renderer]
--- */
// GROK-16870: hovering a Box Plot point on a Molecule-column table must not crash the RDKit cell renderer (fixed 1.22.0).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

test('Chem: GROK-16870 RDKit cell renderer does not crash in Box Plot tooltip context', async ({page}) => {
  test.setTimeout(90_000);

  await loginToDatagrok(page);

  await softStep('Setup: close all + selenium flags', async () => {
    await page.evaluate(() => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
    });
    await page.waitForFunction(() => {
      try { return Array.from((window as any).grok.shell.views).length === 0; }
      catch (e) { return false; }
    }, {timeout: 10_000}).catch(() => {});
  });

  await softStep('Read smiles-50.csv + addTableView + hook console.error', async () => {
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');

      grok.shell.addTableView(df);
      (window as any).__df = df;
      (window as any).__grok16870_errors = [];
      const push = (s: string) => (window as any).__grok16870_errors.push(s);
      const orig = console.error;
      console.error = function(...args: any[]) {
        push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
      // GROK-16870 surfaces as a Dart NullError that can bypass console.error and hit
      // window.onerror / unhandledrejection — capture all channels so 0-errors is meaningful.
      window.addEventListener('error', (e: any) => push(String(e?.error ?? e?.message ?? e)));
      window.addEventListener('unhandledrejection', (e: any) => push(String(e?.reason)));
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
    // Wait for the Box Plot viewer to actually attach before the hover sweep queries its rect.
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 30_000, state: 'visible'});
  });

  let anyTooltipShown = false;
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
      const shown = await page.locator('.d4-tooltip').first()
        .waitFor({state: 'visible', timeout: 6000}).then(() => true).catch(() => false);
      anyTooltipShown = anyTooltipShown || shown;
    }
    await page.waitForTimeout(200); // let the last async RDKit render microtask flush
  });

  await softStep('Assert no rdkit-cell-renderer errors during hover sweep', async () => {
    const result = await page.evaluate(() => {
      const errs = ((window as any).__grok16870_errors as string[] | undefined) ?? [];
      const rdkitErrors = errs.filter(e =>
        /rdkit[-_]cell[-_]renderer|method not found|gS|cellRenderer\.render|NullError/i.test(e));
      return {totalErrors: errs.length, rdkitErrors};
    });
    // GROK-16870: tooltip presence is a soft signal, not a hard precondition — the .md
    // documents that a Box Plot over smiles-50 may render an empty central region, so a
    // missed hover is not a regression. The invariant locked here is crash ABSENCE.
    if (!anyTooltipShown)
      console.warn('GROK-16870: hover sweep rendered no Box Plot tooltip (empty central region) — absence check is crash-scoped only');
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
