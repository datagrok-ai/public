/* ---
sub_features_covered: [chem.notation, chem.notation.action, chem.notation.convert-mol]
--- */
// GROK-17964: Convert Notation column-action must register exactly once across cancel/commit/repeat invocations.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

// Focus the original molecule column on the Context Panel, poll until the Convert Notation
// action link renders, then return how many are attached — replaces the duplicated
// set-current-object + blind-sleep + count blocks.
async function countConvertNotationOnMolCol(page: Page): Promise<number> {
  await page.evaluate(() => { grok.shell.o = grok.shell.t.col((window as any).__grok17964_origMolCol); });
  await expect.poll(async () => page.evaluate(() =>
    Array.from(document.querySelectorAll('label.d4-link-action'))
      .some(l => (l.textContent ?? '').trim().startsWith('Convert Notation')),
  ), {timeout: 15_000, intervals: [250, 500, 1000]}).toBe(true);
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('label.d4-link-action'))
      .filter(l => (l.textContent ?? '').trim().startsWith('Convert Notation')).length);
}

test('Chem: GROK-17964 Convert Notation column-action registration is exactly-once', async ({page}) => {
  test.setTimeout(120_000);

  await loginToDatagrok(page);

  await softStep('Setup: close all + selenium flags', async () => {
    await page.evaluate(() => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
    });
    await expect.poll(async () => page.evaluate(() => Array.from(grok.shell.tableViews).length),
      {timeout: 5000, intervals: [100, 250, 500]}).toBe(0);
  });

  await softStep('Read smiles-50.csv + addTableView', async () => {
    await page.evaluate(async () => {
      const df = await grok.dapi.files.readCsv('System:AppData/Chem/tests/smiles-50.csv');

      grok.shell.addTableView(df);
      (window as any).__df = df;
    });
  });

  await softStep('Wait for Chem menu registration (Molecule semType + action @autostart ready)', async () => {
    await waitForChemMenu(page);
  });

  await softStep('Find molecule column + focus column on Context Panel + expand panes', async () => {
    const result = await page.evaluate(async () => {
      for (let i = 0; i < 30; i++) {
        const df = grok.shell.t;
        const molColName = df?.columns.toList().find((c: any) => c.semType === 'Molecule')?.name;
        if (molColName) {
          (window as any).__df = df;
          grok.shell.o = df.col(molColName);
          (window as any).__grok17964_origMolCol = molColName;
          return {ok: true, molColName};
        }
        await new Promise(r => setTimeout(r, 1000));
      }
      const df = grok.shell.t;
      const allCols = df?.columns.toList().map((c: any) => ({name: c.name, semType: c.semType})) ?? [];
      return {ok: false, molColName: null, allCols};
    });
    if (!result.ok)
      throw new Error(`Setup failed: no Molecule column detected on smiles-50.csv after 30s poll. cols=${JSON.stringify(result.allCols)}`);
    await page.locator('.d4-accordion-pane').first().waitFor({state: 'attached', timeout: 10_000});
    // Expand all accordion panes — chem action labels render only when the Actions pane is expanded.
    await page.evaluate(async () => {
      const panes = Array.from(document.querySelectorAll('.d4-accordion-pane'));
      for (const p of panes) {
        const h = p.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (h && !h.classList.contains('expanded')) {
          h.click();
          await new Promise(r => setTimeout(r, 100));
        }
      }
    });
    // Poll for the Convert Notation action link to render on the expanded Actions pane instead of a blind sleep.
    await expect.poll(async () => page.evaluate(() =>
      Array.from(document.querySelectorAll('label.d4-link-action'))
        .some(l => (l.textContent ?? '').trim().startsWith('Convert Notation')),
    ), {timeout: 15_000, intervals: [250, 500, 1000]}).toBe(true);
  });

  await softStep('Baseline: assert exactly 1 Convert Notation entry on the column Actions pane', async () => {
    const baseline = await page.evaluate(() => {
      const entries = Array.from(document.querySelectorAll('label.d4-link-action'))
        .filter(l => (l.textContent ?? '').trim().startsWith('Convert Notation'));
      return {count: entries.length, sample: entries.slice(0, 3).map(e => (e.textContent ?? '').trim())};
    });
    expect(
      baseline.count,
      `GROK-17964 baseline regression: initial Convert Notation registration count expected 1, got ${baseline.count}. samples=${JSON.stringify(baseline.sample)}`,
    ).toBe(1);
  });

  await softStep('Cancellation path: open Convert Notation dialog, assert controls, CANCEL, recount', async () => {
    await page.evaluate(() => {
      const link = Array.from(document.querySelectorAll('label.d4-link-action'))
        .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation')) as HTMLElement;
      if (!link) throw new Error('Convert Notation link not found pre-cancel');
      link.click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 8000});
    const controls = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog')!;
      const q = (n: string) => dlg.querySelector(`[name="${n}"]`) as HTMLInputElement | null;
      const overwrite = q('input-Overwrite'); const join = q('input-Join'); const kekulize = q('input-Kekulize');
      return {
        hasTarget: !!q('input-Target-Notation'), hasOverwrite: !!overwrite, hasJoin: !!join, hasKekulize: !!kekulize,
        overwriteChecked: !!overwrite?.checked, joinChecked: !!join?.checked, kekulizeChecked: !!kekulize?.checked,
      };
    });
    expect(controls.hasTarget && controls.hasOverwrite && controls.hasJoin && controls.hasKekulize,
      `GROK-17964: Convert Notation dialog must expose Target-Notation/Overwrite/Join/Kekulize controls, got ${JSON.stringify(controls)}`).toBe(true);
    expect(controls.joinChecked, 'Convert Notation: Join must default checked').toBe(true);
    expect(controls.overwriteChecked, 'Convert Notation: Overwrite must default unchecked').toBe(false);
    expect(controls.kekulizeChecked, 'Convert Notation: Kekulize must default unchecked').toBe(false);
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
    await page.locator('.d4-dialog').waitFor({state: 'detached', timeout: 8000});
    const afterCancel = await countConvertNotationOnMolCol(page);
    expect(
      afterCancel,
      `GROK-17964 regression: registration count after CANCEL expected 1, got ${afterCancel}.`,
    ).toBe(1);
  });

  await softStep('Successful completion path: Convert Notation → molblock, OK, verify new molblock column', async () => {
    const baseColCount = await page.evaluate(() => grok.shell.t.columns.length);
    await page.evaluate(() => {
      const link = Array.from(document.querySelectorAll('label.d4-link-action'))
        .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation')) as HTMLElement;
      if (!link) throw new Error('Convert Notation link not found pre-commit');
      link.click();
    });
    await page.locator('.d4-dialog').waitFor({timeout: 8000});
    await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      const targetSelect = dlg?.querySelector('[name="input-Target-Notation"]') as HTMLSelectElement | null;
      if (!dlg || !targetSelect) throw new Error('Convert Notation dialog / Target-Notation select not found');
      targetSelect.value = 'molblock';
      targetSelect.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // Completion is verified, not slept: the dialog detaches and (defaults overwrite=false/join=true)
    // a new converted Molecule column is added to the frame.
    await page.locator('.d4-dialog').waitFor({state: 'detached', timeout: 30_000});
    await expect.poll(async () => page.evaluate(() => grok.shell.t.columns.length),
      {timeout: 30_000, intervals: [500, 1000, 2000]}).toBe(baseColCount + 1);
    const added = await page.evaluate(() => {
      const cols = grok.shell.t.columns.toList();
      const last = cols[cols.length - 1];
      return {name: last.name, semType: last.semType, units: last.meta?.units};
    });
    expect(
      added.semType,
      `GROK-17964: Convert Notation must add a new Molecule column (before=${baseColCount}, added=${JSON.stringify(added)}).`,
    ).toBe('Molecule');
  });

  await softStep('Exactly-once on original column post-commit', async () => {
    const onOriginal = await countConvertNotationOnMolCol(page);
    expect(
      onOriginal,
      `GROK-17964 regression: registration count on ORIGINAL column post-commit expected 1, got ${onOriginal}.`,
    ).toBe(1);
  });

  await softStep('Multi-invocation hardening: open + CANCEL twice on original column, recount', async () => {
    await page.evaluate(() => { grok.shell.o = grok.shell.t.col((window as any).__grok17964_origMolCol); });
    await expect.poll(async () => page.evaluate(() =>
      Array.from(document.querySelectorAll('label.d4-link-action'))
        .some(l => (l.textContent ?? '').trim().startsWith('Convert Notation')),
    ), {timeout: 15_000, intervals: [250, 500, 1000]}).toBe(true);
    for (let i = 0; i < 2; i++) {
      await page.evaluate(() => {
        const link = Array.from(document.querySelectorAll('label.d4-link-action'))
          .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation')) as HTMLElement;
        if (!link) throw new Error('Convert Notation link not found during multi-cancel');
        link.click();
      });
      await page.locator('.d4-dialog').waitFor({timeout: 8000});
      await page.locator('.d4-dialog [name="button-CANCEL"]').click();
      await page.locator('.d4-dialog').waitFor({state: 'detached', timeout: 8000});
    }
    const afterMulti = await countConvertNotationOnMolCol(page);
    expect(
      afterMulti,
      `GROK-17964 regression: registration count after multi-cancel expected 1, got ${afterMulti}.`,
    ).toBe(1);
  });

  await softStep('Global final assertion: only one panel-attached Convert Notation entry visible', async () => {
    const globalCount = await page.evaluate(() => {
      const entries = Array.from(document.querySelectorAll('label.d4-link-action'))
        .filter(l => (l.textContent ?? '').trim().startsWith('Convert Notation'));
      return entries.length;
    });
    expect(
      globalCount,
      `GROK-17964 regression: final global panel-attached Convert Notation count expected 1, got ${globalCount}.`,
    ).toBe(1);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
