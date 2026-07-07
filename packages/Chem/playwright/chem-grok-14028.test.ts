/* ---
sub_features_covered: [chem.actions.copy-smiles, chem.search.substructure, chem.search.substructure.filter, chem.sketcher]
--- */
// GROK-14028: substructure-filter Clear must clear all 3 layers — L1 BitSet, L2 sketcher UI/summary, L3 leaked tags.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForMolecule} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

test('Chem: GROK-14028 Filter Panel Clear 3-layer cleanup invariant', async ({page}) => {
  test.setTimeout(120_000);

  await loginToDatagrok(page);

  await softStep('Setup: close all + selenium flags', async () => {
    await page.evaluate(() => {
      document.body.classList.add('selenium');
      try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { grok.shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
    });
  });

  await softStep('Read spgi-100.csv + addTableView + hook console.error', async () => {
    await page.evaluate(async (path) => {
      const df = await grok.dapi.files.readCsv(path);

      grok.shell.addTableView(df);
      (window as any).__df = df;
      (window as any).__grok14028_errors = [];
      const orig = console.error;
      console.error = function(...args: any[]) {
        (window as any).__grok14028_errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    }, datasetPath);
  });

  await softStep('Wait for Molecule semType (Chem autostart cascade)', async () => {
    await waitForMolecule(page);
  });

  await softStep('Verify Molecule semType detected + Grid renders', async () => {
    const result = await page.evaluate(() => {
      const df = (window as any).__df;
      const cols = df.columns.toList().map((c: any) => ({name: c.name, semType: c.semType}));
      const hasMol = cols.some((c: any) => c.semType === 'Molecule');
      return {hasMol, cols};
    });
    if (!result.hasMol)
      throw new Error(`Setup failed: no Molecule column on spgi-100.csv. cols=${JSON.stringify(result.cols)}`);
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 10000});
  });

  await softStep('Open Filter Panel and wait for Structure filter sketch-link', async () => {
    await page.evaluate(() => grok.shell.tv.getFiltersGroup());
    await page.locator('[name="viewer-Filters"] .d4-filter').first().waitFor({timeout: 30000});
    const probeResult = await page.evaluate(async () => {
      for (let i = 0; i < 90; i++) {
        const filters = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
        const headers = Array.from(filters).map(f => f.querySelector('.d4-filter-header')?.textContent?.trim());
        for (const f of filters) {
          const header = f.querySelector('.d4-filter-header');
          if (!header || !/Structure|struct/i.test(header.textContent || '')) continue;
          const sketchLink = f.querySelector('.sketch-link');
          if (sketchLink) return {ready: true, atMs: (i+1)*500};
        }
        if (i === 89) {
          return {
            ready: false,
            totalFilters: filters.length,
            headers: headers.slice(0, 20),
            sketchLinkInDom: !!document.querySelector('.sketch-link'),
            chemFilterInDom: !!document.querySelector('.chem-filter'),
            grokSketcherInDom: !!document.querySelector('.grok-sketcher'),
          };
        }
        await new Promise(r => setTimeout(r, 500));
      }
      return {ready: false};
    });
    expect(
      (probeResult as any).ready,
      `Structure filter sketch-link not rendered within 45s. Diagnostic: ${JSON.stringify(probeResult)}`,
    ).toBe(true);
  });

  await softStep('Apply benzene substructure filter via Structure filter sketch-link', async () => {
    const sketcherOpened = await page.evaluate(async () => {
      const filters = document.querySelectorAll('[name="viewer-Filters"] .d4-filter');
      let clicked = false;
      for (const f of filters) {
        const header = f.querySelector('.d4-filter-header');
        if (header && /Structure|struct/i.test(header.textContent || '')) {
          const sketchLink = f.querySelector('.sketch-link');
          if (sketchLink) { (sketchLink as HTMLElement).click(); clicked = true; }
          break;
        }
      }
      // Sketcher dialog is position:fixed (offsetParent null but visible) — detect via getBoundingClientRect.
      for (let i = 0; i < 25; i++) {
        await new Promise(r => setTimeout(r, 400));
        const dlg = document.querySelector('.d4-dialog') as HTMLElement | null;
        if (dlg) {
          const rect = dlg.getBoundingClientRect();
          if (rect.width > 100 && rect.height > 100) return clicked;
        }
      }
      return false;
    });
    expect(sketcherOpened, 'Sketcher dialog did not open via Structure filter sketch-link').toBe(true);

    // Fill SMILES input via DOM event sequence — native fill() may not reach Dart-side state.
    // Then poll for the OK button becoming clickable (replaces a fixed parse-wait) and click it.
    await page.evaluate(async () => {
      const smilesInput = Array.from(document.querySelectorAll('.d4-dialog input'))
        .find((i: any) => /smiles/i.test(i.placeholder || '')) as HTMLInputElement | undefined;
      if (!smilesInput) throw new Error('SMILES input not found in sketcher dialog');
      smilesInput.focus();
      smilesInput.value = 'c1ccccc1';
      smilesInput.dispatchEvent(new Event('input', {bubbles: true}));
      smilesInput.dispatchEvent(new Event('change', {bubbles: true}));
      smilesInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true}));
      for (let i = 0; i < 40; i++) {
        const okBtn = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement | null;
        if (okBtn && !okBtn.classList.contains('disabled')) { okBtn.click(); return; }
        await new Promise(r => setTimeout(r, 250));
      }
      throw new Error('OK button never became clickable in sketcher dialog');
    });
    const total = await page.evaluate(() => grok.shell.t.rowCount);
    await expect.poll(
      () => page.evaluate(() => grok.shell.t.filter.trueCount),
      {timeout: 30000, message: 'benzene substructure filter never reduced the visible row count'},
    ).toBeLessThan(total);
    const filterApplied = await page.evaluate(() => ({
      filtered: grok.shell.t.filter.trueCount,
      total: grok.shell.t.rowCount,
      summary: document.querySelector('[name="viewer-Filters"] .d4-filter .d4-filter-summary')?.textContent?.trim() ?? '',
      clearBtnPresent: !!document.querySelector('[name="viewer-Filters"] .d4-filter .chem-clear-sketcher-button'),
    }));
    expect(filterApplied.filtered).toBeLessThan(filterApplied.total);
    expect(filterApplied.filtered).toBeGreaterThan(0);
    // Summary must show an actual substructure token (aromatic-ring SMILES or SMARTS atom),
    // not merely any uppercase letter — [A-Z] catch-all removed.
    expect(filterApplied.summary, 'Filter summary should contain the benzene substructure pattern after apply').toMatch(/\[#\d+\]|c1[a-z0-9]/i);
    expect(filterApplied.clearBtnPresent, 'Clear button should appear after substructure filter applied').toBe(true);
  });

  await softStep('Click Clear button (Reset substructure filter)', async () => {
    // UI-level reset via the inline clear-sketcher button — a real DOM click, not a JS-API
    // BitSet substitution (df.filter.setAll(true)), which the .md forbids for the trigger.
    await page.evaluate(() => {
      const clearBtn = document.querySelector('[name="viewer-Filters"] .d4-filter .chem-clear-sketcher-button') as HTMLElement | null;
      if (!clearBtn) throw new Error('Clear button not found — filter may not have been applied successfully');
      clearBtn.click();
    });
    // Gate the 3-layer read on L1 (BitSet) actually clearing, instead of a fixed sleep.
    await expect.poll(
      () => page.evaluate(() => grok.shell.t.filter.trueCount === grok.shell.t.rowCount),
      {timeout: 15000, message: 'BitSet did not return to all-true after Reset'},
    ).toBe(true);
  });

  await softStep('Assert 3-layer cleanup invariant (L1 BitSet + L2 Sketcher UI + L3 No leaked tags)', async () => {
    const result = await page.evaluate(() => {
      const df = grok.shell.t;
      // L1: BitSet all-true (no filter applied)
      const L1_bitsetCleared = df.filter.trueCount === df.rowCount;
      // L2: sketcher UI cleared — this is the exact GROK-14028 symptom (stale input line).
      //   Primary signal: the atlas-documented chem-substructure-filter column tag is empty/null.
      //   Corroboration: .chem-clear-sketcher-button gone AND .d4-filter-summary has no pattern.
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const substructTag = molCol ? molCol.getTag('chem-substructure-filter') : null;
      const L2_tagCleared = substructTag == null || String(substructTag).length === 0;
      const clearBtnGone = !document.querySelector('[name="viewer-Filters"] .d4-filter .chem-clear-sketcher-button');
      const summary = document.querySelector('[name="viewer-Filters"] .d4-filter .d4-filter-summary')?.textContent?.trim() ?? '';
      const summaryEmptyOfSmiles = !/\[#\d+\]|c1[a-zA-Z0-9]/.test(summary);
      const L2_sketcherCleared = L2_tagCleared && clearBtnGone && summaryEmptyOfSmiles;
      // L3: no leaked virtual / tag columns
      const cols = df.columns.names();
      const leakedTags = cols.filter((n: string) => /^~.*(substructure|highlight)/i.test(n));
      const L3_noLeakedTags = leakedTags.length === 0;
      const consoleErrors = ((window as any).__grok14028_errors as string[] | undefined) ?? [];
      return {
        L1_bitsetCleared,
        L2_sketcherCleared,
        L3_noLeakedTags,
        L2_tagCleared,
        substructTag,
        clearBtnGone,
        summary,
        summaryEmptyOfSmiles,
        leakedTags,
        consoleErrors,
        filteredCount: df.filter.trueCount,
        rowCount: df.rowCount,
      };
    });
    const allLayersCleared = result.L1_bitsetCleared && result.L2_sketcherCleared && result.L3_noLeakedTags;
    expect(
      allLayersCleared,
      `GROK-14028 regression: 3-layer cleanup incomplete. ` +
      `L1(BitSet=${result.filteredCount}/${result.rowCount})=${result.L1_bitsetCleared} ` +
      `L2(tag='${result.substructTag}' tagCleared=${result.L2_tagCleared} clearBtnGone=${result.clearBtnGone} summary='${result.summary}' summaryEmptyOfSmiles=${result.summaryEmptyOfSmiles})=${result.L2_sketcherCleared} ` +
      `L3(leaked=${JSON.stringify(result.leakedTags)})=${result.L3_noLeakedTags}`,
    ).toBe(true);
    expect(
      result.consoleErrors.filter((e: string) => /TypeError|null|undefined|cannot read/i.test(e)).length,
      `Unexpected console errors during Clear flow: ${JSON.stringify(result.consoleErrors)}`,
    ).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
