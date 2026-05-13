/* ---
sub_features_covered: [chem.notation, chem.notation.action, chem.notation.convert-mol]
--- */
// Frontmatter extraction (Section 1 of automator-prompt):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: [chem.notation, .action, .convert-mol]
//   ui_coverage_responsibility: [] (delegated_to: null) — bug-focused slice, no specialty flow ownership
//   related_bugs: [GROK-17964]
//
// Bug invariant (regression-lock per references/bug-library/chem.yaml :: GROK-17964):
//   The `Convert Notation...` action in the column-level Context Panel is
//   registered EXACTLY ONCE per applicable molecule column. The registration
//   count remains 1 across the action lifecycle — cancellation, successful
//   completion (which creates a new molecule column), repeated invocations.
//   Duplicate registration after a handler invocation is the GROK-17964 surface.
// Parallel-coverage with info-panels.md (Context Panel walk, coverage_type=regression)
// and sketcher.md (notation round-trip via clipboard, coverage_type=regression).
//
// Paired scenario: chem-grok-17964.md
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, waitForChemMenu} from '../spec-login';

// storageState consumes auth.json captured via `npx playwright codegen --save-storage=auth.json`
// (cycle bootstrap 2026-05-11; DATAGROK_LOGIN/PASSWORD env not set in this session).
// Future cycles set credentials in env and remove the storageState override.
test.use({...specTestOptions, storageState: 'auth.json'});

test('Chem: GROK-17964 Convert Notation column-action registration is exactly-once', async ({page}) => {
  test.setTimeout(180_000);

  // GROK-17964 is status: fixed per bug-library. Spec exercises the post-fix
  // invariant as a regression-lock. Chem package autostart fires AS A
  // SIDE-EFFECT of opening a dataset with a Molecule column — wait 10s
  // AFTER addTableView (not before).
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
    // Use grok.shell.t (active table) instead of cached __df — semType detection
    // may complete on the shell-tracked df ref even when the __df reference was
    // captured earlier without the detector having fired yet.
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
    await page.waitForTimeout(2000);
    // Expand all accordion panes — chem action labels render only when the
    // containing Actions pane is expanded. In fresh Playwright sessions panes
    // default to collapsed; MCP warm sessions have them pre-expanded.
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
    await page.waitForTimeout(1500);
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

  await softStep('Cancellation path: open Convert Notation dialog, CANCEL, recount', async () => {
    await page.evaluate(async () => {
      const link = Array.from(document.querySelectorAll('label.d4-link-action'))
        .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation')) as HTMLElement;
      if (!link) throw new Error('Convert Notation link not found pre-cancel');
      link.click();
      await new Promise(r => setTimeout(r, 1500));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 8000});
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
    await page.waitForTimeout(1500);
    const afterCancel = await page.evaluate(() => {
      const molColName = (window as any).__grok17964_origMolCol;
      grok.shell.o = grok.shell.t.col(molColName);
      return new Promise(resolve => {
        setTimeout(() => {
          const entries = Array.from(document.querySelectorAll('label.d4-link-action'))
            .filter(l => (l.textContent ?? '').trim().startsWith('Convert Notation'));
          resolve({count: entries.length});
        }, 1800);
      });
    }) as {count: number};
    expect(
      afterCancel.count,
      `GROK-17964 regression: registration count after CANCEL expected 1, got ${afterCancel.count}.`,
    ).toBe(1);
  });

  await softStep('Successful completion path: Convert Notation → molblock, OK, wait for completion', async () => {
    await page.evaluate(async () => {
      // Defensive: ensure no stale dialog open before reopening
      const stale = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
      if (stale && (stale.closest('.d4-dialog') as HTMLElement | null)?.offsetParent !== null) stale.click();
      await new Promise(r => setTimeout(r, 500));
      const link = Array.from(document.querySelectorAll('label.d4-link-action'))
        .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation')) as HTMLElement;
      if (!link) throw new Error('Convert Notation link not found pre-commit');
      link.click();
      await new Promise(r => setTimeout(r, 1500));
      const dlg = document.querySelector('.d4-dialog');
      const targetSelect = dlg?.querySelector('[name="input-Target-Notation"]') as HTMLSelectElement;
      if (targetSelect) {
        // Valid target options (MCP-validated 2026-05-11): smiles, cxsmiles, smarts, cxsmarts, molblock, v3Kmolblock
        targetSelect.value = 'molblock';
        targetSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // Wait for conversion to complete (smiles → molblock on 50 rows).
    await page.waitForTimeout(15_000);
  });

  await softStep('Exactly-once on original column post-commit', async () => {
    const onOriginal = await page.evaluate(() => {
      const molColName = (window as any).__grok17964_origMolCol;
      grok.shell.o = grok.shell.t.col(molColName);
      return new Promise(resolve => {
        setTimeout(() => {
          const entries = Array.from(document.querySelectorAll('label.d4-link-action'))
            .filter(l => (l.textContent ?? '').trim().startsWith('Convert Notation'));
          resolve({count: entries.length});
        }, 1800);
      });
    }) as {count: number};
    expect(
      onOriginal.count,
      `GROK-17964 regression: registration count on ORIGINAL column post-commit expected 1, got ${onOriginal.count}.`,
    ).toBe(1);
  });

  await softStep('Multi-invocation hardening: open + CANCEL twice on original column, recount', async () => {
    // Refocus original column + verify link is present before starting
    const ready = await page.evaluate(async () => {
      const molColName = (window as any).__grok17964_origMolCol;
      grok.shell.o = grok.shell.t.col(molColName);
      await new Promise(r => setTimeout(r, 2500));
      const link = Array.from(document.querySelectorAll('label.d4-link-action'))
        .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation'));
      return {present: !!link};
    });
    expect(ready.present, 'Pre-multi-invocation: Convert Notation link not visible after re-focus on original column').toBe(true);
    for (let i = 0; i < 2; i++) {
      await page.evaluate(async () => {
        const link = Array.from(document.querySelectorAll('label.d4-link-action'))
          .find(l => (l.textContent ?? '').trim().startsWith('Convert Notation')) as HTMLElement;
        if (!link) throw new Error('Convert Notation link not found during multi-cancel');
        link.click();
        await new Promise(r => setTimeout(r, 1500));
      });
      await page.locator('.d4-dialog').waitFor({timeout: 8000});
      await page.locator('.d4-dialog [name="button-CANCEL"]').click();
      await page.waitForTimeout(1500);
    }
    const afterMulti = await page.evaluate(() => {
      const molColName = (window as any).__grok17964_origMolCol;
      grok.shell.o = grok.shell.t.col(molColName);
      return new Promise(resolve => {
        setTimeout(() => {
          const entries = Array.from(document.querySelectorAll('label.d4-link-action'))
            .filter(l => (l.textContent ?? '').trim().startsWith('Convert Notation'));
          resolve({count: entries.length});
        }, 1800);
      });
    }) as {count: number};
    expect(
      afterMulti.count,
      `GROK-17964 regression: registration count after multi-cancel expected 1, got ${afterMulti.count}.`,
    ).toBe(1);
  });

  await softStep('Global final assertion: only one panel-attached Convert Notation entry visible', async () => {
    // Active Context Panel renders a single column's Actions pane; per-column
    // entries do not stack. Hidden `.d4-menu-item-label` "Convert Notation..." in
    // the menu fragment (rect.w === 0) is excluded by the `label.d4-link-action`
    // selector. Global count = panel-attached count.
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

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
