// GROK-17964: Convert Notation column-action must register exactly once across cancel/commit/repeat invocations.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

test('Chem: GROK-17964 Convert Notation column-action registration is exactly-once', async ({page}) => {
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
    await page.waitForTimeout(2000);
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
        targetSelect.value = 'molblock';
        targetSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
    });
    await page.locator('.d4-dialog [name="button-OK"]').click();
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
