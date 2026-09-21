import {expect, Page} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {settleContextPanes, waitQuiet} from './chem-fast-helpers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

// `grok.shell.o = column` does not stay put: the grid sets its own current cell as it finishes
// rendering, and the platform then makes THAT the current object — on smiles-50 the molregno
// cell wins and the panel fills with ChEMBL link panes instead of the column's. So the column is
// re-asserted until the panel is the column's.
const ACTION_PANE_CAP = 60_000;

// The link is what every assertion below reads, so when it does not arrive the state that
// decides it is reported with the failure: what the panel is showing, whether the action
// function is registered at all, and what the current object actually is.
async function expectConvertNotationLink(page: Page, where: string): Promise<void> {
  const present = async () => page.evaluate(() =>
    Array.from(document.querySelectorAll('label.d4-link-action'))
      .some((l) => (l.textContent ?? '').trim().startsWith('Convert Notation')));
  const reassert = async () => page.evaluate(() => {
    const name = (window as any).__grok17964_origMolCol;
    const col = name ? grok.shell.t?.col(name) : null;
    if (col && grok.shell.o !== col) { grok.shell.o = null; grok.shell.o = col; return true; }
    return false;
  });
  const deadline = Date.now() + ACTION_PANE_CAP;
  let reasserted = 0;
  while (Date.now() < deadline) {
    if (await present()) {
      if (reasserted) console.log(`[17964] ${where}: the column had to be re-asserted ${reasserted}x`);
      return;
    }
    // A re-asserted object rebuilds the accordion collapsed, and a collapsed pane holds no
    // content, so the panes are re-expanded on every turn of this loop.
    if (await reassert()) reasserted++;
    await page.evaluate(async () => {
      for (const p of Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane'))) {
        const h = p.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        if (h && !h.classList.contains('expanded')) { h.click(); await new Promise((r) => setTimeout(r, 100)); }
      }
    });
    await page.waitForTimeout(500);
  }
  const diag = await page.evaluate(() => {
    const o: any = grok.shell.o;
    const pane = document.querySelector('.grok-prop-panel .d4-accordion-pane[d4-title="Actions"]');
    return {
      currentObject: o ? `${o.constructor?.name} ${o.name ?? ''} semType=${o.semType ?? ''}` : null,
      contextPanelShown: grok.shell.windows.showContextPanel, simpleMode: grok.shell.windows.simpleMode,
      panes: Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane'))
        .map((x) => x.getAttribute('d4-title')),
      actionsPaneText: pane ? (pane as HTMLElement).innerText.replace(/\s+/g, ' ').slice(0, 300) : null,
      linkActions: Array.from(document.querySelectorAll('label.d4-link-action'))
        .map((l) => (l.textContent ?? '').trim()),
      actionFuncs: DG.Func.find({meta: {action: 'Convert Notation...'}}).map((f: any) => f.nqName),
    };
  });
  expect(false, `${where}: the Chem "Convert Notation..." column action never rendered on the Context ` +
    `Panel within ${ACTION_PANE_CAP} ms. state=${JSON.stringify(diag)}`).toBe(true);
}


test('Chem: GROK-17964 Convert Notation column-action registration is exactly-once', async ({page}) => {
  test.setTimeout(600_000);

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
      // with the context panel hidden `grok.shell.o = column` is ignored and no Actions pane is
      // built, which is how a neighbour that hides it starves the reads below
      try { grok.shell.windows.simpleMode = false; } catch (e) {}
      try { grok.shell.windows.showContextPanel = true; } catch (e) {}
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
    await settleContextPanes(page, 2000);

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

    await expectConvertNotationLink(page, 'after focusing the molecule column');
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
      const dlgDeadline = Date.now() + 1500;
      while (!document.querySelector('.d4-dialog') && Date.now() < dlgDeadline)
        await new Promise(r => setTimeout(r, 25));
    });
    await page.locator('.d4-dialog').waitFor({timeout: 8000});
    await page.locator('.d4-dialog [name="button-CANCEL"]').click();
    await waitQuiet(page.locator('.d4-dialog').waitFor({state: 'detached', timeout: 1500}));
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
      const dlgDeadline = Date.now() + 1500;
      while (!document.querySelector('[name="input-Target-Notation"]') && Date.now() < dlgDeadline)
        await new Promise(r => setTimeout(r, 25));
      const dlg = document.querySelector('.d4-dialog');
      const targetSelect = dlg?.querySelector('[name="input-Target-Notation"]') as HTMLSelectElement;
      if (targetSelect) {
        targetSelect.value = 'molblock';
        targetSelect.dispatchEvent(new Event('change', {bubbles: true}));
      }
    });
    const colsBeforeCommit = await page.evaluate(() => grok.shell.t.columns.length);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // The commit is done when the converted column is in the frame and the dialog is gone; a
    // conversion that never lands still burns the 15 s cap and fails the count below.
    await waitQuiet(page.waitForFunction((before: number) =>
      grok.shell.t.columns.length > before && !document.querySelector('.d4-dialog'),
    colsBeforeCommit, {timeout: 15_000, polling: 100}));
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
      await waitQuiet(page.locator('.d4-dialog').waitFor({state: 'detached', timeout: 1500}));
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
