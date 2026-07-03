import { test, expect } from './helpers';
import {
  clickTopMenuLeaf, inputHost, openDemoCsv, resetShell, selectAllColumnsInPicker,
  setInputValue, topDialog, waitForDialog,
} from './helpers';

// Test Track scenario: EDA/pls.md
// 1. Open cars.csv from Demo files.
// 2. Top Menu > ML > Analyze > PLS...
// 3. Select all available columns in Using; set Components = 3.
// 4. Click OK; expect PLS1, PLS2, PLS3 added.
//
// Scope limitation
// ----------------
// Step 4 cannot be automated end-to-end on dev today:
//   * The scenario's literal "select all available columns" puts the Predict column (`price`)
//     into Using, which the PLS validator silently rejects (RUN becomes disabled with only a
//     red border on Predict — no balloon, no tooltip; see pls-run.md retrospective).
//   * Working around the bug from the UI requires deselecting `price` from Using, but the
//     "Select columns..." sub-dialog uses a canvas-rendered grid with no DOM-accessible
//     per-row checkboxes (linear-regression-run.md notes this same gap for Train Model).
//   * `DG.Dialog.getOpenDialogs()` returns the analysis dialog handle with 0 inputs — the
//     inputs live entirely on the Dart side and are not exposed to JS.
//   * Calling `EDA:PLS` through the function registry (`functions.call` / `Func.apply` /
//     `Func.prepare`) strips the JS wrapper from the ColumnList `features` parameter on the
//     way through the Dart bridge, and the function body fails with
//     "TypeError: t.byIndex is not a function".
//
// The test therefore covers what the UI deterministically permits and stops short of
// asserting PLS1/PLS2/PLS3 columns are added. When the platform fixes the RUN-disable UX
// (auto-exclude Predict from Using, or surface the reason inline), the trailing
// assertion block can be uncommented to verify the full scenario.

test.describe.serial('EDA / PLS', () => {
  test.afterEach(async ({ page }) => { await resetShell(page); });

  test('PLS dialog opens on cars.csv with all expected inputs and accepts Components=3', async ({ page }) => {
    test.setTimeout(180_000);

    await openDemoCsv(page, 'cars.csv');

    await clickTopMenuLeaf(page, 'div-ML---Analyze---PLS...');
    await waitForDialog(page, 'PLS');

    // Dialog has the documented input set: Predict, Using, Components, Quadratic.
    const dlg = topDialog(page);
    for (const host of ['Predict', 'Using', 'Components', 'Quadratic']) {
      await expect(dlg.locator(inputHost(host)).first()).toBeVisible({ timeout: 5_000 });
    }

    // Using picker — clicking "All" selects every numeric column, mirroring the scenario step.
    await selectAllColumnsInPicker(page, 'Using');
    await expect(dlg.locator(`${inputHost('Using')} .ui-input-editor`)).toContainText(/all|\(\d+\)/i);

    await setInputValue(page, 'Components', '3');
    // The Components input reports its value through the editor's text/value attribute —
    // accept either to stay tolerant to the underlying input element type. Wait in-browser
    // (no per-attempt report errors), then assert once.
    await page.waitForFunction(() => {
      const d = Array.from(document.querySelectorAll('.d4-dialog')).slice(-1)[0];
      const el = d?.querySelector('[name="input-host-Components"] .ui-input-editor') as HTMLInputElement | null;
      return ((el?.value ?? el?.textContent?.trim() ?? '')) === '3';
    }, undefined, { timeout: 5_000 });
    const componentsValue = await dlg.locator(`${inputHost('Components')} .ui-input-editor`)
      .first().evaluate((el) => (el as HTMLInputElement).value ?? el.textContent?.trim() ?? '');
    expect(componentsValue).toBe('3');

    // Cancel the dialog cleanly; do not click RUN — see scope limitation in the file header.
    await dlg.locator('[name="button-CANCEL"]').click();
    await expect(page.locator('.d4-dialog .d4-dialog-title', { hasText: /^PLS$/i }))
      .toHaveCount(0, { timeout: 10_000 });
  });
});
