/* ---
realizes: [chem.int.empty-input-analyses]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu} from '../spec-login';
import {finishSpec, closeAllAndWait} from '../helpers/viewers';
import {knownOpenBug} from '../helpers/known-open-bug';
import * as chem from '../helpers/chem';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser chem.md):
//   [name="button-MCS"] — R-Groups Analysis dialog (title "R-Groups Analysis", buttons
//     MCS/OK/CANCEL); derives the MCS core from data. Observed live 2026-08-20 via chrome-devtools MCP.
//   .d4-balloon — the rejection message raised after OK (class-1, chem.md:783). The 2026-08-20
//     note here read "No R-Groups were found ~1s after OK once the MCS core has settled", with a
//     no-balloon path for an OK clicked before the settle. Measured on dev 2026-08-21 (own run,
//     4/4 attempts) supersedes that: on an all-empty Molecule column MCS derives NO core at all,
//     OK stays ENABLED, and confirming raises "No core was provided" — there is no core to settle,
//     so the "core settles" theory that originally justified the 4-attempt retry below no longer
//     describes the mechanism. What the retry is empirically worth: 4/4 observed attempts raised
//     the balloon, so the extra attempts have no measured value and are kept as unproven margin,
//     not as a settle-wait — a wait removed on a corrected theory alone is how flakes come back.
//     The load-bearing wait remains the balloon-DOM waitForFunction (rung-2), never a fixed sleep.
//   Chem Space dialog: title "Chem Space", buttons OK/CANCEL; on an all-empty Molecule column
//     Embed_X_1/Embed_Y_1 columns are appended with no rejection (GROK-18407, open). Observed live 2026-08-20 via MCP.

// Canonical grid-cell click path — chem.md:1938-1948 ("How to click a molecular grid cell — copy
// the green one, do not compute geometry"): the point is resolved by asking the grid's own hitTest,
// ported from Chem/chem-filters-spec.ts#L295 (Gate B green 3/3, cycle 2026-08-19-chem-filters-consolidation-01).
async function gridCellPoint(page: Page, column: string, gridRow: number):
  Promise<{x: number; y: number; tableRow: number} | null> {
  return page.evaluate(async ({col, row}) => {
    const grid = grok.shell.tv.grid;
    grid.scrollToCell(col, row);
    await new Promise((r) => setTimeout(r, 500));
    const mainGrid = [...document.querySelectorAll('[name="viewer-Grid"]')].find((g) => !g.closest('.d4-filter'));
    const overlay = mainGrid?.querySelector('[name="overlay"]') as HTMLElement | null;
    const gc = grid.columns.byName(col);
    if (!overlay || !gc) return null;
    const rect = overlay.getBoundingClientRect();
    const xLocal = Math.round((gc.left + gc.right) / 2);
    const band: number[] = [];
    let tableRow: number | null = null;
    for (let y = Math.round(grid.colHeaderHeight) + 2; y < rect.height; y += 4) {
      const hit = grid.hitTest(xLocal, y);
      if (hit && hit.gridRow === row && hit.gridColumn && hit.gridColumn.name === col && hit.isTableCell) {
        band.push(y);
        if (tableRow === null) tableRow = hit.tableRowIndex;
      }
      else if (band.length > 0)
        break;
    }
    if (band.length === 0 || tableRow === null) return null;
    const yLocal = band[Math.floor(band.length / 2)];
    return {x: Math.round(rect.left + xLocal), y: Math.round(rect.top + yLocal), tableRow};
  }, {col: column, row: gridRow});
}

async function openRGroupsDialog(page: Page): Promise<void> {
  await chem.openChemMenuItem(page, 'R-Groups Analysis...', {delayMs: 700});
  await page.locator('.d4-dialog:has([name="button-MCS"])').waitFor({timeout: 15000});
}

test('Chem: Empty-input boundary — R-Groups no-decomposition + Chemical Space empty-column', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // An all-empty Molecule column makes the decomposition empty BY CONSTRUCTION. Running
  // R-Groups on smiles.csv instead made the outcome depend on whether a core happens to
  // match that data — it does, so the balloon never appeared and the step failed on
  // correct product behaviour (cycles 2026-08-19-chem-new-07, 2026-08-20-chem-new-15).
  await softStep('Setup: open a table whose Molecule column is entirely empty', async () => {
    await page.evaluate(async () => {
      document.body.classList.add('selenium');
      try { (grok as any).shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
      try { (grok as any).shell.windows.simpleMode = true; } catch (e) {}
      grok.shell.closeAll();
      (window as any).__errors = [];
      const orig = console.error;
      console.error = function (...args: any[]) {
        (window as any).__errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
      let csv = 'id,structure\n';
      for (let i = 1; i <= 10; i++) csv += i + ',\n';
      const df = (DG as any).DataFrame.fromCsv(csv);
      df.col('structure').semType = 'Molecule';
      grok.shell.addTableView(df);
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await waitForChemMenu(page);
    const allEmpty = await page.evaluate(() => {
      const col = grok.shell.t.col('structure');
      return Array.from({length: grok.shell.t.rowCount}, (_: any, i: number) => col.get(i))
        .every((v: any) => v === '' || v == null);
    });
    expect(allEmpty,
      'the structure column must be entirely empty, or the decomposition is not guaranteed empty')
      .toBe(true);
  });

  let colsBefore: string[] = [];
  await softStep('S1.1: Chem → Analyze → R-Groups Analysis → dialog opens', async () => {
    colsBefore = await page.evaluate(() => grok.shell.t.columns.toList().map((c: any) => c.name));
    console.log(`[empty-input] S1.1 columns before the R-Groups run = ${JSON.stringify(colsBefore)}`);
    await openRGroupsDialog(page);
  });

  let balloonText: string | null = null;
  let okWasEnabled: boolean | null = null;
  await softStep('S1.2-3: choose MCS strategy (derive core from data), then OK to start decomposition', async () => {
    for (let attempt = 0; attempt < 4 && balloonText === null; attempt++) {
      if (attempt > 0)
        await openRGroupsDialog(page);
      await page.locator('.d4-dialog [name="button-MCS"]').click();
      await page.locator('.d4-dialog [name="button-OK"]').waitFor({state: 'visible', timeout: 5000});
      // Whether MCS can derive a core from an all-empty column is itself a product fact.
      // Record it so a missing balloon can be told apart from a disabled OK.
      okWasEnabled = await page.evaluate(() => {
        const ok = document.querySelector('.d4-dialog [name="button-OK"]') as HTMLElement | null;
        if (!ok) return null;
        if ((ok as HTMLButtonElement).disabled) return false;
        if (ok.classList.contains('disabled')) return false;
        return ok.classList.contains('enabled');
      });
      console.log(`[empty-input] attempt ${attempt}: MCS pressed, OK enabled = ${okWasEnabled}`);
      await page.locator('.d4-dialog [name="button-OK"]').click();
      // Record EVERY balloon the run raises, not only the one the step hopes for: which
      // message an all-empty column produces is the product fact this step is here to learn.
      const allBalloons = await page.evaluate(async () => {
        await new Promise((r) => setTimeout(r, 1200));
        return Array.from(document.querySelectorAll('.d4-balloon'))
          .map((b) => (b.textContent ?? '').trim()).filter(Boolean);
      });
      if (allBalloons.length > 0)
        console.log(`[empty-input] attempt ${attempt}: balloons = ${JSON.stringify(allBalloons)}`);
      // Read the MATCHING balloon's own text. A `[class*="balloon"]` sweep also matches the
      // container, whose textContent is every balloon glued together, so an unrelated
      // platform balloon ends up concatenated into the recorded evidence.
      const seen = await page.waitForFunction(() => {
        const hit = Array.from(document.querySelectorAll('.d4-balloon'))
          .map((b) => b.textContent?.trim() ?? '')
          .find((t) => /No core was provided|No R.?Groups were found/i.test(t));
        return hit ?? false;
      }, {timeout: 8000}).catch(() => null);
      if (seen)
        balloonText = String(await seen.jsonValue());
    }
  });

  await softStep('S1.4: "No core was provided" balloon; no null-ref crash; dataset unfiltered', async () => {
    const result = await page.evaluate(() => {
      const t = grok.shell.t;
      const errs = ((window as any).__errors ?? []) as string[];
      return {
        nullRefErr: errs.find((e) => /null.*reference|cannot read prop|cannot set prop/i.test(e)) ?? null,
        rowCount: t.rowCount,
        filterTrueCount: t.filter.trueCount,
      };
    });
    console.log(`[empty-input] S1.4 probe = ${JSON.stringify(result)}, balloon = ${JSON.stringify(balloonText)}, okEnabled = ${okWasEnabled}`);
    expect(result.nullRefErr,
      `GROK-16329 regression: null-reference crash on empty R-Groups decomposition. err=${result.nullRefErr}`)
      .toBeNull();
    expect(result.filterTrueCount,
      'dataset must remain fully unfiltered after an empty-result R-Groups run')
      .toBe(result.rowCount);
    // Observed 2026-08-21 (own run, 4/4 attempts): an all-empty column makes MCS derive no
    // core, OK stays ENABLED, and confirming raises "No core was provided" — the empty-RESULT
    // branch ("No R-Groups were found") is unreachable because rGroupDecomp short-circuits on
    // an empty core before decomposing. Both are clean rejections; this scenario is the
    // empty-INPUT boundary, so it asserts the one empty input actually produces.
    expect(balloonText,
      'empty input must be rejected with a visible message and no silent close. ' +
      `OK enabled after MCS on the all-empty column: ${okWasEnabled}`)
      .toMatch(/No core was provided/i);
  });

  await softStep('S1.5: no R-group columns appended; no trellis; app remains interactive', async () => {
    const result = await page.evaluate(() => {
      const t = grok.shell.t;
      return {
        cols: t.columns.toList().map((c: any) => c.name),
        rCols: t.columns.toList().map((c: any) => c.name).filter((n: string) => /^R\d+$/.test(n)),
        hasTrellis: Array.from(grok.shell.tv.viewers).some((v: any) => v.type === 'Trellis plot'),
      };
    });
    const appended = result.cols.filter((n: string) => !colsBefore.includes(n));
    console.log(`[empty-input] S1.5 probe = ${JSON.stringify(result)}, colsBefore = ${JSON.stringify(colsBefore)}, appended = ${JSON.stringify(appended)}`);
    // /^R\d+$/ below only sees the dialog's DEFAULT 'R' column prefix (chem.md:622), which this
    // spec neither sets nor asserts; the appended-set diff is prefix-independent, as in S2.
    expect(appended,
      'no result columns of any name may be appended when decomposition yields nothing').toEqual([]);
    expect(result.rCols,
      'no R-group result columns may be appended when decomposition yields nothing').toEqual([]);
    expect(result.hasTrellis,
      'no Trellis plot may be created for an empty R-Groups result').toBe(false);
    const gridInteractive = await page.locator('.d4-grid[name="viewer-Grid"]').first().isVisible();
    expect(gridInteractive, 'grid must remain visible/interactive after the failed run').toBe(true);

    // Visibility alone is satisfied by a frozen grid, so drive the scenario's own gesture. Grid row 3,
    // not 0: a fresh view can already sit on row 0, and an assertion a frozen grid passes proves nothing.
    const pt = await gridCellPoint(page, 'structure', 3);
    expect(pt, 'a structure cell in grid row 3 must be locatable through the grid\'s own hit test').not.toBeNull();
    const before = await page.evaluate(() =>
      ({row: grok.shell.t.currentRowIdx, col: grok.shell.t.currentCol?.name ?? null}));
    expect({row: before.row, col: before.col},
      'the current cell must not ALREADY be the one the click aims at, or the post-click read cannot tell ' +
      'a live grid from a frozen one')
      .not.toEqual({row: pt!.tableRow, col: 'structure'});

    await page.mouse.move(pt!.x, pt!.y, {steps: 4});
    await page.mouse.click(pt!.x, pt!.y);
    await page.waitForFunction((tr) => grok.shell.t.currentRowIdx === tr &&
      grok.shell.t.currentCol?.name === 'structure', pt!.tableRow, {timeout: 15000}).catch(() => {});
    // No cell tooltip is asserted here: this fixture has two columns and both are visible in the grid,
    // so grid_core.dart:357-359 filters every column out of the cell tooltip (showVisibleColumnsInTooltip
    // defaults to false, grid_look.dart:299), the row tooltip gets no rows, and tooltip.dart:449 returns
    // before _show — TOOLTIP_SHOWN (:365) never fires. Not showing it is the product's designed behaviour
    // on this table, so a tooltip assertion could only be made to pass by changing the fixture.
    const after = await page.evaluate(() => ({
      row: grok.shell.t.currentRowIdx,
      col: grok.shell.t.currentCol?.name ?? null,
    }));
    console.log(`[empty-input] S1.5 cell click: aimed grid row 3 = table row ${pt!.tableRow}, ` +
      `before = ${JSON.stringify(before)}, after = ${JSON.stringify(after)}`);
    expect(after.row,
      'the click must make current the very cell the grid\'s hit test resolved — a grid that stopped ' +
      `handling input leaves the current row at ${before.row}`).toBe(pt!.tableRow);
    expect(after.col,
      'the click must make the structure column current, not merely move the row').toBe('structure');
  });

  await softStep('S2.setup: open a table with an all-empty Molecule column', async () => {
    await closeAllAndWait(page);
    await page.evaluate(async () => {
      (window as any).__errors = [];
      // Balloons survive closeAllAndWait, so Scenario 1's linger into Scenario 2's readings.
      // Harmless while the load-bearing assert is embedCount, but a balloon-based rejection
      // assert here would be false-green from the first run.
      document.querySelectorAll('.d4-balloon').forEach((b) => b.remove());
      let csv = 'id,structure\n';
      for (let i = 1; i <= 10; i++) csv += i + ',\n';
      const df = (DG as any).DataFrame.fromCsv(csv);
      const col = df.col('structure');
      col.semType = 'Molecule';
      grok.shell.addTableView(df);
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});
    await waitForChemMenu(page);
    const allEmpty = await page.evaluate(() => {
      const col = grok.shell.t.col('structure');
      return Array.from({length: grok.shell.t.rowCount}, (_: any, i: number) => col.get(i))
        .every((v: any) => v === '' || v == null);
    });
    expect(allEmpty, 'the structure column must be entirely empty to exercise the empty-input boundary').toBe(true);
  });

  await softStep('S2.2: Chem → Analyze → Chemical Space → dialog opens with a title', async () => {
    await chem.openChemMenuItem(page, 'Chemical Space...', {delayMs: 700});
    await page.locator('.d4-dialog:has-text("Chem Space")').first().waitFor({timeout: 15000});
    const title = await page.evaluate(() =>
      document.querySelector('.d4-dialog .d4-dialog-header, .d4-dialog .d4-dialog-title')
        ?.textContent?.trim() ?? '');
    expect(title, `Chemical Space dialog title expected, got "${title}"`).toMatch(/Chem\s*Space/i);
  });

  await softStep('S2.3-5: Run on empty column must be REJECTED — no Embed_X column, no silent success', async () => {
    const embedBefore = await page.evaluate(() =>
      grok.shell.t.columns.toList().map((c: any) => c.name).filter((n: string) => /^Embed_X/i.test(n)).length);
    await page.locator('.d4-dialog [name="button-OK"]').click();

    await page.waitForFunction((before) => {
      const cols = grok.shell.t.columns.toList().map((c: any) => c.name);
      const embed = cols.filter((n: string) => /^Embed_X/i.test(n)).length;
      const balloons = Array.from(document.querySelectorAll('.d4-balloon-container .d4-balloon, .d4-balloon, [class*="balloon"]'))
        .map((b) => b.textContent?.trim() ?? '').filter(Boolean);
      const dlgOpen = !!document.querySelector('.d4-dialog');
      return embed !== before || balloons.length > 0 || !dlgOpen;
    }, embedBefore, {timeout: 20000}).catch(() => {});

    const outcome = await page.evaluate(() => {
      const cols = grok.shell.t.columns.toList().map((c: any) => c.name);
      return {
        embedCount: cols.filter((n: string) => /^Embed_X/i.test(n)).length,
        balloons: Array.from(document.querySelectorAll('.d4-balloon-container .d4-balloon, .d4-balloon, [class*="balloon"]'))
          .map((b) => b.textContent?.trim() ?? '').filter(Boolean),
        dlgOpen: !!document.querySelector('.d4-dialog'),
      };
    });

    const errs = await page.evaluate(() => ((window as any).__errors ?? []) as string[]);
    console.log(`[empty-input] S2.3-5 outcome = ${JSON.stringify(outcome)}, embedBefore = ${embedBefore}`);
    expect(errs.find((e) => /null.*reference|cannot read prop|uncaught/i.test(e)) ?? null,
      `no JS crash may fire on empty-column Chemical Space. errs=${JSON.stringify(errs.slice(0, 3))}`).toBeNull();

    await knownOpenBug('GROK-18407', () => {
      expect(outcome.embedCount,
        `GROK-18407: empty-column Chemical Space must be rejected — no Embed_X embedding column may be appended. outcome=${JSON.stringify(outcome)}`)
        .toBe(embedBefore);
    });
  });

  await page.evaluate(() => grok.shell.closeAll());
  finishSpec();
});
