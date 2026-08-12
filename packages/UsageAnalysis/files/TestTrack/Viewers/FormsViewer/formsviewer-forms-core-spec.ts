/* ---
realizes: [formsviewer.cp.forms-core, formsviewer.int.sort-mirrors-grid, formsviewer.int.pinned-rows-persist-by-value, formsviewer.int.selection-intersects-filter, formsviewer.edge.pin-non-unique-value-warns, formsviewer.edge.pinned-row-absent-from-ordinary-cards]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {knownOpenBug} from '../../helpers/known-open-bug';
import * as projects from '../../helpers/projects';
import {
  HOST, ORDINARY, CURRENT, PINNED, PINNED_PANE,
  cardFieldValue, cardIndexByValue, balloonCount, drawnLabelNames, waitForOrderStable,
  fieldValuesByPosition, withConsoleErrorCount,
} from '../../helpers/forms';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Both reads slice by the leading offset FIRST, then drop the nulls (order matters).
async function tailUsubjids(page: Page, offset: number): Promise<string[]> {
  return (await fieldValuesByPosition(page, 'USUBJID')).slice(offset).filter((x): x is string => x !== null);
}

async function ordinaryUsubjids(page: Page): Promise<string[]> {
  return (await fieldValuesByPosition(page, 'USUBJID')).filter((x): x is string => x !== null);
}

async function ordinaryHeights(page: Page, offset: number): Promise<number[]> {
  return page.evaluate(({sel, off}) => Array.from(document.querySelectorAll(sel)).slice(off)
    .map((c) => parseFloat(((c as HTMLElement).querySelector('[column="HEIGHT"]') as HTMLInputElement)?.value))
    .filter((x) => !Number.isNaN(x)), {sel: ORDINARY, off: offset});
}

async function sortIndicatorLabels(page: Page): Promise<string[]> {
  return page.evaluate((host) => Array.from(document.querySelectorAll(`${host} .d4-multi-form-header .d4-multi-form-column-name`))
    .filter((l) => l.querySelector('.d4-multi-form-column-sort-indicator'))
    .map((l) => (l.querySelector('div[name^="div-"]') as HTMLElement)?.getAttribute('name') ?? '')
    .filter((n) => n.length > 0), HOST);
}

async function sortArrow(page: Page, column: string): Promise<string | null> {
  return page.evaluate(({host, col}) => {
    const label = Array.from(document.querySelectorAll(`${host} .d4-multi-form-header .d4-multi-form-column-name`))
      .find((l) => l.querySelector(`div[name="div-${col}"]`));
    const ind = label?.querySelector('.d4-multi-form-column-sort-indicator');
    return ind ? (ind.textContent ?? '').trim() : null;
  }, {host: HOST, col: column});
}

/**
 * Right-click a card (real trusted gesture), then invoke a top-level context-menu item.
 * `columnField` lands the right-click on that field's `[column]` element, load-bearing for Pin Row:
 * the viewer records the pin as (column, stringified value) from `e.target.closest('[column]')`. A
 * centre click would pin whatever field sits at the midpoint (HEIGHT, a non-unique float); USUBJID
 * gives a unique key the by-value layout restore resolves back to the row.
 */
async function cardContextMenu(
  page: Page, cardSelector: string, cardIndex: number, itemName: string, columnField?: string,
): Promise<void> {
  const card = page.locator(cardSelector).nth(cardIndex);
  const target = columnField ? card.locator(`[column="${columnField}"]`).first() : card;
  await target.click({button: 'right'});
  await page.locator(`[name="${itemName}"]`).first().waitFor({timeout: 5000});
  await page.locator(`[name="${itemName}"]`).first().click();
}

/**
 * Sustained zero-balloon check: a unique-value pin must raise no warning — the negative control for
 * 6c's positive warning. Mirrors 6c's poll budget so a late balloon has the same window to appear.
 */
async function expectNoBalloonSustained(page: Page, windowMs = 15_000): Promise<void> {
  const deadline = Date.now() + windowMs;
  let seen = 0;
  do {
    seen = await balloonCount(page);
    if (seen > 0) break;
    await page.waitForTimeout(500);
  } while (Date.now() < deadline);
  expect(seen).toBe(0);
}

test('Forms viewer — core ladder (p0)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // #### Step 1
  await softStep('Step 1 — Add the Forms viewer; default field set is the first 20 visible columns', async () => {
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    const expectedFields = await page.evaluate(() =>
      grok.shell.t.columns.names().filter((n: string) => !n.startsWith('~')).slice(0, 20));

    const drawn = await drawnLabelNames(page);
    expect(drawn).toEqual(expectedFields);
    expect(drawn.some((n) => n.startsWith('~'))).toBe(false);
    expect(await balloonCount(page)).toBe(0);
  });

  // #### Step 2
  await softStep('Step 2 — The current-row card shows the grid value and follows the current row', async () => {
    // A fresh headless mount has no current row (idx -1) → the current-row card builds empty and
    // grid.cell('HEIGHT', -1) reads "". Establish a row first. CURRENT selects by green indicator.
    const startRow = 12;
    await page.evaluate((r) => { grok.shell.t.currentRowIdx = r; }, startRow);
    const gridStart = await page.evaluate((r) => grok.shell.tv.grid.cell('HEIGHT', r).cell.valueString, startRow);
    await expect.poll(() => cardFieldValue(page, 0, 'HEIGHT', CURRENT), {timeout: 10_000}).toBe(gridStart);

    await page.evaluate(() => { grok.shell.t.currentRowIdx = 77; });
    const grid77 = await page.evaluate(() => grok.shell.tv.grid.cell('HEIGHT', 77).cell.valueString);
    await expect.poll(() => cardFieldValue(page, 0, 'HEIGHT', CURRENT), {timeout: 10_000}).toBe(grid77);
  });

  // #### Step 3
  await softStep('Step 3 — Show Selected Rows renders one card per selected row beyond the two leading cards', async () => {
    // Confirm the ship default BEFORE any write — a set-then-read would be a tautology.
    const defaultOn = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return vw.props.showSelectedRows;
    });
    expect(defaultOn).toBe(true);

    // Three selected rows with mixed SEX so Step 4's filter can exclude one; a current row keeps the
    // leading offset a stable 2.
    const picked = await page.evaluate(() => {
      const df = grok.shell.t;
      df.currentRowIdx = 0;
      let f = -1; let m = -1; let m2 = -1;
      for (let i = 0; i < df.rowCount; i++) {
        const s = df.col('SEX').get(i);
        if (s === 'F' && f < 0) f = i;
        else if (s === 'M' && m < 0) m = i;
        else if (s === 'M' && m2 < 0 && i !== m) m2 = i;
        if (f >= 0 && m >= 0 && m2 >= 0) break;
      }
      df.selection.setAll(false);
      df.selection.set(f, true); df.selection.set(m, true); df.selection.set(m2, true);
      return {f, m, m2, usubjids: [f, m, m2].map((r) => df.col('USUBJID').get(r))};
    });

    // Leading offset 2 = current-row card + the always-built (headless-empty) mouse-over card, so
    // three selected rows make five card elements.
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(5);
    const extras = await tailUsubjids(page, 2);
    expect(extras.length).toBe(3);
    for (const usub of extras) expect(picked.usubjids).toContain(usub);
  });

  // #### Step 4
  await softStep('Step 4 — The selected-row cards equal selection ∩ filter after a filter change', async () => {
    await v.applyCategoricalFilter(page, 'SEX', ['M']);

    // Recompute selection ∩ filter INSIDE the poll (same evaluate as the DOM tail), so the expected
    // set is never captured before the 50 ms-debounced re-render lands. `excludedASelectedRow`
    // requires the filter to actually drop a selected row — the load-bearing claim, kept as an
    // invariant not a hardcoded count.
    await expect.poll(() => page.evaluate((sel) => {
      const df = grok.shell.t;
      const inter: string[] = [];
      let selCount = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (df.selection.get(i)) selCount++;
        if (df.selection.get(i) && df.filter.get(i)) inter.push(df.col('USUBJID').get(i));
      }
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({
        excludedASelectedRow: inter.length < selCount,
        match: JSON.stringify(inter) === JSON.stringify(tail),
      });
    }, ORDINARY), {timeout: 20_000})
      .toBe(JSON.stringify({excludedASelectedRow: true, match: true}));

    await v.resetFilters(page);
    await expect.poll(async () => (await tailUsubjids(page, 2)).length, {timeout: 15_000}).toBe(3);
  });

  // #### Step 5a
  await softStep('Step 5a — Sorting the grid by HEIGHT mirrors the card order and marks the HEIGHT label', async () => {
    await page.evaluate(() => grok.shell.tv.grid.sort(['HEIGHT'], [true]));
    // Card order must EQUAL df.getSortedOrder(...) at runtime, NOT mere monotonicity of three HEIGHT
    // values (passes by luck ~1/6). Expected and DOM tail read in one evaluate; `nonEmpty` guards a
    // vacuous empty-vs-empty match.
    await expect.poll(() => page.evaluate((sel) => {
      const df = grok.shell.t;
      const selFilter = df.selection.clone().and(df.filter);
      const order = df.getSortedOrder(['HEIGHT'], [true], selFilter);
      const expected = Array.from(order).map((r: number) => df.col('USUBJID').get(r));
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({nonEmpty: expected.length > 0, match: JSON.stringify(expected) === JSON.stringify(tail)});
    }, ORDINARY), {timeout: 20_000}).toBe(JSON.stringify({nonEmpty: true, match: true}));
    expect(await sortArrow(page, 'HEIGHT')).not.toBeNull();
  });

  // #### Step 5b
  await softStep('Step 5b — sortByColumnName overrides the grid sort; the indicator moves to WEIGHT', async () => {
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: 'WEIGHT'});
    });
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(['div-WEIGHT']);

    const weightArrow = await sortArrow(page, 'WEIGHT');
    expect(weightArrow).not.toBeNull();
    await expect.poll(() => page.evaluate(({sel, asc}) => {
      const df = grok.shell.t;
      const selFilter = df.selection.clone().and(df.filter);
      const order = df.getSortedOrder(['WEIGHT'], [asc], selFilter);
      const expected = (Array.from(order).map((r: number) => df.col('USUBJID').get(r))) as string[];
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({nonEmpty: expected.length > 0, match: JSON.stringify(expected) === JSON.stringify(tail)});
    }, {sel: ORDINARY, asc: weightArrow === '↑'}), {timeout: 20_000}).toBe(JSON.stringify({nonEmpty: true, match: true}));
    expect(await page.evaluate(() => grok.shell.tv.grid.sortByColumns.map((c: any) => c.name)))
      .toEqual(['HEIGHT']);
  });

  // #### Step 5c — GROK-20380 (OPEN), known-red
  await softStep('Step 5c — Turning Use Grid Sort OFF stops mirroring the grid sort (GROK-20380 known-red)', async () => {
    // Clear the viewer's own sort so the grid mirror is the ordering authority.
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: null});
    });
    // READINESS BARRIER, not the order assertion: waits for the grid HEIGHT sort to land in the card
    // order before the useGridSort toggle flips below. The real order check is exact getSortedOrder
    // equality (5a/5b/5c); this monotonicity is only a settle gate, never a verdict channel.
    await expect.poll(async () => {
      const heights = await ordinaryHeights(page, 2);
      return heights.length >= 2 && heights.every((h, i, a) => i === 0 || a[i - 1] <= h);
    }, {timeout: 20_000}).toBe(true);

    // Turn Use Grid Sort OFF through the property panel (Context Panel > Misc), not vw.props. The
    // category can ship collapsed, so expand it first, then flip the checkbox.
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'use-grid-sort');
    await v.setPropertyGridCheckbox(page, 'use-grid-sort', false, 'misc');
    await page.waitForTimeout(3500);

    // With useGridSort OFF the card order must NOT mirror the grid sort — exact sequence equality
    // against df.getSortedOrder(grid.sortByColumns, ...), never weak monotonicity (could pass by
    // accident on three rows and lose the self-flip signal). Open → strict assert fails, logged
    // green; fixed → wrapper throws [KNOWN_BUG_FIXED].
    const mirror = JSON.parse(await page.evaluate((sel) => {
      const df = grok.shell.t;
      const grid = grok.shell.tv.grid;
      const selFilter = df.selection.clone().and(df.filter);
      const order = df.getSortedOrder(grid.sortByColumns.map((c: any) => c.name), grid.sortTypes, selFilter);
      const mirrored = (Array.from(order).map((r: number) => df.col('USUBJID').get(r))) as string[];
      const tail = (Array.from(document.querySelectorAll(sel)).slice(2)
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value ?? null)
        .filter((x) => x !== null)) as string[];
      return JSON.stringify({expectedLen: mirrored.length, tailLen: tail.length,
        mirrors: JSON.stringify(mirrored) === JSON.stringify(tail)});
    }, ORDINARY)) as {expectedLen: number; tailLen: number; mirrors: boolean};
    // Non-emptiness is a HARD precondition asserted OUTSIDE the wrapper: an empty compared set must
    // fail LOUDLY, not fold into `mirrors === false` and let the wrapper report FIXED on a broken
    // build. Both sides stay non-empty in either bug state, so this cannot false-fail on fix day.
    expect(mirror.expectedLen).toBeGreaterThan(0);
    expect(mirror.tailLen).toBeGreaterThan(0);
    await knownOpenBug('GROK-20380', () => { expect(mirror.mirrors).toBe(false); });
  });

  // #### Step 5d
  await softStep('Step 5d — Double-clicking the sort label cycles the indicator; a different label does not move it', async () => {
    // Re-enable Use Grid Sort through the panel (symmetric with the Step 5c OFF
    // toggle above — the prose frames the reset as a panel action too). Sort By is
    // the column-selector string property, addressed as sortByColumnName; it stays
    // a setOptions write (not a boolean, and no boolean-panel control applies).
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'use-grid-sort');
    await v.setPropertyGridCheckbox(page, 'use-grid-sort', true, 'misc');
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: 'AGE'});
    });
    await expect.poll(() => sortArrow(page, 'AGE'), {timeout: 20_000}).not.toBeNull();

    const ageLabel = page.locator(`${HOST} .d4-multi-form-header [name="div-AGE"]`).first();
    // Three double-clicks produce three DISTINCT indicator states (one the cleared/null) — NOT an
    // exact order: entry direction (sortAscending) is uncontrolled here, so a fixed arrow sequence
    // would go red on the incoming state. Do not "strengthen" this to a fixed sequence.
    const seq: (string | null)[] = [];
    for (let i = 0; i < 3; i++) {
      await ageLabel.dblclick();
      await page.waitForTimeout(1500);
      seq.push(await sortArrow(page, 'AGE'));
    }
    expect(new Set(seq).size).toBe(3);
    expect(seq).toContain(null);
    expect(seq.filter((s) => s !== null).length).toBe(2);

    // The "different-column dblclick leaves the indicator on the sort column" claim is CONDITIONAL on
    // entry direction: the ondblclick handler (forms-viewer.ts:268-277) branches only on
    // sortByColumnName set/unset and sortAscending, never on WHICH label. From ascending a dblclick
    // on any label clears the sort; the claim holds ONLY from DESCENDING. Establish it
    // deterministically: clear the viewer's sort, then a single AGE dblclick always lands descending.
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: null});
    });
    await page.waitForTimeout(1500);
    await ageLabel.dblclick();
    await page.waitForTimeout(1500);
    await expect.poll(() => sortArrow(page, 'AGE'), {timeout: 15_000}).toBe('↓');
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 15_000}).toEqual(['div-AGE']);

    // From descending, a dblclick on a DIFFERENT column's label leaves the sort (and indicator) on
    // AGE — the descending branch only flips AGE to ascending, never reassigns or clears.
    await page.locator(`${HOST} .d4-multi-form-header [name="div-HEIGHT"]`).first().dblclick();
    await page.waitForTimeout(1500);
    expect(await sortIndicatorLabels(page)).toEqual(['div-AGE']);
  });

  // #### Step 6a
  await softStep('Step 6a — Pin Row moves the card to the pinned pane and removes it from the ordinary set', async () => {
    // Show Mouse Over Row OFF via the property panel (real UI actuation, not a props write). Puts the
    // leading offset at 1 and keeps the pinned row from re-showing at the mouse-over position; a
    // current row outside the selection keeps that single leading card populated.
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', false, 'misc');
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: null});
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = 0;
      df.selection.setAll(false);
      df.selection.set(5, true); df.selection.set(10, true); df.selection.set(20, true);
    });
    // Offset 1 + three selected rows = 4 ordinary cards.
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(4);
    await waitForOrderStable(page);

    // Pin target = first selected-row card, addressed by POSITION in the full positional array
    // (leading empties kept as null; offset 1 here). Never reuse a null-dropped value index as a DOM
    // position — an empty leading card collapses the list and misaligns it.
    const byPos6a = await fieldValuesByPosition(page, 'USUBJID');
    const targetPos6a = byPos6a.findIndex((val, i) => i >= 1 && val !== null);
    expect(targetPos6a).toBeGreaterThanOrEqual(0);
    const beforeCount = await page.locator(ORDINARY).count();
    const anchor = await cardFieldValue(page, targetPos6a, 'USUBJID');
    expect(anchor).not.toBeNull();

    // Right-click the USUBJID field so the pin is recorded by a unique key — see cardContextMenu.
    await cardContextMenu(page, ORDINARY, targetPos6a, 'div-Pin-Row', 'USUBJID');
    await expectNoBalloonSustained(page);

    await expect.poll(async () =>
      page.evaluate((sel) => getComputedStyle(document.querySelector(sel) as HTMLElement).display, PINNED_PANE),
    {timeout: 15_000}).not.toBe('none');
    const pinnedValues = await page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
      .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED);
    expect(pinnedValues).toEqual([anchor]);

    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(beforeCount - 1);
    expect(await ordinaryUsubjids(page)).not.toContain(anchor);
    expect(await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return vw.props.pinnedRowValues;
    })).toEqual([anchor]);
    // The pinned row is STILL selected — otherwise "absent although still selected" is no claim.
    expect(await page.evaluate((usub) => {
      const df = grok.shell.t;
      for (let i = 0; i < df.rowCount; i++)
        if (df.col('USUBJID').get(i) === usub) return df.selection.get(i);
      return false;
    }, anchor)).toBe(true);
  });

  // #### Step 6b
  await softStep('Step 6b — Unpin Row returns the row to the ordinary set and hides the pinned pane', async () => {
    const anchor = await page.evaluate((sel) =>
      ((document.querySelector(sel) as HTMLElement)?.querySelector('[column="USUBJID"]') as HTMLInputElement)?.value,
    PINNED);
    const beforeCount = await page.locator(ORDINARY).count();

    await cardContextMenu(page, PINNED, 0, 'div-Unpin-Row');

    await expect.poll(async () =>
      page.evaluate((sel) => getComputedStyle(document.querySelector(sel) as HTMLElement).display, PINNED_PANE),
    {timeout: 15_000}).toBe('none');
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(beforeCount + 1);
    expect(await ordinaryUsubjids(page)).toContain(anchor);

    // Re-pin so the persistence peak carries a pinned row.
    await page.evaluate(() => {
      const df = grok.shell.t;
      df.selection.setAll(false);
      df.selection.set(5, true); df.selection.set(10, true); df.selection.set(20, true);
    });
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(4);
    // Re-pin the first selected-row card by POSITION in the full positional array (offset 1). Do not
    // derive the DOM position from a null-dropped value index.
    const byPos6b = await fieldValuesByPosition(page, 'USUBJID');
    const rePinPos = byPos6b.findIndex((val, i) => i >= 1 && val !== null);
    expect(rePinPos).toBeGreaterThanOrEqual(0);
    await cardContextMenu(page, ORDINARY, rePinPos, 'div-Pin-Row', 'USUBJID');
    await expectNoBalloonSustained(page);
    await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel)).length, PINNED),
      {timeout: 15_000}).toBe(1);
  });

  // #### Step 6c — formsviewer.edge.pin-non-unique-value-warns
  await softStep('Step 6c — Pinning through a NON-UNIQUE field raises the exact warning; the single pin is preserved', async () => {
    // 6a/6b leave one USUBJID pin (unique value → no warning). This step adds a SECOND pin through a
    // card's SEX field — a value repeating earlier in the column — so findRowByValue resolves to an
    // earlier row and pinRow warns. The SEX pin is removed at the end, leaving the USUBJID pin.
    const pinnedBefore = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return Array.from(vw.props.pinnedRowValues as string[]);
    });
    expect(pinnedBefore.length).toBe(1);

    // A row whose SEX value is NOT the first occurrence (so the pin is non-unique). Select it, keep
    // the current row off it, so it gets its own ordinary card to right-click.
    const target = await page.evaluate(() => {
      const df = grok.shell.t;
      const sex = df.col('SEX');
      const firstOf: Record<string, number> = {};
      for (let i = 0; i < df.rowCount; i++) {
        const val = String(sex.get(i));
        if (!(val in firstOf)) firstOf[val] = i;
      }
      let row = -1;
      for (let i = 0; i < df.rowCount; i++)
        if (firstOf[String(sex.get(i))] !== i) { row = i; break; }
      df.currentRowIdx = 0;
      if (row >= 0) df.selection.set(row, true);
      return {row, usub: df.col('USUBJID').get(row), sex: String(sex.get(row))};
    });
    expect(target.row).toBeGreaterThanOrEqual(0);

    await waitForOrderStable(page);
    const cardIdx = await cardIndexByValue(page, 'USUBJID', target.usub);
    expect(cardIdx).toBeGreaterThanOrEqual(0);

    // Right-click the SEX field so the pin records the non-unique SEX value.
    await cardContextMenu(page, ORDINARY, cardIdx, 'div-Pin-Row', 'SEX');

    // The exact, product-computed warning text (pinRow, forms-viewer.ts#L561).
    await expect.poll(() => page.evaluate(() =>
      document.querySelector('.d4-balloon.warning .d4-balloon-content')?.textContent ?? null),
    {timeout: 15_000}).toBe("You have pinned a non-unique value. It won't be applied from the layout.");

    // The pin still happened for the session — pinnedRowValues now also carries the SEX value.
    expect(await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return Array.from(vw.props.pinnedRowValues as string[]);
    })).toContain(target.sex);

    // Restore: unpin the SEX row (found by value in the pinned pane), leaving the USUBJID pin
    // Scenario 7 relies on.
    const pinnedSexIdx = await cardIndexByValue(page, 'USUBJID', target.usub, PINNED);
    expect(pinnedSexIdx).toBeGreaterThanOrEqual(0);
    await cardContextMenu(page, PINNED, pinnedSexIdx, 'div-Unpin-Row');
    await expect.poll(() => page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return Array.from(vw.props.pinnedRowValues as string[]);
    }), {timeout: 15_000}).toEqual(pinnedBefore);
  });

  // #### Step 7a / Step 7b
  await softStep('Step 7a / Step 7b — Re-applying the saved layout over a deliberately corrupted view restores the field set, sort-label identity and pinned row by value, and drops a foreign viewer not in the layout', async () => {
    // Move the field set OFF its default BEFORE saving, so 7a/7c prove a REAL round-trip not a
    // default-vs-default compare (a build ignoring the saved set and rebuilding the default would
    // pass both sides). Chosen set is a proper, reordered subset: USUBJID (pin anchor + by-value key)
    // and AGE (sort indicator) stay; the rest trimmed and reversed so DRAWN ORDER also differs.
    const chosenFields = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      const all = Array.from(vw.props.fieldsColumnNames as string[]);
      const must = ['USUBJID', 'AGE'];
      const rest = all.filter((n) => !must.includes(n));
      const picked = [...must, ...rest.slice(0, Math.max(1, rest.length - 2))].reverse();
      vw.setOptions({fieldsColumnNames: picked});
      return picked;
    });
    await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(chosenFields);

    // Establish an accumulated sort so its label identity can be checked after restore.
    await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      vw.setOptions({sortByColumnName: 'AGE'});
    });
    await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(['div-AGE']);

    // Baseline is the DRAWN label sequence; the property is an additional channel, not the sole one.
    const preLabels = await drawnLabelNames(page);
    const pre = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return {
        fields: Array.from(vw.props.fieldsColumnNames as string[]),
        pinnedValues: Array.from(vw.props.pinnedRowValues as string[]),
      };
    });
    const preIndicator = await sortIndicatorLabels(page);
    // Assert every baseline non-empty AT CAPTURE: the softStep shell keeps running past a failure, so
    // an empty baseline must not degrade the post-restore compare into empty-vs-empty.
    expect(preLabels.length).toBeGreaterThan(0);
    expect(pre.pinnedValues.length).toBeGreaterThan(0);
    expect(preIndicator.length).toBeGreaterThan(0);

    // Save through the REAL user route the atlas names — View | Layout | Save to Gallery — not
    // programmatic grok.dapi.layouts.save. The save is SILENT (auto-named), so locate the saved
    // layout deterministically by diffing the table's applicable-layout ids before/after: each run
    // claims only its own fresh id, any leaked earlier layout is already in the before-set.
    const currentUserId = await page.evaluate(() => String(grok.shell.user.id));
    const applicableLayouts = async (): Promise<{id: string; authorId: string | null}[]> =>
      page.evaluate(async () => ((await grok.dapi.layouts.getApplicable(grok.shell.t)) ?? [])
        .map((l: any) => ({id: String(l.id), authorId: l.author && l.author.id ? String(l.author.id) : null})));
    const beforeIds = (await applicableLayouts()).map((l) => l.id);

    // Drive the leaf with the REAL mouse; the menu mechanic (real pointer move to expand the subgroup
    // — a synthetic mouseover does NOT — size-not-DOM readiness, bundle retry) lives in
    // v.driveTopMenuLeaf. Keep only the hard actuation assertion; NO programmatic-save fallback.
    expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery'])).toBe(true);

    // The silent save reaches the gallery a second or two later. Require EXACTLY ONE new
    // applicable-layout id; if more appear (a concurrent save on this shared server), fail LOUDLY
    // and delete NOTHING rather than guess which is ours and delete a stranger's entity.
    let fresh: string[] = [];
    await expect.poll(async () => {
      // Author filter first, then exactly-one — never apply/delete a stranger's concurrent layout as ours.
      fresh = (await applicableLayouts())
        .filter((l) => !beforeIds.includes(l.id) && l.authorId === currentUserId).map((l) => l.id);
      return fresh.length;
    }, {timeout: 20_000, intervals: [500, 1000, 2000, 3000]}).toBeGreaterThanOrEqual(1);
    if (fresh.length !== 1) {
      throw new Error(
        `Layout save produced ${fresh.length} new applicable layouts authored by the current user ` +
        `(${fresh.join(', ')}), expected exactly 1 — refusing to guess which is ours; deleting none. ` +
        `A concurrent layout creation on the same dataset is the likely cause.`);
    }
    const layoutId = fresh[0];

    try {
      // Corrupt the current view rather than closeAll + reopen: close Forms and add a FOREIGN
      // Histogram the saved layout lacks. Applying the layout must then BOTH bring Forms back AND
      // drop the Histogram — the "layout removes a widget not in it" half a fresh reopen can't prove.
      await page.evaluate(() => {
        const tv = grok.shell.tv;
        const forms = tv.viewers.find((x: any) => x.type === 'FormsViewer');
        if (forms) forms.close();
        tv.addViewer('Histogram');
      });
      // Confirm the corruption landed before applying, else "Forms returned / Histogram gone" is vacuous.
      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 20_000}).toBe(1);
      await expect.poll(() => page.locator(HOST).count(), {timeout: 20_000}).toBe(0);

      // Wrap the layout re-apply (the section's riskiest op) in the zero-console-error floor: it must
      // not emit a single console.error. Texts collected for diagnostics only; threshold stays zero.
      const applyErrTexts: string[] = [];
      const applyErrCount = await withConsoleErrorCount(page, async () => {
        await page.evaluate(async (id) => {
          const saved = await grok.dapi.layouts.find(id);
          grok.shell.tv.loadLayout(saved);
        }, layoutId);
        await page.locator(HOST).first().waitFor({timeout: 30_000});
      }, undefined, applyErrTexts);
      expect(applyErrCount, `layout-apply console errors: ${JSON.stringify(applyErrTexts)}`).toBe(0);

      // Applying the layout drops the foreign Histogram it does not contain — Step 7a's second claim.
      await expect.poll(() => page.locator('[name="viewer-Histogram"]').count(), {timeout: 20_000}).toBe(0);
      expect(await page.evaluate(() => grok.shell.tv.viewers
        .filter((x: any) => x.type === 'Histogram').length)).toBe(0);
      await expect.poll(() => drawnLabelNames(page), {timeout: 20_000}).toEqual(preLabels);
      await expect.poll(() => page.evaluate(() => {
        const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
        return vw ? Array.from(vw.props.fieldsColumnNames as string[]) : null;
      }), {timeout: 20_000}).toEqual(pre.fields);

      // Label identity only, NOT direction — GROK-20666.
      await expect.poll(() => sortIndicatorLabels(page), {timeout: 20_000}).toEqual(preIndicator);

      // Restored BY VALUE: resolvePinnedRows maps the persisted (column, value) pair back to a row.
      await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED),
      {timeout: 20_000}).toEqual(pre.pinnedValues);
    } finally {
      await page.evaluate(async (id) => {
        const saved = await grok.dapi.layouts.find(id);
        if (saved) await grok.dapi.layouts.delete(saved);
      }, layoutId);
    }
  });

  // #### Step 7c
  await softStep('Step 7c — A project round-trip preserves the field set and the pinned row across a session boundary', async () => {
    const preLabels = await drawnLabelNames(page);
    const pre = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      // Return plain JS arrays, never raw prop objects: a Dart-boxed value walks a deep proxy chain
      // and Playwright throws "Cannot serialize result: object reference chain is too long".
      return {
        fields: Array.from(vw.props.fieldsColumnNames as string[]),
        pinnedValues: Array.from(vw.props.pinnedRowValues as string[]),
      };
    });
    const preIndicator = await sortIndicatorLabels(page);
    // Guard every baseline non-empty at capture (same reason as Step 7a) against empty-vs-empty.
    expect(preLabels.length).toBeGreaterThan(0);
    expect(pre.fields.length).toBeGreaterThan(0);
    expect(pre.pinnedValues.length).toBeGreaterThan(0);
    expect(preIndicator.length).toBeGreaterThan(0);
    // Lowercase unique name per saveProjectViaUI: server stores a PascalCased `name`, lookup by
    // friendlyName still works.
    const projName = `zz-formsviewer-core-${Date.now()}`;
    let projectId: string | null = null;
    try {
      // Only the REAL ribbon Save captures the view layout: grok.dapi.projects.save assembles an
      // empty project that reopens as a bare Grid — a false green.
      const saved = await projects.saveProjectViaUI(page, projName);
      projectId = saved.projectId;

      // No console-error floor around the project round-trip: the ribbon save/publish path emits a
      // message-less Dart stack trace (ProjectMeta.publish → Project.upload → ProjectClient.save) in
      // ~40% of runs — permanent platform noise (operator ruling), so a zero floor would be untrue.
      // The layout half (7a) keeps its floor; that path is clean.
      await projects.reopenAndAssertProvenance(page, projectId as string);
      await page.locator(HOST).first().waitFor({timeout: 30_000});

      await expect.poll(() => drawnLabelNames(page), {timeout: 30_000}).toEqual(preLabels);
      await expect.poll(() => page.evaluate(() => {
        const vw = grok.shell.tv?.viewers?.find((x: any) => x.type === 'FormsViewer');
        return vw ? Array.from(vw.props.fieldsColumnNames as string[]) : null;
      }), {timeout: 30_000}).toEqual(pre.fields);

      await expect.poll(() => page.evaluate((sel) => Array.from(document.querySelectorAll(sel))
        .map((c) => ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement)?.value), PINNED),
      {timeout: 30_000}).toEqual(pre.pinnedValues);

      // Sort-indicator label identity survives the session boundary too (not direction — GROK-20666).
      await expect.poll(() => sortIndicatorLabels(page), {timeout: 30_000}).toEqual(preIndicator);
    } finally {
      if (projectId)
        await projects.deleteProjectWithCleanup(page, {projectId});
    }
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
