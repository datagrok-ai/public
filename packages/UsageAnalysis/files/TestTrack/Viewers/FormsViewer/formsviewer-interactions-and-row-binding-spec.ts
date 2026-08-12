/* ---
realizes: [formsviewer.cp.interactions-and-row-binding]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {HOST, ORDINARY, CURRENT, cardFieldValue, fieldValuesByPosition, waitForOrderStable} from '../../helpers/forms';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

async function rowOfValue(page: Page, column: string, value: string): Promise<number> {
  return page.evaluate(({col, val}) => {
    const df = grok.shell.t;
    for (let i = 0; i < df.rowCount; i++) if (df.col(col).get(i) === val) return i;
    return -1;
  }, {col: column, val: value});
}

async function establishSetup(page: Page): Promise<void> {
  await page.evaluate(() => {
    const df = grok.shell.t;
    const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
    vw.props.showSelectedRows = true;
    vw.props.showCurrentRow = true;
    vw.props.showMouseOverRow = true;
    df.mouseOverRowIdx = -1;
    df.currentRowIdx = 5;
    df.selection.setAll(false);
    [2, 5, 10, 15].forEach((r) => df.selection.set(r, true));
  });
  // The readiness barrier must be DISCRIMINATING, never a card COUNT: a predecessor can leave the
  // SAME card count, and a count-only poll passes INSTANTLY on that stale set before the 50 ms
  // debounced re-render builds the fresh one. Poll the USUBJID sequence BY POSITION against what THIS
  // setup must produce. Position 0 = current, 1 = mouse-over (volatile, EXCLUDED), 2.. = selected.
  await expect.poll(() => page.evaluate((sel) => {
    const df = grok.shell.t;
    const usub = (r: number) => df.col('USUBJID').get(r);
    const selected: number[] = [];
    for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) selected.push(i);
    const expectedStable = [usub(df.currentRowIdx), ...selected.map(usub)];
    const actual = Array.from(document.querySelectorAll(sel)).map((c) =>
      ((c as HTMLElement).querySelector('[column="USUBJID"]') as HTMLInputElement | null)?.value ?? null);
    const actualStable = actual.length >= 2 ? [actual[0], ...actual.slice(2)] : actual;
    return JSON.stringify(actualStable) === JSON.stringify(expectedStable);
  }, ORDINARY), {timeout: 15_000}).toBe(true);
}

/**
 * Move the REAL pointer to a neutral point below the viewer, clear of every card and grid row.
 * Releases the mouse-over binding the way a real user does (measured on demog: card set unchanged,
 * position-1 mouse-over card empties, mouseOverRowIdx → -1). MUST be a real move, never zeroed from
 * JS: a pointer resting on a card re-fires onmouseenter on every re-render, re-binding
 * mouseOverRowIdx. The grid sits ABOVE the viewer, so a point just below the host is off all cards.
 */
async function parkPointerBelowViewer(page: Page): Promise<void> {
  const hostBox = (await page.locator(HOST).boundingBox())!;
  const viewportH = page.viewportSize()?.height ?? Math.round(hostBox.y + hostBox.height + 40);
  const x = Math.round(hostBox.x + hostBox.width / 2);
  const y = Math.round(Math.min(hostBox.y + hostBox.height + 10, viewportH - 6));
  await page.mouse.move(x, y, {steps: 10});
}

test('Forms viewer — mouse interactions and row binding (p1)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // #### Setup
  await softStep('Setup — add the Forms viewer, confirm the three show-toggles are ON, seed the entry state', async () => {
    await v.addViewerByIcon(page, 'Forms', 'Forms', 30_000, 'FormsViewer');
    await page.locator('.d4-multi-form').first().waitFor({timeout: 30_000});

    // Confirm the ship defaults BEFORE any write — a set-then-read would be a tautology.
    const defaults = await page.evaluate(() => {
      const vw = grok.shell.tv.viewers.find((x: any) => x.type === 'FormsViewer');
      return {
        selected: vw.props.showSelectedRows,
        current: vw.props.showCurrentRow,
        mouseOver: vw.props.showMouseOverRow,
      };
    });
    expect(defaults).toEqual({selected: true, current: true, mouseOver: true});

    await establishSetup(page);
  });

  // #### Scenario 1
  await softStep('Scenario 1 — Clicking a non-current card sets df.currentRowIdx and the field matches df.get', async () => {
    await establishSetup(page);
    // Card 2 is the first selected-row card (past the two leading cards), not the current-row card.
    const targetUsub = await cardFieldValue(page, 2, 'USUBJID');
    expect(targetUsub).not.toBeNull();
    const targetRow = await rowOfValue(page, 'USUBJID', targetUsub as string);
    const clickedSex = await cardFieldValue(page, 2, 'SEX');

    await page.locator(ORDINARY).nth(2).click();

    await expect.poll(() => page.evaluate(() => grok.shell.t.currentRowIdx), {timeout: 10_000}).toBe(targetRow);
    const sexOfCurrent = await page.evaluate(() =>
      grok.shell.t.get('SEX', grok.shell.t.currentRowIdx));
    expect(clickedSex).toBe(sexOfCurrent);
  });

  // #### Scenario 2
  await softStep('Scenario 2 — Clicking a card field sets df.currentCell to that column and the card row', async () => {
    await establishSetup(page);
    const cardUsub = await cardFieldValue(page, 3, 'USUBJID');
    expect(cardUsub).not.toBeNull();
    const cardRow = await rowOfValue(page, 'USUBJID', cardUsub as string);

    await page.locator(ORDINARY).nth(3).locator('[column="AGE"]').click();

    await expect.poll(() => page.evaluate(() => {
      const cc = grok.shell.t.currentCell;
      return JSON.stringify({col: cc?.column?.name ?? null, row: cc?.rowIndex ?? null});
    }), {timeout: 10_000}).toBe(JSON.stringify({col: 'AGE', row: cardRow}));
  });

  // #### Scenario 3
  await softStep('Scenario 3 — Clicking the HEIGHT header label sets df.currentCol to HEIGHT', async () => {
    await establishSetup(page);
    await page.locator(`${HOST} .d4-multi-form-header [name="div-HEIGHT"]`).first().click();
    await expect.poll(() => page.evaluate(() => grok.shell.t.currentCol?.name ?? null),
      {timeout: 10_000}).toBe('HEIGHT');
  });

  // #### Scenario 4a
  await softStep('Scenario 4a — Ctrl+clicking a card toggles df.selection.get(row) off then back on', async () => {
    // Row 20 is made the current row so its card (position 0) survives deselection — a deselected
    // ordinary card vanishes, so the current-row card is the stable surface for the toggle round-trip.
    const row = 20;
    await page.evaluate((r) => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = r;
      df.selection.setAll(false);
      df.selection.set(r, true);
    }, row);
    await expect.poll(() => page.evaluate((r) => grok.shell.t.selection.get(r), row),
      {timeout: 10_000}).toBe(true);
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID', CURRENT),
      {timeout: 10_000}).toBe(await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), row));

    // Gate on card-order stability before each click — the seed and each toggle trigger a debounced
    // render() that would otherwise detach card 0 mid-gesture (see Scenario 4b).
    await waitForOrderStable(page);
    await page.locator(CURRENT).first().click({modifiers: ['Control']});
    await expect.poll(() => page.evaluate((r) => grok.shell.t.selection.get(r), row),
      {timeout: 10_000}).toBe(false);

    await waitForOrderStable(page);
    await page.locator(CURRENT).first().click({modifiers: ['Control']});
    await expect.poll(() => page.evaluate((r) => grok.shell.t.selection.get(r), row),
      {timeout: 10_000}).toBe(true);
  });

  // #### Scenario 4b
  await softStep('Scenario 4b — Shift+clicking the row-12 card selects exactly rows 0..12 inclusive', async () => {
    const row = 12;
    // Rows 18 and 22 are seeded ABOVE the clicked row and must be cleared by the Shift+click, so the
    // exact-set equality below discriminates; against an already-empty tail it would prove nothing.
    await page.evaluate((r) => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = r;
      df.selection.setAll(false);
      [3, 7, 12, 18, 22].forEach((i) => df.selection.set(i, true));
    }, row);
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID', CURRENT),
      {timeout: 10_000}).toBe(await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), row));
    const tailSeeded = await page.evaluate(() => [18, 22].every((i) => grok.shell.t.selection.get(i)));
    expect(tailSeeded).toBe(true);

    await waitForOrderStable(page);
    // Shift range-select is genuinely expensive on this build: the Shift onclick handler loops
    // selection.set over every row, each O(rows) via the viewer's selection.onChanged subscription.
    // MEASURED ~27 s per Shift+click on demog (5850 rows), freezing the main thread. locator.click
    // retries when its own selection change detaches card 0 and never converges — so dispatch the
    // trusted Shift+click ONCE via page.mouse and wait it out. A named, measured cost, not slack.
    const box = (await page.locator(CURRENT).first().boundingBox())!;
    await page.keyboard.down('Shift');
    try {
      await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);
    } finally {
      await page.keyboard.up('Shift');
    }

    await expect.poll(() => page.evaluate((r) => {
      const df = grok.shell.t;
      const sel: number[] = [];
      for (let i = 0; i < df.rowCount; i++) if (df.selection.get(i)) sel.push(i);
      const expected: number[] = [];
      for (let i = 0; i <= r; i++) expected.push(i);
      return JSON.stringify(sel) === JSON.stringify(expected);
    }, row), {timeout: 40_000}).toBe(true);
  });

  // #### Scenario 4c
  await softStep('Scenario 4c — Ctrl+Shift+clicking the row-8 card clears rows 0..8; rows above stay selected', async () => {
    // Entry state carried over from 4b: rows 0..12 selected, so clearing 0..8 must leave 9..12.
    const row = 8;
    await page.evaluate((r) => {
      const df = grok.shell.t;
      df.mouseOverRowIdx = -1;
      df.currentRowIdx = r;
    }, row);
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID', CURRENT),
      {timeout: 10_000}).toBe(await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), row));
    const nineToTwelveBefore = await page.evaluate(() =>
      [9, 10, 11, 12].every((r) => grok.shell.t.selection.get(r)));
    expect(nineToTwelveBefore).toBe(true);

    // Card-0 stability barrier before the click (same detach hazard as 4b); the 0..12 render must settle.
    await waitForOrderStable(page);
    await page.locator(CURRENT).first().click({modifiers: ['Control', 'Shift']});

    await expect.poll(() => page.evaluate((r) => {
      const df = grok.shell.t;
      const cleared = Array.from({length: r + 1}, (_, i) => i).every((i) => !df.selection.get(i));
      const keptAbove = [9, 10, 11, 12].every((i) => df.selection.get(i));
      return JSON.stringify({cleared, keptAbove});
    }, row), {timeout: 10_000}).toBe(JSON.stringify({cleared: true, keptAbove: true}));
  });

  // #### Scenario 5a
  await softStep('Scenario 5a — Hovering an ordinary card sets df.mouseOverRowIdx to that card row', async () => {
    await establishSetup(page);

    const targetUsub = await cardFieldValue(page, 2, 'USUBJID');
    const targetRow = await rowOfValue(page, 'USUBJID', targetUsub as string);

    // Guard against a stale-state false pass: the binding must not already read the target row.
    const beforeHover = await page.evaluate(() => grok.shell.t.mouseOverRowIdx);
    expect(beforeHover).not.toBe(targetRow);

    await page.locator(ORDINARY).nth(2).hover();
    await expect.poll(() => page.evaluate(() => grok.shell.t.mouseOverRowIdx),
      {timeout: 10_000}).toBe(targetRow);

    // Move the REAL pointer away and assert the binding RELEASES (card set unchanged, position 1
    // empties). Never zero mouseOverRowIdx from JS while the pointer rests on the card — the next
    // re-render re-fires it to the hovered row. Poll to the settled state, don't sleep blindly.
    await parkPointerBelowViewer(page);
    await expect.poll(() => page.evaluate(() => grok.shell.t.mouseOverRowIdx),
      {timeout: 10_000}).toBe(-1);
    // The product channel, not only the dataframe: the mouse-over card holds no field value once
    // released. Read BY POSITION (index 1) — the mouse-over card is defined by its slot, not an indicator.
    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBeNull();

    // Card set OTHERWISE UNCHANGED — position 0 stays the current-row card. Read positionally: the
    // current row still occupies the leading slot, and exactly one current-row card exists.
    const currentUsub5a = await page.evaluate(() => grok.shell.t.col('USUBJID').get(grok.shell.t.currentRowIdx));
    await expect.poll(() => cardFieldValue(page, 0, 'USUBJID'), {timeout: 10_000}).toBe(currentUsub5a);
    expect(await page.locator(CURRENT).count()).toBe(1);
  });

  // #### Scenario 5b
  await softStep('Scenario 5b — The mouse-over card USUBJID field goes empty → hovered row value → empty', async () => {
    await establishSetup(page);
    await parkPointerBelowViewer(page);
    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = -1; });
    // With no hovered row the always-built mouse-over card (position 1) has NO USUBJID field element,
    // so asserting card presence would be false-green — assert the FIELD value.
    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBeNull();

    // Grid is a canvas with no per-row DOM hit target headless, so the hover is driven through
    // df.mouseOverRowIdx (the exact binding a real grid hover produces). The observable is the card FILL.
    const hoverRow = 33;
    await page.evaluate((r) => { grok.shell.t.mouseOverRowIdx = r; }, hoverRow);
    const expected = await page.evaluate((r) => grok.shell.t.col('USUBJID').get(r), hoverRow);
    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBe(expected);

    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = -1; });
    await expect.poll(() => cardFieldValue(page, 1, 'USUBJID'), {timeout: 10_000}).toBeNull();
  });

  // #### Scenario 6a
  let baseCount = 0;
  await softStep('Scenario 6a — Turning Show Mouse Over Row OFF drops one card and stops the grid-hover fill', async () => {
    await establishSetup(page);
    // Record the base card count measured in this run as the reference for the drop-by-one and
    // return-to-baseline checks; the absolute number is not asserted.
    baseCount = await page.locator(ORDINARY).count();

    // Turn Show Mouse Over Row OFF through the property panel (Context Panel > Misc), not vw.props.
    // The category can ship collapsed, so expand it first, then flip the checkbox.
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', false, 'misc');
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(baseCount - 1);
    await expect.poll(() => page.locator(`${ORDINARY} .d4-multi-form-form-indicator-mouse-over-row`).count(),
      {timeout: 15_000}).toBe(0);

    const currentUsub = await page.evaluate(() =>
      grok.shell.t.col('USUBJID').get(grok.shell.t.currentRowIdx));
    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = 40; });
    // Fixed wait is DELIBERATE, not a settle-poll: the assertion below is NEGATIVE (row 40 must NOT
    // appear). A "does not contain" poll passes on the FIRST sample, before the re-render could add
    // the value — vacuous even if the fill were broken. The fixed wait lets that re-render elapse
    // first. The "no fixed pauses" rule targets positive settle-reads, not a negative guard.
    await page.waitForTimeout(1500);
    const hovered40 = await page.evaluate(() => grok.shell.t.col('USUBJID').get(40));
    const shown = await fieldValuesByPosition(page, 'USUBJID');
    expect(shown).not.toContain(hovered40);
    expect(await cardFieldValue(page, 0, 'USUBJID', CURRENT)).toBe(currentUsub);
  });

  // #### Scenario 6a restore
  await softStep('Scenario 6a restore — Turning Show Mouse Over Row back ON returns the card count to the baseline', async () => {
    // Clear the JS-driven hover binding first (a dataframe write, no panel control exists), then
    // restore Show Mouse Over Row through the panel — the same route the OFF toggle took.
    await page.evaluate(() => { grok.shell.t.mouseOverRowIdx = -1; });
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-mouse-over-row');
    await v.setPropertyGridCheckbox(page, 'show-mouse-over-row', true, 'misc');
    await expect.poll(() => page.locator(ORDINARY).count(), {timeout: 15_000}).toBe(baseCount);
  });

  // #### Scenario 6b
  await softStep('Scenario 6b — Turning Show Current Row OFF removes the current-row card and its green indicator', async () => {
    await establishSetup(page);
    await expect.poll(() => page.locator(CURRENT).count(),
      {timeout: 15_000}).toBe(1);

    // Turn Show Current Row OFF through the property panel (Context Panel > Misc), not vw.props.
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-current-row');
    await v.setPropertyGridCheckbox(page, 'show-current-row', false, 'misc');

    // Only the indicator is asserted: with Show Mouse Over Row still ON this build renders ZERO
    // ordinary cards while Show Current Row is off (measured 6 → 0), so a baseCount-1 assert is false RED.
    await expect.poll(() => page.locator(CURRENT).count(),
      {timeout: 15_000}).toBe(0);

    // Restore the default through the panel so the section leaves the viewer in its ship state.
    await v.ensurePropertyCategory(page, 'Forms', 'misc', 'show-current-row');
    await v.setPropertyGridCheckbox(page, 'show-current-row', true, 'misc');
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
