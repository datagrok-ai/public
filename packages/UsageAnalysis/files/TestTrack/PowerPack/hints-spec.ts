/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: ui-smoke
//   sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
//   ui_coverage_responsibility: [add-new-column-function-hint-tooltip] (delegated_to: add-new-column.md)
//   related_bugs: []
//   produced_from: migrated
//
// Bug-library cross-reference:
//   Consulted bug-library/powerpack.yaml — no curated bug intersects this
//   scenario's sub_features_covered. GROK-17109 (calculated column persistence
//   across save+datasync+reopen) and GROK-17004 (paste handler crash on complex
//   formulas) touch the same sub-features but their reproduction surfaces do
//   not exercise the hover/tooltip flow; they are covered by sibling specs
//   (Powerpack/add-new-column-spec.ts, Powerpack/formula-refreshing-spec.ts,
//   Powerpack/highlight-spec.ts).
//
// Delegation note:
//   The basic dialog UI surface (open/close, preview grid, OK/Cancel,
//   Recent Activities autofill) is owned by Powerpack/add-new-column-spec.ts.
//   This spec owns the specialty add-new-column-function-hint-tooltip flow.
//   Dataset open + dialog open + function insertion are PRECONDITIONS for the
//   hover assertion — they are not in this spec's ui_coverage_responsibility.
//   The hover step (Step 4) MUST be DOM-driven (pyramid_layer ui-smoke);
//   substituting CodeMirror dispatch or direct hoverTooltip invocation would
//   bypass the owned UI flow.
//
// Reference template: Powerpack/add-new-column-spec.ts (same directory, same
//   ui-smoke pyramid_layer, same powerpack feature). Same login + setup
//   pattern; same dialog-open mechanic via [name="icon-add-new-column"];
//   same CodeMirror surface (.add-new-column-dialog-cm-div .cm-content);
//   same autocomplete tooltip selector (.cm-tooltip-autocomplete) for Step 3.
//
// Source citations for selectors:
//   - Toolbar icon: [name="icon-add-new-column"] — verified in
//     grok-browser/references/dialogs-menus.md "Toolbar Ribbon Icons" table
//     (precedent: Powerpack/add-new-column-spec.ts, autocomplete-spec.ts).
//   - Dialog root: .d4-dialog filtered by hasText 'Add New Column' (precedent:
//     same sibling specs).
//   - CodeMirror surface: .cm-content inside .add-new-column-dialog-cm-div —
//     PowerPack/src/dialogs/add-new-column.ts:143/175 attaches that class to
//     the codeMirrorDiv host; inner editor exposes .cm-content (CodeMirror 6
//     default).
//   - Autocomplete tooltip: .cm-tooltip-autocomplete — rendered by the
//     @codemirror/autocomplete package (precedent: autocomplete-spec.ts).
//   - Function-signature hover tooltip: .cm-tooltip.cm-tooltip-hover — the
//     CodeMirror 6 hoverTooltip extension's standard render class. The dialog
//     wires the extension via _hoverTooltip / hoverTooltipCustom at
//     PowerPack/src/dialogs/add-new-column.ts:142/424/1004-1021. The tooltip
//     body is a <div> with textContent = res.signature (e.g. "Abs(num)").
//   - Cancel button: [name="button-Add-New-Column---CANCEL"] — set by
//     prepareForSeleniumTests in
//     PowerPack/src/dialogs/add-new-column.ts:344-349 (called from init() at
//     line 250; unconditional once selenium-mode is on).
//
// Atlas provenance (derived_from):
//   feature-atlas/powerpack.yaml#sub_features[powerpack.dialogs.add-new-column].interactions
//     derived_from: public/packages/PowerPack/src/tests/add-new-column.ts#L68
//     (the apitests-layer 'hints' test that asserts the in-dialog hintDiv
//     signature surface; THIS spec exercises the parallel hover-tooltip
//     surface, not the hintDiv, per the scenario body's explicit "Hover the
//     mouse pointer over the inserted function name" wording in Step 4.)
//   feature-atlas/powerpack.yaml#sub_features[powerpack.dialogs.add-new-column-func].interactions
//     derived_from: public/packages/PowerPack/src/package.ts#L405
//     (toolbar icon + Edit | Add New Column top-menu both open AddNewColumnDialog).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('PowerPack: Add new column — hover inserted function name surfaces signature tooltip (demog)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  // ---- Login + setup phase ----
  await loginToDatagrok(page);

  // ---- Step 1: open the demog dataset (precondition; delegated to parent) ----
  // The scenario's Setup paragraph is satisfied by the parent add-new-column.md
  // (per ui_coverage_delegated_to). The dataset is still required as a
  // precondition for the dialog to open; using readCsv via JS API for the
  // setup phase is consistent with sibling Powerpack specs.
  await softStep('Step 1: open demog dataset; verify grid renders', async () => {
    await page.evaluate(async () => {
      const grok = (window as any).grok;
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      try { grok.shell.closeAll(); } catch (_) { /* best effort */ }
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    await page.waitForTimeout(1000);
    // Sanity: demog has the standard numerical columns the dialog will list.
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? df.columns.names() : [];
    });
    expect(cols.length).toBeGreaterThan(0);
  });

  // ---- Step 2: open the Add New Column dialog via toolbar icon ----
  // Precondition for the hover assertion; the dialog-open mechanic itself is
  // owned by the parent smoke (add-new-column-spec.ts Step 1).
  await softStep('Step 2: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    // Verify the formula editor (CodeMirror) is mounted alongside the dialog.
    const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
    await cm.waitFor({timeout: 15_000, state: 'visible'});
  });

  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();

  // ---- Step 3: insert a function into the formula editor via type + autocomplete ----
  // Scenario body offers three insertion paths (type + autocomplete, click
  // entry in functions panel, drag function onto field). Type + autocomplete
  // is the canonical example in the scenario ("type `a` and select `Abs` from
  // the autocomplete tooltip; the function is inserted as `Abs(num)`"); it is
  // also the path with the most stable DOM (sibling autocomplete-spec.ts
  // exercises the same mechanic). Insertion is a precondition for hover —
  // not in this spec's ui_coverage_responsibility — so any supported path is
  // acceptable per the scenario's "using any supported approach" wording.
  await softStep('Step 3: type "a", accept Abs from autocomplete tooltip; verify "Abs(num)" inserted', async () => {
    await cm.click();
    await page.waitForTimeout(150);
    // Clear any pre-existing content defensively.
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    // Type "a" to trigger the autocomplete tooltip.
    await page.keyboard.type('a', {delay: 60});
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    // Accept the first match — the autocomplete list is alphabetically
    // anchored on the typed prefix; "Abs" is the canonical first match
    // for "a" per the scenario example. Try clicking the "Abs" entry
    // explicitly; if textContent matching fails (registry drift), fall
    // back to Enter on the highlighted entry.
    const clickedAbs = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return false;
      const items = Array.from(tip.querySelectorAll('li'));
      const abs = items.find((li) => /\bAbs\b/.test(li.textContent || ''));
      if (!abs) return false;
      abs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      abs.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      (abs as HTMLElement).click();
      return true;
    });
    if (!clickedAbs)
      await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
    // Verify the editor doc now contains "Abs(" — the function name plus
    // open-paren parameter list. The exact param-placeholder text depends on
    // the function's signature (e.g. "num"); we assert structurally.
    const doc = await page.evaluate(() => {
      const cmDiv = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmDiv) return '';
      const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
      return view ? view.state.doc.toString() : (cmDiv.innerText || '');
    });
    const normalized = doc.replace(/\s+/g, ' ').trim();
    expect(/^Abs\([^)]*\)/.test(normalized)).toBe(true);
    // Dismiss any lingering autocomplete tooltip so the hover tooltip in
    // Step 4 has an unobstructed render surface.
    await page.keyboard.press('Escape');
    await page.waitForTimeout(150);
  });

  // ---- Step 4: OWNED FLOW — hover the inserted function name; verify signature tooltip ----
  // ui_coverage_responsibility flow: add-new-column-function-hint-tooltip.
  // pyramid_layer ui-smoke → MUST be DOM-driven. We use page.hover() on the
  // CodeMirror text range covering the "Abs" function name. The CodeMirror 6
  // hoverTooltip extension (wired at PowerPack/src/dialogs/add-new-column.ts:
  // 142/424/1004-1021) responds to a real mouse hover on the function name
  // by rendering .cm-tooltip.cm-tooltip-hover whose <div> textContent is the
  // function's signature ("Abs(num)").
  await softStep('Step 4: hover over inserted function name; verify signature tooltip surfaces', async () => {
    // Locate the bounding box of the "Abs" text inside the CodeMirror
    // contenteditable surface. CM6 renders text inside .cm-line spans, so
    // we walk the DOM to find the text node carrying "Abs" and compute the
    // pixel rect at its mid-point. Playwright's page.mouse.move({steps: ...})
    // then drives a real mouse-move sequence over that point — this is the
    // input the hoverTooltip extension listens for.
    const rect = await page.evaluate(() => {
      const cmDiv = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmDiv) return null;
      // Walk text nodes inside CM6 .cm-line containers and find the one
      // that contains the substring "Abs". The first text node is the
      // inserted function name.
      const walker = document.createTreeWalker(cmDiv, NodeFilter.SHOW_TEXT);
      let node: Node | null;
      while ((node = walker.nextNode())) {
        const t = (node as Text).data;
        const idx = t.indexOf('Abs');
        if (idx >= 0) {
          const range = document.createRange();
          range.setStart(node, idx);
          range.setEnd(node, idx + 3); // length of "Abs"
          const r = range.getBoundingClientRect();
          if (r.width > 0 && r.height > 0)
            return {x: r.left + r.width / 2, y: r.top + r.height / 2};
        }
      }
      return null;
    });
    expect(rect).not.toBeNull();
    // Move the mouse away first (off the editor) so the subsequent move into
    // the function name is a fresh hover gesture that the extension picks up.
    await page.mouse.move(10, 10);
    await page.waitForTimeout(120);
    // Drive a real mouse-move sequence into the function-name pixel rect.
    // Using `steps` makes the move a continuous gesture rather than a
    // teleport — CM6's hoverTooltip listens to mousemove events with a
    // brief debounce, so a continuous gesture is more reliable.
    await page.mouse.move(rect!.x, rect!.y, {steps: 10});
    // The hoverTooltip extension has an internal hover delay (default ~300ms
    // in CM6 view package). Wait up to 5s for the .cm-tooltip-hover element
    // to appear and become visible.
    const hover = page.locator('.cm-tooltip.cm-tooltip-hover').first();
    await hover.waitFor({timeout: 5_000, state: 'visible'});
    await expect(hover).toBeVisible();
    // Verify the tooltip body carries the function's signature.
    // hoverTooltipCustom (add-new-column.ts:1014-1018) sets
    //   const dom = document.createElement('div');
    //   dom.textContent = res.signature;
    // The signature for Abs is "Abs(num)" — function name plus its input
    // parameter list. We assert the tooltip text contains the function name
    // ("Abs") AND a parenthesised parameter list — the scenario's
    // verification wording is "function name plus its input parameter list,
    // e.g. Abs(num)". A param-name-agnostic structural match avoids coupling
    // to a specific registry snapshot.
    const tooltipText = (await hover.textContent() ?? '').trim();
    expect(tooltipText.length).toBeGreaterThan(0);
    // Function name present.
    expect(/\bAbs\b/.test(tooltipText)).toBe(true);
    // Parenthesised parameter list present (e.g. "(num)", "(num: number)",
    // etc.). The assertion deliberately allows for signature-format drift
    // while still proving the tooltip is the signature surface, not, say,
    // an empty popup or an error tooltip.
    expect(/Abs\s*\([^)]*\)/.test(tooltipText)).toBe(true);
  });

  // ---- Cleanup: dismiss the dialog (CANCEL) and close the table view ----
  // The scenario does not commit the new column; CANCEL leaves df untouched.
  await page.evaluate(() => {
    const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => { /* best effort */ });
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best effort */ }
  }).catch(() => { /* best effort */ });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
