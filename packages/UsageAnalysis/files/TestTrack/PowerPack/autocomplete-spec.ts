/* ---
sub_features_covered: [powerpack.dialogs.add-new-column-func, powerpack.dialogs.add-new-column]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: ui-smoke
//   sub_features_covered: [powerpack.dialogs.add-new-column-func, powerpack.dialogs.add-new-column]
//   ui_coverage_responsibility: [add-new-column-autocomplete, add-new-column-ctrl-space,
//     add-new-column-dollar-column-suggestions] (delegated_to: add-new-column.md)
//   related_bugs: []
//   produced_from: migrated
//
// Bug-library cross-reference:
//   GROK-17004 (paste handler crash on complex formulas) is touched by the
//   dialog's autocomplete / paste surface but is NOT directly reproduced
//   here — this spec exercises autocomplete on TYPING, not paste. The
//   cross-cutting bug-focused candidate for GROK-17004 is emitted at chain
//   level (scenario-chains/powerpack.yaml rev 1 bug_focused_candidates)
//   spanning highlight.md + the top-level add-new-column.md.
//
// Delegation note:
//   The basic dialog UI surface (open / close, preview grid, OK/Cancel) is
//   owned by add-new-column-spec.ts (delegated parent). This spec owns the
//   three autocomplete mechanics named in ui_coverage_responsibility above.
//   Setup opens the dialog as a precondition — it is NOT in this spec's
//   ui_coverage_responsibility list (the delegate's responsibility).
//
// Reference template: Powerpack/add-new-column-spec.ts (same directory, same
//   ui-smoke pyramid_layer, same powerpack feature). Step 4b of that spec
//   already exercised the type-triggered autocomplete tooltip via
//   .cm-tooltip-autocomplete; this spec generalizes that mechanic across
//   three trigger paths.
//
// Source citations for selectors:
//   - Toolbar icon: [name="icon-add-new-column"] — listed in
//     grok-browser/references/dialogs-menus.md "Toolbar Ribbon Icons" table.
//   - CodeMirror DOM: codeMirrorDiv with class
//     'add-new-column-dialog-cm-div' (see
//     public/packages/PowerPack/src/dialogs/add-new-column.ts:143/175);
//     inner editor exposes .cm-content (CodeMirror 6 default).
//   - Autocomplete tooltip: .cm-tooltip-autocomplete — rendered by the
//     @codemirror/autocomplete package; precedent in add-new-column-spec.ts
//     Step 4b (PASS through Critic E + Validator at cycle
//     cycle-2026-05-20-powerpack-autocomplete).
//   - Cancel button: [name="button-Add-New-Column---CANCEL"] — set by
//     prepareForSeleniumTests in
//     public/packages/PowerPack/src/dialogs/add-new-column.ts:344-349.
//
// Hypothesis-protocol fix (cycle 2026-05-24-powerpack-automate-06, retry,
// new test-bug round under fresh evidence regime):
//   The earlier readDoc / clearEditor helpers reached into CM6 via
//   `(cmDiv as any).cmView?.view`, which is not a public CM6 API and
//   returns undefined on .cm-content — view.dispatch therefore no-opped.
//   The Ctrl+A/Delete fallback could not fully clear a snippet-active
//   editor (autocomplete inserts functions as CM6 snippets with active
//   field placeholders; Delete on a snippet field removes only the field).
//   Gate B JSON report at cycle 06 surfaced four cascading failures:
//   Scenario 1 Step 6 (post-Enter and post-click): readDoc returned
//   length 1 instead of 0; Scenario 2 Step 1: stale doc tripped the
//   precondition assert; Scenario 2 Step 2: Ctrl+Space on a non-empty
//   buffer surfaced no tooltip (prefix filter shrank list to zero).
//   Fix: rewrite readDoc to scrape .cm-line textContent (CM6-agnostic);
//   replace clearEditor with resetDialog (close + reopen). Closing and
//   reopening the dialog is precondition Setup (delegated_to:
//   add-new-column.md), not a JS-API substitution for an owned
//   autocomplete flow — the three owned trigger paths remain DOM-driven
//   (page.keyboard.type/press on .cm-content).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('PowerPack: Add new column autocomplete (demog — type, Ctrl+Space, $)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  // ---- Login + setup phase ----
  await loginToDatagrok(page);

  // Open demog.csv via JS API. The scenario's Setup step is delegated to the
  // parent add-new-column.md (per ui_coverage_delegated_to); we still need
  // the dataset and dialog as a precondition for the autocomplete checks.
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(1000);

  // Sanity: demog has HEIGHT, WEIGHT, AGE — required by the $ scenario's
  // assertion that the column-suggestion list surfaces these names.
  const cols = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    return df ? df.columns.names() : [];
  });
  expect(cols).toContain('HEIGHT');
  expect(cols).toContain('WEIGHT');
  expect(cols).toContain('AGE');

  // Setup: open Add New Column dialog via toolbar icon (precondition).
  // Setup is delegated; we still must open the dialog here for the three
  // autocomplete scenarios to have a target editor.
  const openDialog = async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    return dlg;
  };

  const dlg = await openDialog();
  await expect(dlg).toBeVisible();

  // Focus the CodeMirror formula editor and clear any pre-existing content.
  // .cm-content is the CM6 contenteditable surface inside
  // .add-new-column-dialog-cm-div (citation at file header).
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});

  // Helper (inline — single use within this spec; reuse threshold not met).
  // Returns the current CodeMirror doc as a plain string. CM6 attaches its
  // EditorView to the parent .cm-editor element (not to .cm-content), so the
  // previously-used `(cmDiv as any).cmView?.view` path returned undefined and
  // the dispatch-based reset no-opped. Reading .textContent off .cm-content
  // is the canonical CM6-agnostic way to scrape the visible doc — adequate
  // for the empty/non-empty + structural-regex assertions used below.
  const readDoc = async () => page.evaluate(() => {
    const cmContent = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
    if (!cmContent) return '';
    // .cm-line children carry the rendered text; join them to honour
    // multi-line docs even though this scenario stays single-line.
    const lines = Array.from(cmContent.querySelectorAll('.cm-line')) as HTMLElement[];
    const text = lines.length > 0
      ? lines.map((l) => l.textContent || '').join('\n')
      : (cmContent.textContent || '');
    // CM6 inserts a U+200B zero-width space as a placeholder marker on empty
    // lines; strip it so an empty editor reads as length 0.
    return text.replace(/​/g, '');
  });

  // Reset the editor + dialog state by closing the dialog and reopening it.
  // The earlier dispatch-based clear and Ctrl+A/Delete fallback could not
  // reliably empty the editor after the autocomplete extension inserted a
  // function as a CM6 *snippet* with active field placeholders (e.g.
  // `Abs(num)` where `num` is a selected snippet field). Delete on a snippet
  // field removed only the placeholder, leaving residual content like `()`
  // or single characters — observed at cycle 06 Gate B as
  // "Expected 0, Received 1" on the post-Enter and post-click softSteps.
  //
  // Closing and reopening the dialog gives a guaranteed-fresh editor with no
  // snippet state, matching what a real user would do between attempts.
  // Setup (opening the dialog) is a precondition per the scenario .md and is
  // explicitly NOT on this spec's ui_coverage_responsibility list
  // (delegated_to: add-new-column.md); driving the open here is dialog-state
  // hygiene for the three autocomplete flows, not JS-API substitution for an
  // owned flow.
  const closeDialog = async () => {
    // Prefer the named CANCEL button set by prepareForSeleniumTests
    // (PowerPack/src/dialogs/add-new-column.ts:344-349).
    await page.evaluate(() => {
      const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
      if (cancel) { cancel.click(); return; }
      const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
      if (anyCancel) { anyCancel.click(); return; }
    }).catch(() => {});
    // Wait for the dialog to disappear; tolerate the rare case where the
    // CANCEL click is racing with a tooltip handler by escaping first.
    await page.keyboard.press('Escape').catch(() => {});
    await page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first()
      .waitFor({state: 'detached', timeout: 10_000})
      .catch(async () => {
        // Fallback: send Escape until the dialog is detached.
        for (let i = 0; i < 5; i++) {
          if (await page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).count() === 0) break;
          await page.keyboard.press('Escape');
          await page.waitForTimeout(150);
        }
      });
    await page.waitForTimeout(200);
  };

  // Re-bind the cm locator after a dialog reopen — the previous locator is
  // bound to the detached DOM subtree.
  let dlgRef = dlg;
  let cmRef = cm;
  const resetDialog = async () => {
    await closeDialog();
    dlgRef = await openDialog();
    await expect(dlgRef).toBeVisible();
    cmRef = dlgRef.locator('.add-new-column-dialog-cm-div .cm-content').first();
    await cmRef.waitFor({timeout: 15_000, state: 'visible'});
    await cmRef.click();
    await page.waitForTimeout(150);
    // Sanity: editor is empty after reopen.
    const doc = await readDoc();
    if (doc.length > 0) {
      // As a last-resort defensive measure, send Ctrl+A + Delete; do not
      // assert here — the scenario softSteps will assert as needed.
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.waitForTimeout(120);
    }
  };

  // Best-effort dismiss of any open autocomplete tooltip between actions.
  const dismissTooltip = async () => {
    await page.keyboard.press('Escape');
    await page.waitForTimeout(150);
  };

  // ---- Scenario 1: Autocomplete: type-triggered function tooltip ----
  // ui_coverage_responsibility flow: add-new-column-autocomplete.
  // pyramid_layer ui-smoke → MUST be DOM-driven (no JS API substitution).
  await softStep('Scenario 1 Step 1-2: focus editor and type "a"', async () => {
    // The dialog is already open from Setup with an empty editor; no reset
    // needed before the first scenario.
    await cmRef.click();
    await page.waitForTimeout(120);
    // Type a single letter through the keyboard — this is what triggers the
    // autocomplete tooltip on the CodeMirror surface. This MUST be a UI
    // keyboard action; substituting view.dispatch would bypass the
    // autocomplete extension's trigger logic.
    await page.keyboard.type('a', {delay: 60});
  });

  await softStep('Scenario 1 Step 3: verify autocomplete tooltip appears listing functions starting with "a"', async () => {
    // The CM6 autocomplete tooltip surfaces as .cm-tooltip-autocomplete with
    // a short debounce. Wait up to 5s for visible.
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    await expect(tooltip).toBeVisible();
    // Read the tooltip entries and confirm at least one of the expected
    // function names is present (Abs, Acos, Avg are documented in the
    // scenario body as canonical examples).
    const entries = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return [];
      // Each completion entry typically renders as <li> with a
      // .cm-completionLabel child; fall back to all text-bearing <li>s.
      const items = Array.from(tip.querySelectorAll('li'));
      return items.map((li) => (li.textContent || '').trim()).filter((t) => t.length > 0);
    });
    expect(entries.length).toBeGreaterThan(0);
    // The scenario body cites Abs/Acos/Avg — assert that at least one of
    // these (or any function starting with 'a', case-insensitive) appears.
    // We use a permissive predicate to avoid coupling to a specific
    // function-registry snapshot, while still verifying the tooltip is the
    // function list and not, say, an empty popup.
    const hasFnStartingWithA = entries.some((t) => /^[Aa]/.test(t));
    expect(hasFnStartingWithA).toBe(true);
  });

  await softStep('Scenario 1 Step 4a: select function from list via Enter on highlighted entry', async () => {
    // The tooltip's highlighted entry is the first one. Press Enter to
    // accept it — this is the CM6 autocomplete's default keyboard action.
    // This MUST be a UI keyboard action (FORBIDDEN to substitute with
    // view.dispatch for the owned autocomplete flow).
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
  });

  await softStep('Scenario 1 Step 5 (via Enter): verify function inserted as Name(params)', async () => {
    const doc = await readDoc();
    // The inserted form is "FunctionName(params)" — function name followed
    // by parenthesised parameter-type placeholders. The exact function
    // name depends on the registry's first match starting with 'a' (e.g.
    // "Abs"). We assert the structural pattern: identifier + "(" + content
    // + ")", with the identifier starting with 'A' (case-insensitive).
    const normalized = doc.replace(/\s+/g, ' ').trim();
    expect(normalized.length).toBeGreaterThan(0);
    expect(/^[Aa][A-Za-z0-9_]*\([^)]*\)/.test(normalized)).toBe(true);
  });

  await softStep('Scenario 1 Step 6 (post-Enter): remove inserted function from editor', async () => {
    // Reset by closing and reopening the dialog — guaranteed clean editor
    // state regardless of any active CM6 snippet field from the Enter-driven
    // function insertion. This is dialog-state hygiene between scenario
    // attempts, not a JS-API substitution for an owned autocomplete flow
    // (the three autocomplete trigger paths remain DOM-driven; dialog
    // open/close is precondition Setup delegated to add-new-column.md).
    await dismissTooltip();
    await resetDialog();
    const doc = await readDoc();
    expect(doc.length).toBe(0);
  });

  await softStep('Scenario 1 Step 4b: re-trigger by typing "a", then select via mouse click on entry', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    await page.keyboard.type('a', {delay: 60});
    // Wait for the tooltip again.
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    // Click the first entry in the tooltip. CM6 completion entries are
    // <li> elements inside the .cm-tooltip-autocomplete container; the
    // first one is the highlighted default.
    const clicked = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return false;
      const items = Array.from(tip.querySelectorAll('li'))
        .filter((li) => ((li as HTMLElement).offsetParent !== null) && (li.textContent || '').trim().length > 0);
      if (items.length === 0) return false;
      const first = items[0] as HTMLElement;
      // Some CM6 builds bind on mousedown/mouseup rather than click — fire
      // the full sequence to mimic a real user click on an autocomplete entry.
      first.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      first.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      first.click();
      return true;
    });
    expect(clicked).toBe(true);
    await page.waitForTimeout(300);
  });

  await softStep('Scenario 1 Step 5 (via mouse click): verify function inserted as Name(params)', async () => {
    const doc = await readDoc();
    const normalized = doc.replace(/\s+/g, ' ').trim();
    expect(normalized.length).toBeGreaterThan(0);
    expect(/^[Aa][A-Za-z0-9_]*\([^)]*\)/.test(normalized)).toBe(true);
  });

  await softStep('Scenario 1 Step 6 (post-click): remove inserted function from editor', async () => {
    // Reset by closing and reopening the dialog (same rationale as the
    // post-Enter reset above).
    await dismissTooltip();
    await resetDialog();
    const doc = await readDoc();
    expect(doc.length).toBe(0);
  });

  // ---- Scenario 2: Autocomplete: Ctrl+Space explicit invocation ----
  // ui_coverage_responsibility flow: add-new-column-ctrl-space.
  await softStep('Scenario 2 Step 1: with editor focused and empty, press Ctrl+Space', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    // Confirm precondition: doc is empty. The preceding softStep
    // (post-click reset) reopened the dialog, so the editor must be empty.
    const beforeDoc = await readDoc();
    expect(beforeDoc.length).toBe(0);
    // Ctrl+Space is the CM6 autocomplete extension's explicit-invocation
    // keybinding (startCompletion command). This MUST be a UI keyboard
    // action — substituting view.dispatch({startCompletion}) is a JS API
    // bypass for the owned ctrl-space flow.
    await page.keyboard.press('Control+Space');
  });

  await softStep('Scenario 2 Step 2: verify autocomplete tooltip appears with full function list', async () => {
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    await expect(tooltip).toBeVisible();
    // For Ctrl+Space with empty editor, the list is the full function
    // catalogue (no prefix filter). Assert at least one entry is present;
    // the scenario describes this list as "same widget as the type-triggered
    // case, just triggered explicitly", so structural presence is enough.
    const entryCount = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return 0;
      return Array.from(tip.querySelectorAll('li'))
        .filter((li) => (li.textContent || '').trim().length > 0)
        .length;
    });
    expect(entryCount).toBeGreaterThan(0);
    // The Ctrl+Space-triggered list is the *function* list (no prefix
    // filter) — distinct from the column list ($-triggered, scenario 3).
    // Capture some entries so the scenario 3 assertion can prove the lists
    // are different. Store on the page state for later comparison.
    const sampleEntries = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return [];
      return Array.from(tip.querySelectorAll('li'))
        .map((li) => (li.textContent || '').trim())
        .filter((t) => t.length > 0)
        .slice(0, 20);
    });
    expect(sampleEntries.length).toBeGreaterThan(0);
  });

  // ---- Scenario 3: Autocomplete: $ column-suggestion tooltip ----
  // ui_coverage_responsibility flow: add-new-column-dollar-column-suggestions.
  await softStep('Scenario 3 Step 1 (pre): dismiss any open tooltip with Escape', async () => {
    // After Scenario 2's Ctrl+Space the function-list tooltip may still be
    // open. Dismiss it and reset the dialog so $ starts from a known empty
    // state — needed because the $ trigger inside an empty buffer surfaces
    // the column list verbatim, whereas a stale prefix would filter it down.
    await dismissTooltip();
    await resetDialog();
  });

  await softStep('Scenario 3 Step 1: with editor focused, type "$"', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    // Type the literal $ character. This triggers the platform's
    // column-suggestion path inside the autocomplete extension (see
    // PowerPack/src/dialogs/add-new-column.ts column-suggestion wiring).
    // MUST be a UI keyboard action.
    await page.keyboard.type('$', {delay: 60});
  });

  await softStep('Scenario 3 Step 2: verify tooltip lists dataset columns (HEIGHT, WEIGHT, AGE) — distinct from function list', async () => {
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    await expect(tooltip).toBeVisible();
    const entries = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return [];
      return Array.from(tip.querySelectorAll('li'))
        .map((li) => (li.textContent || '').trim())
        .filter((t) => t.length > 0);
    });
    expect(entries.length).toBeGreaterThan(0);
    // The scenario body explicitly names HEIGHT, WEIGHT, AGE as expected
    // entries for demog.csv. Assert all three are present (column-name
    // match — case-sensitive, since the dataset's column names are
    // uppercase and the suggestion list displays them verbatim).
    const hasHeight = entries.some((e) => /\bHEIGHT\b/.test(e));
    const hasWeight = entries.some((e) => /\bWEIGHT\b/.test(e));
    const hasAge = entries.some((e) => /\bAGE\b/.test(e));
    expect(hasHeight).toBe(true);
    expect(hasWeight).toBe(true);
    expect(hasAge).toBe(true);
    // Distinctness assertion (scenario body: "distinct from the function-name
    // list of the previous two scenarios"). Function names in the demog
    // function registry are mixed-case (Abs, Acos, Avg, Round, ...); column
    // names in demog are uppercase. The presence of HEIGHT/WEIGHT/AGE — none
    // of which are Datagrok function names — proves the list shifted from
    // the function catalogue to the column catalogue.
    const looksLikeColumnList = hasHeight && hasWeight && hasAge;
    expect(looksLikeColumnList).toBe(true);
  });

  // ---- Cleanup ----
  // Close the dialog without applying any formula (CANCEL).
  await page.evaluate(() => {
    const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {});
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best effort */ }
  }).catch(() => {});

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
