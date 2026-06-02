/* ---
sub_features_covered:
  - bio.search.subsequence
  - bio.search.subsequence.top-menu
  - bio.search.subsequence.editor
  - bio.search.subsequence.filter
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: integration
//   sub_features_covered: [bio.search.subsequence, .subsequence.top-menu,
//     .subsequence.editor, .subsequence.filter]
//   ui_coverage_responsibility: ["Bio | Search | Subsequence Search",
//     subsequence-search-filter, filter-panel-sequence-input,
//     filter-panel-reset-filter] (delegated_to: null)
//   related_bugs: [] (scenario carries no related_bugs — GROK-16111
//     covers bio.search.similarity surface, not subsequence per chain
//     bug_match_attempts_skipped no_affecting_scenarios entry)
//   produced_from: migrated
//   coverage_type: regression
//
// Atlas provenance: realizes critical path bio.cp.subsequence-search
// (priority p0) end-to-end — top-menu dispatch -> filter widget docking
// -> filter-panel input -> grid-filter callback -> Reset Filter
// restores full row count. Multi-subsystem integration (atlas
// pyramid_layer: integration; ui-smoke slot is held by manage.md per
// Rule 1 cardinality).
//
// Selector sources:
//   Class 1 — bio.md (grok-browser reference) documented:
//     [name="div-Bio"] (bio.md L606)
//     [name="div-Bio---Search---Subsequence-Search-..."] (bio.md L246/L608,
//       trailing "-" before "..." REQUIRED)
//     [name="viewer-Filters"] (bio.md L252)
//     [name="viewer-Filters"] [data-source="Bio:Bio Substructure Filter"]
//       (bio.md L253/L618)
//     input[placeholder="Substructure"] inside the filter widget
//       (bio.md L256/L618)
//     [name="div-View---Reset-Filter"] (bio.md L260/L582/L621 — top-menu
//       Reset; documented but NOT used here, see Selector recon-notes)
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in bio.md):
//   [name="viewer-Filters"] .d4-filter-group-header [name="icon-arrow-rotate-left"]
//     — Filters-panel filter-group header reset icon (the rotate-left
//     "Reset filter" affordance local to the filter panel itself, not
//     the top-menu View | Reset Filter). Reached via Setup → Bio | Search
//     | Subsequence Search dispatch → filter docks into Filters panel
//     → filter-group-header is the row above the docked filter widgets.
//     Observed live 2026-06-02 via chrome-devtools MCP take_snapshot +
//     evaluate_script on dev.datagrok.ai/filter_FASTA.csv: the
//     filter-group-header carries five icons, the fourth being
//     `<i class="grok-icon fal fa-arrow-rotate-left"
//          name="icon-arrow-rotate-left" aria-label="Reset filter">`
//     — i.e. the same element carries BOTH the `[name=icon-...]` selector
//     AND the FontAwesome `.fa-arrow-rotate-left` class. The `[name=]`
//     form is preferred (project name=-first selector convention). The
//     bio.md "Reset path" note at L269 documents only the top-menu
//     `[name="div-View---Reset-Filter"]` path; the scenario step 4
//     wording "Click Reset Filter" reads as the filter-panel-local
//     reset, which this `[name="icon-arrow-rotate-left"]` icon is.
//     Empirical observations (2026-06-02): clicking the filter-group-
//     header reset icon (a) does NOT raise a confirmation dialog (no
//     .d4-dialog appears), (b) drops df.filter.trueCount back to
//     df.rowCount synchronously, and (c) tears down + rebuilds the
//     filter widget; the substructure input is briefly absent then
//     re-mounts with .value === '' within ~2 s. Earlier dispatch
//     emitted a defensive OK-button click branch; this dispatch removes
//     it because no dialog is raised on this surface.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Bio Subsequence Search filter and reset on filter_FASTA', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Setup phase: open filter_FASTA.csv (14 rows x 1 col "fasta" with
  // semType=Macromolecule, units=fasta, aligned=SEQ per bio.md "Test
  // datasets"). Wait for Macromolecule detector + Bio package init
  // (cell renderer + filter widget registration).
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv('System:AppData/Bio/tests/filter_FASTA.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 4000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
    const hasMacromolecule = cols.some((c: any) => c.semType === 'Macromolecule');
    if (hasMacromolecule) {
      for (let i = 0; i < 60; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Bio top-menu readiness poll + Bio init-completion probe. Mirrors
  // the analyze-spec.ts cold-start stabilization: the Bio top-menu entry
  // appears only once Bio package functions register against the
  // Macromolecule TableView; Bio:getSeqHelper resolves only after
  // initBio completes (bio.md "Init-order invariant").
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
  await page.evaluate(async () => {
    const probes = ['Bio:getSeqHelper', 'Bio:getMonomerLibHelper', 'Bio:getBioLib'];
    for (const fn of probes) {
      try { await (grok as any).functions.call(fn, {}); return; } catch { /* try next */ }
    }
    await new Promise((r) => setTimeout(r, 3000));
  });

  // Per-leaf function-registration probe — guarantees Bio:bioSubstructureFilter
  // / Bio:subsequenceSearchTopMenu (per atlas bio.search.subsequence.top-menu /
  // bio.search.subsequence.filter) is findable before the menu dispatch click.
  await page.evaluate(async () => {
    const names = [
      'Bio:bioSubstructureFilter', 'Bio:bioSubstructureFilterPanel',
      'Bio:subsequenceSearchTopMenu', 'Bio:subsequenceSearch',
    ];
    const findAny = (): boolean => {
      for (const n of names) {
        try {
          if ((grok as any).functions.find && (grok as any).functions.find(n)) return true;
        } catch { /* try next */ }
      }
      return false;
    };
    const deadline = Date.now() + 15_000;
    while (Date.now() < deadline) {
      if (findAny()) return;
      await new Promise((r) => setTimeout(r, 300));
    }
    // Defensive settle if no candidate is findable (function rename across
    // Bio versions); per-step 30s click tolerance below is the ceiling.
    await new Promise((r) => setTimeout(r, 1500));
  });

  // Scenario step "Setup" — verifies the open dataset shape matches the
  // bio.md "Test datasets" contract (14 rows, single Macromolecule
  // column). Scenario does not assert a specific row count explicitly,
  // but Step 5 "Verify all rows should be present" anchors the post-
  // reset assertion to the open-time count — capture and reuse it.
  let baseRows = -1;
  await softStep('Open filter_FASTA.csv (14 rows, Macromolecule fasta column)', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i));
      const macro = cols.find((c: any) => c.semType === 'Macromolecule');
      return {
        rows: df.rowCount,
        macroName: macro ? (macro as any).name : null,
        units: macro ? (macro as any).getTag('units') : null,
      };
    });
    baseRows = info.rows;
    expect(info.rows).toBeGreaterThan(0);
    expect(info.macroName).not.toBeNull();
    expect(info.units).toBe('fasta');
  });

  // Scenario step 1 — On menu ribbon, open Bio > Search > Subsequence
  // Search. Per bio.md (L246) clicking this leaf docks the Bio
  // substructure filter widget into the Filters panel; there is NO
  // modal dialog. Trailing "-" before "..." is REQUIRED in the
  // name-attribute (source caption "Subsequence Search " carries a
  // trailing space).
  await softStep('Open Bio > Search > Subsequence Search', async () => {
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      document.querySelector('[name="div-Bio---Search"]')!
        .dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Bio---Search---Subsequence-Search-..."]') as HTMLElement).click();
    });
    // Wait for the Bio substructure filter widget to dock into the
    // Filters panel (atlas bio.search.subsequence.filter realises
    // bioSubstructureFilter widget — bio.md "Apply subsequence flow"
    // step 1: docks filter widget, no dialog).
    await page.locator('[name="viewer-Filters"] input[placeholder="Substructure"]')
      .waitFor({timeout: 30_000});
    const info = await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      return {
        filterCount: fg.filters.length,
        column: (fg.filters[0] as any).columnName ?? null,
      };
    });
    expect(info.filterCount).toBe(1);
    expect(info.column).not.toBeNull();
  });

  // Scenario step 2 — On the filter panel, set Sequence to
  // "RTDEVSNHTHDKPTLTWFEEIFEEYHSP". Scenario step 3 expects "1 row"
  // (single match in filter_FASTA.csv). Native setter + input/change
  // events drive the Dart-side filter compute; the .grok-loader
  // returns to display:none when the compute completes (bio.md
  // "Apply subsequence flow" step 4).
  const subsequence = 'RTDEVSNHTHDKPTLTWFEEIFEEYHSP';
  await softStep(`Set Sequence filter to ${subsequence} -> 1 row`, async () => {
    await page.evaluate((seq) => {
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement;
      const nativeSetter = Object.getOwnPropertyDescriptor(
        window.HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(input, seq);
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      input.dispatchEvent(new KeyboardEvent('keyup', {key: 'Enter', code: 'Enter', bubbles: true}));
    }, subsequence);
    // Filter compute is async. Poll until trueCount converges (no
    // simple "loader hidden" probe — bio.md notes the .grok-loader
    // toggles transiently). Cap at 30s for cold filter widget init.
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      return df.filter.trueCount < df.rowCount;
    }, null, {timeout: 30_000});
    const state = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement;
      return {value: input.value, filtered: df.filter.trueCount, total: df.rowCount};
    });
    expect(state.value).toBe(subsequence);
    // Scenario step 3: "Verify there is 1 row."
    expect(state.filtered).toBe(1);
    expect(state.total).toBe(baseRows);
  });

  // Scenario step 4 — Click Reset Filter. The filter-group-header
  // rotate-left icon ([name="icon-arrow-rotate-left"] — live-MCP-observed
  // 2026-06-02, see Selector recon-notes above) resets all filters in
  // the panel. Empirical: NO confirmation dialog is raised on this
  // surface; df.filter.trueCount returns to df.rowCount synchronously;
  // the filter widget is briefly torn down + rebuilt and the
  // substructure input.value settles back to '' within ~2 s.
  await softStep('Click Reset Filter -> all rows present', async () => {
    await page.locator('[name="viewer-Filters"] .d4-filter-group-header [name="icon-arrow-rotate-left"]')
      .click();
    // Wait for the BitSet to reset first (synchronous on click).
    await page.waitForFunction(() => {
      const df = grok.shell.tv.dataFrame;
      return df.filter.trueCount === df.rowCount;
    }, null, {timeout: 15_000});
    // Then wait for the widget DOM to re-mount with an empty input.
    // The reset path tears down + rebuilds the filter widget; querying
    // input.value too early returns null (input not in DOM yet).
    await page.waitForFunction(() => {
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement | null;
      return input !== null && input.value === '';
    }, null, {timeout: 15_000});
    const state = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const input = document.querySelector(
        '[name="viewer-Filters"] input[placeholder="Substructure"]',
      ) as HTMLInputElement | null;
      return {value: input?.value ?? null, filtered: df.filter.trueCount, total: df.rowCount};
    });
    // Scenario step 5: "Verify all rows should be present."
    expect(state.filtered).toBe(state.total);
    expect(state.total).toBe(baseRows);
    // Reset clears the substructure input — empirical 2026-06-02 MCP
    // observation; the filter widget re-mounts with .value === ''
    // within ~2 s of the reset click.
    expect(state.value).toBe('');
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
