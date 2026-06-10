/* ---
sub_features_covered: [models.validators.class-imbalance, models.validators.string-features, models.validators.too-many-unique-categories, models.validators.highly-correlated]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (breadth-extension edge coverage; no pyramid claim)
//   sub_features_covered: [models.validators.class-imbalance,
//                          models.validators.string-features,
//                          models.validators.too-many-unique-categories,
//                          models.validators.highly-correlated]
//   ui_coverage_responsibility: [] (no owned UI flow; section ui-smoke stays
//                                   with predictive-models.md)
//   related_bugs: []
//   coverage_type: edge
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/models.yaml#sub_features[models.validators.class-imbalance] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L68
//   feature-atlas/models.yaml#sub_features[models.validators.string-features] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L134
//   feature-atlas/models.yaml#sub_features[models.validators.too-many-unique-categories] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L31
//   feature-atlas/models.yaml#sub_features[models.validators.highly-correlated] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L44
//
// Spec target: edge — four pre-train validator surfaces in the
// PredictiveModelingView. Each scenario builds a small in-memory dataframe
// that deterministically triggers one validator predicate; the spec asserts
// the corresponding tip text surfaces in the "Insights & Tips" widget and
// that no error-balloon / no validation block was raised (`isError = false`
// is the validator-source contract for all four predicates per
// predictive_modeling_validators.dart#L37, L54, L106, L143).
//
// Paradigm: playwright DOM-driving for the train-UI plumbing (top-menu,
// Predict picker via column-selector mousedown, Features picker via canvas
// overlay clicks) so the empirical end-to-end validator-pipeline trigger is
// exercised (parameters form → runDataValidators → modelPreview.addTip).
// JS API used for in-memory dataset construction (DG.DataFrame.fromColumns)
// + warnings tray inspection — non-owned setup paradigm per
// `pyramid_layer: absent` constraint defaults.

// Cycle-02 hypothesis-retry recon-notes (2026-06-10, cycle
// 2026-06-09-models-automate-02):
//   Carried-forward Gate B failure from cycle 01: 3 attempts × ~123s,
//   failure_keys [B-RUN-PASS, B-STAB-01]. Cycle 01's round-2 spec used
//   trusted Playwright `page.mouse.click` at sweep coords led by
//   (dx=40, dy=28) on the Select-columns dialog overlay canvas. That fix
//   was correct in *coordinates* but appears unreliable in *mechanism*:
//   the trusted-click path empirically stabilises only some attempts
//   under cold-start, matching the B-STAB-01 (instability across attempts)
//   flake signature.
//
//   This-dispatch (cycle-02 retry-1) live MCP recon on dev.datagrok.ai
//   2026-06-10 (PredictiveModelingView opened against the same 200-row
//   3-column class-imbalance dataframe, Predict = target set, Features
//   editor opened) confirmed empirically that **synthetic dispatchEvent
//   on the overlay canvas at (box.right - 40, box.y + 28) toggles the
//   row 0 checkbox deterministically** (1→0→1→0 across 4 consecutive
//   single-clicks; the sweep across dx ∈ {15..50}, dy ∈ {20..40} found
//   the entire checkbox-column hit zone responds — see live recon log).
//   The synthetic-events path is *not* refuted on this surface (the
//   round-1 "synthetic refuted" claim was about the chemprop overlay, a
//   different widget). This makes synthetic the more deterministic
//   primary path for the row-0 checkbox toggle here.
//
//   Same-paradigm tactical fix (NOT a paradigm pivot — both synthetic
//   dispatchEvent and Playwright page.mouse.click are canvas-click
//   mechanisms on the same overlay; switching the *preference* between
//   the two is tactical within the canvas-click paradigm): make synthetic
//   dispatchEvent the primary toggle path at (dx=40, dy=28), then fall
//   back to page.mouse.click sweep if synthetic somehow misses. Also
//   tighten the post-`label-None` wait so the sweep starts from a known
//   "0 checked" baseline (avoid the readChecked race that the prior
//   dispatch may have hit cold-start).
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser
// reference) — covers all four scenarios:
//   .d4-pm-model-widget — top-level wrapper div emitted by
//     PredictiveModelPreviewWidget.makeWidget (predictive_modeling_view.dart#L38).
//     Each widget hosts a header (.d4-dialog-title) + body. Multiple instances
//     coexist inside the PMV root (one per preview section: "Insights & Tips",
//     "Performance", "Features", "Predicted X vs Actual", "Residuals", ...).
//     Observed live 2026-06-09 via chrome-devtools MCP on dev.datagrok.ai
//     build 1.28.0 — root.querySelectorAll('.d4-pm-model-widget') returned 8
//     widgets on a trained PMV and 1 widget (the tips host) on a Predict+
//     Features-set untrained PMV. NOT in grok-browser/references/models.md —
//     that file L487-499 incorrectly states "models.validators.* — NOT
//     surfaced as a discrete UI element". The actual surface is this widget,
//     selected by `.d4-pm-model-widget` whose `.d4-dialog-title` text equals
//     "Insights & Tips" — that scoping is load-bearing because the unscoped
//     `.d4-pm-model-widget` selector matches multiple widgets.
//   .d4-pm-model-widget .d4-dialog-title — header span inside each widget
//     carrying the widget's title (e.g. "Insights & Tips", "Performance",
//     "Features", "Residuals"). Observed live 2026-06-09 — querySelector on
//     each widget returned the title span verbatim, allowing widget filtering
//     by text content. The "Insights & Tips" widget is the validator-tips
//     host per PredictiveModelPreviewWidget.clear() at
//     predictive_modeling_view.dart#L24.
//   .d4-pm-model-widget li (under the "Insights & Tips" host) — each tip
//     renders as a top-level <li> inside the widget's <ul>; nested <li>s
//     carry sub-bullets (e.g. class-imbalance: outer <li> = column name,
//     inner <li>s = imbalanced categories with ratios). Markup is produced by
//     ValidatorHelper.messageToHtml + Markup.render at
//     predictive_modeling_view.dart#L46-49. Asserting on innerText of the
//     widget host catches the validator phrase regardless of nested-list
//     structure. Observed live 2026-06-09: class-imbalance tip rendered as
//     "Some columns contain class imbalance:\n\ntarget\nA (1.800)\nB (0.200)"
//     in widget.innerText on the class-imbalance test dataframe; string-
//     features tip rendered as "Column 'cat_feature' is categorical. Most
//     models require converting them to numerical." on the string-features
//     test dataframe.
//   [name="dialog-Select-columns..."] [name="viewer-Grid"] [name="overlay"]
//     — canvas overlay for the column-picker (Select-columns... dialog) inside
//     the PMV Features picker. Empirically verified live 2026-06-10 via
//     chrome-devtools MCP on dev.datagrok.ai: the row-0 checkbox cell hit-zone
//     spans dx ∈ {15..55} from box.right and dy ∈ {20..40} from box.y for the
//     3-column dataframes this spec builds; synthetic dispatchEvent('mousedown'
//     / 'mouseup' / 'click') with bubbles:true and clientX/clientY in the hit
//     zone toggles "N checked" deterministically (4 consecutive single-clicks
//     observed: 0→1→0→1→0). The empirical center is (dx=40, dy=28). NOT in
//     grok-browser/references/models.md.

import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// ────────────────────────────── helpers ──────────────────────────────────

// Open ML > Models > Train Model... and wait for the PredictiveModelingView
// to mount. The Dart top-menu ignores Playwright .hover() — synthesize the
// mouseover/mouseenter/mousemove sequence directly on the ML > Models entry
// so the submenu reveals "Train Model..." (pattern from train-spec.ts).
async function openTrainView(page: Page) {
  await page.locator('[name="div-ML"]').click();
  await page.evaluate(() => {
    const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
    if (!models) throw new Error('ML > Models submenu not found');
    const r = models.getBoundingClientRect();
    const ev = (type: string) => new MouseEvent(type, {
      bubbles: true, cancelable: true, view: window,
      clientX: r.left + 5, clientY: r.top + 5,
    });
    models.dispatchEvent(ev('mouseover'));
    models.dispatchEvent(ev('mouseenter'));
    models.dispatchEvent(ev('mousemove'));
  });
  await page.locator('[name="div-ML---Models---Train-Model..."]').click();
  await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel',
    null, {timeout: 15_000});
  // Wait for the PMV form to be fully reconciled — Predict editor mounted, the
  // "Insights & Tips" widget host present. Without this wait, subsequent
  // setPredict() raced the form mount on cold-start (retry attempts under
  // --repeat-each surface the race as B-STAB-01: scenario 1 PASSes warm,
  // scenario 1 of attempt-2/3 FAILs because the PMV root has not finished
  // mounting before setPredict opens the column-selector).
  await page.waitForFunction(() => {
    const root = (window as any).grok.shell.v?.root;
    if (!root) return false;
    const predictHost = root.querySelector('[name="input-host-Predict"] .d4-column-selector');
    const tipsWidget = Array.from(root.querySelectorAll('.d4-pm-model-widget'))
      .find((w: Element) => w.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Insights & Tips');
    return !!predictHost && !!tipsWidget;
  }, null, {timeout: 15_000});
  await page.waitForTimeout(500);
}

// Set Predict via the .d4-column-selector mousedown popup. Same pattern as
// train-spec.ts:116-133 + chemprop-spec.ts:62-76; programmatic
// `view.predict = ...` does not visibly update the form per the
// grok-browser/references/models.md gotcha at L202-208.
async function setPredict(page: Page, columnName: string) {
  // Skip if already on the target column (idempotent — useful when scenarios
  // happen to land on the desired column after PMV open).
  const current = await page.evaluate(() =>
    (window as any).grok.shell.v.root
      .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim());
  if (current === columnName) return;
  await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const sel = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
    sel.dispatchEvent(new MouseEvent('mousedown', {
      bubbles: true, cancelable: true, view: window, button: 0,
    }));
  });
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout: 5_000});
  await page.evaluate(() =>
    (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
  await page.keyboard.type(columnName);
  await page.waitForTimeout(400);
  await page.keyboard.press('Enter').catch(() => {});
  await page.waitForTimeout(700);
  // If the backdrop is still up (Enter race), close it via Escape and retry
  // by clicking the suggestion list entry directly.
  const stillOpen = await page.evaluate(() => !!document.querySelector('.d4-column-selector-backdrop'));
  if (stillOpen) {
    await page.keyboard.press('Escape').catch(() => {});
    await page.waitForTimeout(300);
  }
  await expect.poll(async () => await page.evaluate(() =>
    (window as any).grok.shell.v.root
      .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
    {timeout: 10_000}).toBe(columnName);
}

// Open the Features Select-columns... dialog and check each named column
// individually via the proven chemprop-spec pattern (search filter + trusted
// Playwright page.mouse.click on the canvas overlay). Round-2 hypothesis-
// protocol MCP recon 2026-06-09 (B-RUN-PASS+B-STAB-01 retry): previous
// dispatch used "All" + OK which selects EVERY column including target —
// predict-conflict validator at predictive_modeling_view.dart#L351-355
// REJECTS (sets d4-invalid, training blocked) but does NOT auto-remove
// target → Features render as "(3) All" (synthetic summary, NOT column-name
// list) with class d4-invalid → assertion ".ui-input-column-names contains
// 'feature1'" fails deterministically across all 3 attempts (B-STAB-01
// is a false-positive flake flag — failure was deterministic same-text
// "(3) All"). Same-paradigm tactical fix: per sibling chemprop-spec.ts
// L209-251 (which already pivoted from synthetic dispatchEvent to trusted
// page.mouse.click after empirically refuting synthetic events on this same
// canvas), use label-None + search + trusted-click + OK. Empirical backing:
// chemprop-spec.ts (already-landed) drives the IDENTICAL Select-columns dialog
// with this pattern; MCP recon 2026-06-09 confirmed [name="label-None"],
// [name="label-All"], and input[placeholder="Search"] are all addressable
// on a 3-column dataframe (the same fixtures this spec builds).
async function selectFeatures(page: Page, columnNames: string[]) {
  await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const editor = root.querySelector('[name="div-Features"]') as HTMLElement;
    const r = editor.getBoundingClientRect();
    const opts = (t: string) => new MouseEvent(t, {
      bubbles: true, cancelable: true, view: window, button: 0,
      clientX: r.left + 10, clientY: r.top + 5,
    });
    editor.dispatchEvent(opts('mousedown'));
    editor.dispatchEvent(opts('mouseup'));
    editor.dispatchEvent(opts('click'));
  });
  const dlg = page.locator('[name="dialog-Select-columns..."]');
  await dlg.waitFor({timeout: 10_000});
  // Clear via label-None (label is a real DOM element — synthetic .click() works).
  await dlg.locator('[name="label-None"]').click();
  // Helper: read "N checked" count from any label inside the dialog.
  const readChecked = async () => {
    return await dlg.evaluate((d) => {
      const lbl = Array.from(d.querySelectorAll('label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      const m = lbl?.textContent?.match(/^(\d+)/);
      return m ? parseInt(m[1], 10) : 0;
    });
  };
  // Wait for "0 checked" baseline after label-None — without this, the
  // first row sweep starts before the click was committed and `before`
  // may read a stale non-zero value, breaking the +1 increment detection.
  await expect.poll(readChecked, {timeout: 5_000, intervals: [100, 200, 300],
    message: 'Select-columns dialog never reached 0 checked after label-None'}).toBe(0);
  // For each desired column: filter the canvas via the search input (so the
  // target column is the only row visible), then toggle the row-0 checkbox
  // cell via synthetic dispatchEvent at the empirical hit-zone (dx=40,
  // dy=28). Live MCP recon 2026-06-10 verified synthetic dispatchEvent
  // toggles the checkbox deterministically on this overlay (different
  // surface from the chemprop overlay which refused synthetic events).
  // Fall back to trusted page.mouse.click sweep if synthetic somehow
  // misses (defensive — kept as a safety net for cold-start layout
  // perturbation).
  const search = dlg.locator('input[placeholder="Search"]');
  for (const name of columnNames) {
    const before = await readChecked();
    await search.fill('');
    await page.waitForTimeout(200);
    await search.fill(name);
    await page.waitForTimeout(500);
    const overlay = dlg.locator('[name="viewer-Grid"] [name="overlay"]');
    const box = await overlay.boundingBox();
    if (!box) throw new Error('column-picker overlay canvas not measurable');
    let toggled = false;
    // Primary path: synthetic dispatchEvent at the recon-verified hit
    // zone (dx=40, dy=28). MCP recon 2026-06-10 single-click toggles
    // deterministic across 4 consecutive presses (1→0→1→0).
    for (const [dx, dy] of [[40, 28], [40, 34], [42, 28], [35, 28]]) {
      const clientX = box.x + box.width - dx;
      const clientY = box.y + dy;
      const synthetic = await dlg.evaluate((d, {cx, cy}) => {
        const ov = d.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement;
        if (!ov) return 'no-overlay';
        const opts = (t: string) => new MouseEvent(t, {
          bubbles: true, cancelable: true, view: window,
          button: 0, clientX: cx, clientY: cy,
        });
        ov.dispatchEvent(opts('mousedown'));
        ov.dispatchEvent(opts('mouseup'));
        ov.dispatchEvent(opts('click'));
        return 'dispatched';
      }, {cx: clientX, cy: clientY});
      if (synthetic !== 'dispatched') break;
      await page.waitForTimeout(300);
      if ((await readChecked()) === before + 1) {
        toggled = true;
        break;
      }
    }
    // Fallback: trusted page.mouse.click sweep across the broader
    // hit-zone if synthetic missed (cold-start layout shift safety net).
    if (!toggled) {
      outer: for (const dx of [40, 35, 45, 30, 50]) {
        for (const dy of [28, 32, 24, 36, 20, 40]) {
          await page.mouse.click(box.x + box.width - dx, box.y + dy);
          await page.waitForTimeout(300);
          if ((await readChecked()) === before + 1) {
            toggled = true;
            break outer;
          }
        }
      }
    }
    if (!toggled)
      throw new Error(`Failed to toggle "${name}" via synthetic dispatchEvent + page.mouse.click fallback on overlay (search-filtered to "${name}")`);
  }
  // Restore the unfiltered view before OK so the dialog commits the actual
  // checked set (defensive — checked state should be preserved through
  // filter changes, but a stale filter narrows what OK pushes back).
  await search.fill('');
  await page.waitForTimeout(200);
  const finalChecked = await readChecked();
  if (finalChecked !== columnNames.length) {
    throw new Error(
      `Expected ${columnNames.length} checked after toggles, got ${finalChecked}`);
  }
  await dlg.locator('[name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('[name="dialog-Select-columns..."]'),
    null, {timeout: 10_000});
  // After OK closes, the Features input should NOT carry d4-invalid (no
  // target overlap) and should render the column-name list rather than the
  // synthetic "(N) All" summary. Assert on absence of d4-invalid — this is
  // the load-bearing check that the predict-conflict validator did not
  // reject the selection.
  await expect.poll(async () => await page.evaluate(() => {
    const inp = document.querySelector('[name="input-host-Features"] .ui-input-column-names');
    return inp?.className.includes('d4-invalid') ? 'invalid' : 'valid';
  }), {timeout: 5_000, message: 'Features input should not be d4-invalid after selection'})
    .toBe('valid');
  // Wait for the action-checkbox host to materialize — confirms the form has
  // reconciled with Predict + Features and the validator pipeline has run.
  await page.waitForFunction(() => {
    const host = document.querySelector('[name="input-host-Impute-missing"]') as HTMLElement | null;
    return !!host && host.offsetParent !== null;
  }, null, {timeout: 10_000}).catch(() => {
    // Action checkboxes may not surface for all engine/target combinations.
    // The validator widget content is the real assertion — proceed without
    // hard-failing here.
  });
  // Brief settle for runDataValidators (debounced) to populate tips.
  await page.waitForTimeout(1_000);
}

// Wait until the "Insights & Tips" widget's text contains the expected
// validator phrase. Returns the widget's full innerText for assertions on
// secondary content (column names, ratios, correlation pairs).
async function getTipsText(page: Page): Promise<string> {
  return await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const widgets = Array.from(root.querySelectorAll('.d4-pm-model-widget')) as HTMLElement[];
    const tips = widgets.find((w) =>
      w.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Insights & Tips');
    return tips ? (tips as HTMLElement).innerText : '';
  });
}

async function waitForTipsContaining(page: Page, phrase: string, timeout = 30_000) {
  await expect.poll(async () => (await getTipsText(page)).toLowerCase(),
    {timeout, intervals: [500, 1000, 1500], message:
      `waiting for "Insights & Tips" widget to contain "${phrase}"`})
    .toContain(phrase.toLowerCase());
}

// Verify the validator is non-blocking: no error-state balloon present, no
// platform error in shell.warnings. Per predictive_modeling_validators.dart
// L37 (too-many-unique-categories), L54 (highly-correlated), L106 (class-
// imbalance), L143 (string-features) — all four predicates set
// `isError = false` so the train UI remains enterable. The SAVE button stays
// d4-disabled (no successful train yet — that's expected pre-train state per
// models.md L259-260), so we do NOT assert SAVE-enabled; we assert (a) the
// SAVE button element exists in DOM (the form is mounted and the validator
// path did not crash it), (b) no error-balloon, (c) no error warnings.
async function assertNonBlocking(page: Page) {
  await expect(page.locator('[name="button-Save"]')).toBeAttached({timeout: 5_000});
  const errorBalloons = await page.locator('.d4-balloon-error').count();
  expect(errorBalloons).toBe(0);
  // grok.shell.warnings accumulates across the session; restrict to fatal
  // platform errors (objects with isError/level=='error'), not any string
  // containing "error" (which catches benign messages like "no errors").
  const errCount = await page.evaluate(() => {
    const w: any[] = (window as any).grok.shell.warnings ?? [];
    return w.filter((x) =>
      x && typeof x === 'object' && (x.isError === true || x.level === 'error')).length;
  });
  expect(errCount).toBe(0);
}

// Build an in-memory dataframe via JS API + open it in a TableView. The
// validator pipeline runs against the dataframe's columns directly, so
// dataset construction shape (types, lengths, content) is the actual
// test input.
async function openInMemoryDataFrame(page: Page, builder: () => unknown) {
  await page.evaluate(async (builderSrc: string) => {
    const g: any = (window as any).grok;
    g.shell.closeAll();
    // Wait for view tear-down (PMV from a prior scenario must close before we
    // mount a new TableView, otherwise the next openTrainView races).
    await new Promise<void>((resolve) => {
      const start = Date.now();
      const poll = () => {
        if (!g.shell.v || g.shell.v.type === undefined || Date.now() - start > 3000) resolve();
        else setTimeout(poll, 100);
      };
      poll();
    });
    await new Promise((r) => setTimeout(r, 400));
    // eslint-disable-next-line @typescript-eslint/no-implied-eval
    const df = (new Function('DG', `return (${builderSrc})(DG)`))((window as any).DG);
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, builder.toString());
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  // Confirm the TableView is the active view before openTrainView races top-menu click.
  await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'TableView',
    null, {timeout: 10_000});
}

// ────────────────────────────── test ─────────────────────────────────────

test('Models validators — pre-train edge surfaces (class-imbalance, string-features, too-many-unique, highly-correlated)', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
  });

  // ─────────── Scenario 1: class-imbalance (categorical target ±20%) ────────

  await softStep('1.1 Build 200-row class-imbalance dataframe (target = 180×A / 20×B)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 200;
      const f1 = new Float32Array(n);
      const f2 = new Float32Array(n);
      const tgt: string[] = [];
      // Deterministic seeds for reproducibility (no Math.random()).
      for (let i = 0; i < n; i++) {
        f1[i] = ((i * 9301 + 49297) % 233280) / 233280;
        f2[i] = ((i * 1597 + 12345) % 65536) / 65536;
        // 180 of 200 = 'A' (count ratio = 1.8), 20 of 200 = 'B' (ratio = 0.2)
        // → both outside 1±0.2 → class-imbalance fires per validator source L83.
        tgt.push(i < 180 ? 'A' : 'B');
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array('feature1', f1),
        DG.Column.fromFloat32Array('feature2', f2),
        DG.Column.fromStrings('target', tgt),
      ]);
    });
  });

  await softStep('1.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('1.3 Set Predict = target', async () => {
    await setPredict(page, 'target');
  });

  await softStep('1.4 Set Features = (feature1, feature2) — exclude target to avoid predict-conflict d4-invalid', async () => {
    await selectFeatures(page, ['feature1', 'feature2']);
    // Per the predict-conflict validator at predictive_modeling_view.dart#L351-355,
    // Features ∩ Predict must be empty for the input to be valid. selectFeatures
    // already asserted absence of d4-invalid after OK; additionally confirm the
    // text rendering carries the chosen column names (not the "(N) All" synthetic
    // summary that prior dispatch surfaced when target was incorrectly included).
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feature1', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feature2');
  });

  await softStep('1.5 Insights & Tips surfaces class-imbalance warning citing target', async () => {
    // sub_feature: models.validators.class-imbalance
    // Validator source (predictive_modeling_validators.dart#L101) prepends
    // "Some columns contain class imbalance:" then lists the offending column
    // with its out-of-range categories + ratios.
    await waitForTipsContaining(page, 'class imbalance');
    const tipsText = await getTipsText(page);
    expect(tipsText).toContain('target');
    // Validator caps at 4 imbalanced columns per source L107; with only one
    // categorical target column, exactly one column is listed.
    expect(tipsText).toMatch(/A \(\d+\.\d+\)|B \(\d+\.\d+\)/);
  });

  await softStep('1.6 Validator is non-blocking (no error-balloon, no error warnings, SAVE button mounted)', async () => {
    await assertNonBlocking(page);
  });

  // ─────────── Scenario 2: string-features (StringColumn feature) ───────────

  await softStep('2.1 Build 60-row string-features dataframe (cat_feature = red/green/blue)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 60;
      const num = new Float32Array(n);
      const tgt = new Float32Array(n);
      const cats: string[] = [];
      for (let i = 0; i < n; i++) {
        num[i] = ((i * 9301 + 49297) % 233280) / 233280;
        tgt[i] = (((i * 1597 + 12345) % 65536) / 65536) * 10;
        cats.push(['red', 'green', 'blue'][i % 3]);
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromStrings('cat_feature', cats),
        DG.Column.fromFloat32Array('num_feature', num),
        DG.Column.fromFloat32Array('target', tgt),
      ]);
    });
  });

  await softStep('2.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('2.3 Set Predict = target (numerical regression target)', async () => {
    await setPredict(page, 'target');
  });

  await softStep('2.4 Set Features = (cat_feature, num_feature) — exclude target', async () => {
    await selectFeatures(page, ['cat_feature', 'num_feature']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('cat_feature', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('num_feature');
  });

  await softStep('2.5 Insights & Tips surfaces string-features hint citing cat_feature', async () => {
    // sub_feature: models.validators.string-features
    // Validator source (predictive_modeling_validators.dart#L137-142) emits
    // "Column 'X' is categorical. Most models require converting them to
    // numerical." for each StringColumn feature with no semType. The phrase
    // "convert" + "numerical" is the actionable hint paired with the
    // one-hot preprocessing action per the atlas binding.
    await waitForTipsContaining(page, 'categorical');
    const tipsText = await getTipsText(page);
    expect(tipsText).toContain('cat_feature');
    expect(tipsText.toLowerCase()).toContain('convert');
    expect(tipsText.toLowerCase()).toContain('numerical');
  });

  await softStep('2.6 Validator is non-blocking', async () => {
    await assertNonBlocking(page);
  });

  // ─────────── Scenario 3: too-many-unique-categories (>0.8) ────────────────

  await softStep('3.1 Build 50-row dataframe with id_like StringColumn (categories/length = 1.0)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 50;
      const ids: string[] = [];
      const f1 = new Float32Array(n);
      const tgt = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        // Each row unique → categories.length / column.length = 1.0 → above
        // the 0.8 threshold (predictive_modeling_validators.dart#L31).
        ids.push(`row_${i}`);
        f1[i] = ((i * 9301 + 49297) % 233280) / 233280;
        tgt[i] = (((i * 1597 + 12345) % 65536) / 65536) * 10;
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromStrings('id_like', ids),
        DG.Column.fromFloat32Array('feature1', f1),
        DG.Column.fromFloat32Array('target', tgt),
      ]);
    });
  });

  await softStep('3.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('3.3 Set Predict = target', async () => {
    await setPredict(page, 'target');
  });

  await softStep('3.4 Set Features = (id_like, feature1) — exclude target', async () => {
    await selectFeatures(page, ['id_like', 'feature1']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('id_like', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feature1');
  });

  await softStep('3.5 Insights & Tips surfaces too-many-unique-categories warning citing id_like', async () => {
    // sub_feature: models.validators.too-many-unique-categories
    // Validator source (predictive_modeling_validators.dart#L37) emits
    // "Column(s) 'X' contain(s) too many unique categories." for any
    // categorical column where categories/length > 0.8. Note: stringFeatures
    // ALSO fires on id_like (it is also a StringColumn with no semType) —
    // this is correct per source: the two predicates are independent and may
    // both list the same column. We assert specifically on the
    // "too many unique categories" phrase to scope to this validator.
    await waitForTipsContaining(page, 'too many unique categories');
    const tipsText = await getTipsText(page);
    expect(tipsText).toContain('id_like');
  });

  await softStep('3.6 Validator is non-blocking', async () => {
    await assertNonBlocking(page);
  });

  // ─────────── Scenario 4: highly-correlated (Pearson > 0.9) ────────────────

  await softStep('4.1 Build 30-row dataframe with feat_a + feat_b (Pearson > 0.9)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 30;
      const a = new Float32Array(n);
      const b = new Float32Array(n);
      const tgt = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        const base = i;
        a[i] = base;
        // Small deterministic noise — < 5% of the signal range so Pearson
        // correlation remains > 0.9 well above the validator threshold
        // (predictive_modeling_validators.dart#L48).
        b[i] = base + (((i * 7919) % 100) / 100) * 0.5;
        tgt[i] = base * 0.3 + (((i * 1597 + 12345) % 17) / 17) * 0.1;
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array('feat_a', a),
        DG.Column.fromFloat32Array('feat_b', b),
        DG.Column.fromFloat32Array('target', tgt),
      ]);
    });
  });

  await softStep('4.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('4.3 Set Predict = target', async () => {
    await setPredict(page, 'target');
  });

  await softStep('4.4 Set Features = (feat_a, feat_b) — exclude target', async () => {
    await selectFeatures(page, ['feat_a', 'feat_b']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feat_a', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feat_b');
  });

  await softStep('4.5 Insights & Tips surfaces highly-correlated warning citing feat_a <> feat_b pair', async () => {
    // sub_feature: models.validators.highly-correlated
    // Validator source (predictive_modeling_validators.dart#L56-60) emits
    // "Columns are highly correlated:\n* A <> B" for each pair with
    // Pearson correlation > 0.9. Validator scans only when 2..100 numerical
    // columns are present (source L46) — our 3 numerical columns satisfy
    // that. The pair text uses " <> " (space-angle-space) as the delimiter.
    await waitForTipsContaining(page, 'highly correlated');
    const tipsText = await getTipsText(page);
    // Pair may be listed as "feat_a <> feat_b" or "feat_b <> feat_a" depending
    // on column iteration order; accept either form.
    expect(tipsText).toMatch(/feat_a\s*<>\s*feat_b|feat_b\s*<>\s*feat_a/);
  });

  await softStep('4.6 Validator is non-blocking', async () => {
    await assertNonBlocking(page);
  });

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
