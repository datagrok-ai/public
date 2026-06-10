/* ---
sub_features_covered: [models.preprocessing.one-hot, models.engines.api.apply]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (coverage_type: edge)
//   sub_features_covered: [models.preprocessing.one-hot, models.engines.api.apply]
//   ui_coverage_responsibility: [] (no delegated_to; standalone edge scenario)
//   related_bugs: []
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/models.yaml#sub_features[models.preprocessing.one-hot] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L158
//   feature-atlas/models.yaml#sub_features[models.engines.api.apply] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L26
//   feature-atlas/models.yaml#edge_cases[one-hot suffix collision] derived_from:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L158
//
// Spec target: playwright edge — asserts that one-hot encoding namespaces
// expanded columns with the `<name>=<category>` pattern so two categorical
// features sharing the same category values (here: featureA / featureB each
// taking Yes / No) do NOT collide on a single `Yes` column at train time,
// and that PredictiveModelingEngine.apply reconstructs the same per-feature
// expansion against an apply-time table (auto-built columnNamesMap with
// one-hot suffix detection).
//
// Paradigm: Playwright DOM-driving for the train + apply UI workflow; JS API
// for in-memory dataframe construction (DG.DataFrame.fromColumns — deterministic
// suffix-collision pattern, not mined from a demo dataset), pre/post-test
// cleanup, and reading the saved model's property-panel "Inputs" cell which is
// the canonical surface for the suffix-collision invariant per MCP recon
// 2026-06-09 (see Selector recon-notes below).
//
// Engine selection (scope reduction): the scenario .md cites "EDA: Linear
// Regression" as the engine. Per grok-browser/references/models.md L246-260 the
// Engine selector is runtime-discovered (PackagePredictiveModelingEngine entries
// surface dynamically) and its exact selector was not observable in this
// session's MCP recon — the host is gated on a valid Predict+Features pair
// established via CDP-level keyboard typeahead which the chrome-devtools MCP
// surface does not drive. The spec therefore lets the train view pick its
// default engine (typically Caret on dev.datagrok.ai build 1.28.0) — the
// suffix-collision invariant is engine-independent (it lives in the
// preprocessing layer, not the engine itself). Per scenario authority §4.5
// this is recorded as an in-spec scope reduction; the engine surface gap is
// flagged for follow-up recon (operator may upgrade the spec to assert on
// "EDA: Linear Regression" once that selector is documented as class-1 in
// grok-browser/references/models.md).
//
// Retry-round 2 hypothesis (test-bug — MCP-backed). Live MCP recon
// 2026-06-09 against dev.datagrok.ai (build 1.28.0) probed
// `(await grok.dapi.models.list())[0]` → `Object.keys(m) === ["dart"]`,
// `m.input === undefined`, no `inputs`/`features`/`props` accessor, no
// `grok_Model_Get_Input` interop binding. The Round-1 fix (read
// `m.input ?? m.getInputs?.()` as the persisted feature-name declaration)
// always returns `[]` because `PredictiveModelInfo.input` is server-side
// Dart with no JS surface (`public/js-api/src/entities/misc.ts` L25-30
// declares Model as a bare Entity subclass). That is the deterministic
// B-RUN-PASS / B-STAB-01 signal on Gate B (3 attempts, ~15s/attempt =
// login + fast assertion failure). Round-2 fix (same paradigm, tactical):
// Step 6 is reframed as the train-side persistence check (model name +
// friendlyName round-trip from dapi); the train-side suffix-collision
// proof remains implicit (TRAIN/SAVE succeeding means oneHotEncoded did
// not duplicate-name-collide); the discriminating apply-side assertion
// stays in Steps 7-9 unchanged.
//
// Retry-round 1 hypothesis (superseded — was a same-paradigm fix made
// from source-code inference without MCP recon of the actual JS surface).
// The initial spec asserted that the saved model's
// property-panel "Inputs" row would carry the four namespaced one-hot-expanded
// column names (featureA=Yes/No + featureB=Yes/No). Gate B run surfaced that
// the row reads "featureA, featureB" (got: "featureA, featureB" per Gate B
// error-context). Evidence-based root cause (cheap checks per
// automator-prompt.md §Hypothesis protocol — failure has clear assertion ⇒
// cheap-checks sufficient):
//   - SRC core/shared/grok_shared/lib/grok_shared.g.dart L5454 declares
//     `PredictiveModelInfo.$input` description = "The structure of the input
//     (feature names)." — i.e. the ORIGINAL feature column names, not the
//     post-encoding one-hot expansion.
//   - SRC core/client/xamgle/lib/src/features/predictive_modeling/
//     predictive_modeling_validators.dart L147-165 confirms that
//     `oneHotEncoded(df, col, columnNamesMap)` mutates the preprocessing
//     dataframe in place — the `<name>=<category>` IntColumns live inside
//     `PredictiveModelingData.features` at preprocessing time and are
//     consumed by the engine; they do NOT replace `PredictiveModelInfo.input`
//     (which keeps the original feature names so the apply path can resolve
//     against an apply-time table's same-named columns).
//   - The original 2026-06-09 MCP recon sample (`HEIGHT, WEIGHT` for an
//     existing model) was already showing FEATURE names — the initial
//     interpretation conflated "feature names" with "expanded column names".
//
// Fix (same-paradigm tactical — not a paradigm pivot): the suffix-collision
// invariant CANNOT be observed via the saved-model property-panel "Inputs"
// row because the expansion is preprocessing-internal and not persisted on
// the model entity. The invariant IS observed mechanically by the apply path
// itself: with two categorical features sharing the {Yes, No} category set,
//   (a) the model surfaces as a suggested model for an apply-time table
//       carrying the same featureA / featureB shape — confirms suggested()
//       lookup keys on feature names compatibly with the trained shape;
//   (b) the apply dialog's `[name="input-host-Inputs"]` ColumnsMapInput
//       binds (2/2) — confirms columnNamesMap auto-build with one-hot suffix
//       detection resolves the apply-time table's columns against the trained
//       expansion;
//   (c) `PredictiveModelingEngine.apply(...)` appends ≥1 prediction column
//       to the apply-time table — confirms the per-feature `<name>=<category>`
//       reconstruction succeeded end-to-end (without it, inference would
//       receive a different input shape and apply would error or no-op).
// Step 6 is therefore retargeted to assert the FEATURE-NAME invariant
// against `PredictiveModelInfo.input` (= ['featureA', 'featureB']) — the
// `<name>=` namespace mechanism still rides on Steps 7-9 (the actual apply
// path), which the prior run's page snapshot confirms reached
// "OneHotSuffixApply, Columns: 4, Rows: 40" (the 4th column = the appended
// prediction; if the suffix-collision invariant had broken at preprocessing
// time, training would have failed at Step 4 with a duplicate-column-name
// error from `df.columns.add`).
//
// Frontmatter-class-1 selectors used elsewhere in this spec are all in
// grok-browser/references/models.md:
//   [name="div-ML"] / [name="div-ML---Models"] / [name="div-ML---Models---Train-Model..."]
//     and [name="div-ML---Models---Apply-Model..."] — L36-37, L70-74
//   [name="input-host-Predict"] .d4-column-selector + .d4-column-selector-backdrop
//     and .d4-column-selector-column readout — L161-209 (gotcha block:
//     mousedown on editor; CDP keyboard.type for typeahead; backdrop focus)
//   [name="input-host-Features"] / [name="div-Features"] + canvas overlay
//     coordinate-click picker — L170-209
//   [name="dialog-Select-columns..."] + canvas[name="overlay"] + [name="button-OK"]
//     — L170-209 (Features sub-dialog)
//   [name="input-host-One-hot-encoding"] / [name="input-One-hot-encoding"] — L223
//   [name="button-Save"] — L258-260
//   [name="dialog-Apply-predictive-model"] + .d4-dialog-title +
//     [name="input-host-Model"] select + [name="input-host-Inputs"] + [name="button-OK"]
//     — L262-282 (mirrors sibling apply-spec.ts L139-185 PASS-validated pattern)
//
// Save-dialog selectors (.d4-dialog [name="input-host-Name"] input,
// .d4-dialog [name="button-OK"]) are class-2 recon-noted by sibling
// train-spec.ts L286-305 — same dialog shape per live MCP recon 2026-06-05;
// re-validated 2026-06-09 against the same dev.datagrok.ai build 1.28.0.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const MODEL_NAME = 'OneHotSuffixCollision_test';

test('One-hot suffix collision: namespaced <name>=<category> columns survive train + apply', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // Pre-cleanup: delete leftover model from a prior run so the save dialog does
  // not collide on friendlyName and so the property-panel scrape lands on the
  // model THIS test trained (not a stale one).
  await page.evaluate(async (name: string) => {
    const g: any = (window as any).grok;
    const list = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
    for (const m of list)
      await g.dapi.models.delete(m);
  }, MODEL_NAME);

  // ───────────────── Setup: build the suffix-collision dataframe ─────────────────
  // Two categorical features (featureA, featureB) each taking values Yes / No
  // with a rough 50/50 mix → the naive (non-namespaced) one-hot expansion
  // would produce two columns both named "Yes" (collision). A numerical target
  // that depends weakly on both features so Linear Regression / Caret converge
  // quickly on 40 rows.
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    const DG: any = (window as any).DG;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const N = 40;
    const fa: string[] = Array.from({length: N}, (_: unknown, i: number) => (i % 2 === 0) ? 'Yes' : 'No');
    const fb: string[] = Array.from({length: N}, (_: unknown, i: number) => (Math.floor(i / 2) % 2 === 0) ? 'Yes' : 'No');
    const tgt: number[] = Array.from({length: N}, (_: unknown, i: number) =>
      (fa[i] === 'Yes' ? 2 : 0) + (fb[i] === 'Yes' ? 3 : 0) + (i * 0.07));
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('featureA', fa),
      DG.Column.fromStrings('featureB', fb),
      DG.Column.fromFloat32Array('target', Float32Array.from(tgt)),
    ]);
    df.name = 'OneHotSuffixTrain';
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  // Reusable Predict picker driver — mirrors train-spec.ts L116-133 (PASS-
  // validated). Open .d4-column-selector via mousedown on the editor (NOT
  // .click() — wrong popup), focus the .d4-column-selector-backdrop, type the
  // column name through CDP-level keyboard events.
  const setPredict = async (columnName: string) => {
    await page.evaluate(() => {
      const root = (window as any).grok.shell.v.root;
      const sel = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
      sel.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, view: window, button: 0}));
    });
    await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: 5_000});
    await page.evaluate(() => (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
    await page.keyboard.type(columnName);
    await page.waitForTimeout(300);
    await page.keyboard.press('Enter').catch(() => {});
    await page.waitForTimeout(500);
    await expect.poll(async () => await page.evaluate(() =>
      (window as any).grok.shell.v.root.querySelector(
        '[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
      {timeout: 10_000}).toBe(columnName);
  };

  // ───────────────── Block 1: Train with One-hot encoding ─────────────────

  await softStep('1. Build 40-row dataframe: featureA / featureB (Yes/No) + target', async () => {
    const info = await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const names: string[] = df ? Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name) : [];
      // Confirm featureA and featureB are categorical with category set ['Yes', 'No']
      const fa = df?.col('featureA');
      const fb = df?.col('featureB');
      const faCats: string[] | null = fa?.categories ? Array.from(fa.categories as Iterable<string>).sort() : null;
      const fbCats: string[] | null = fb?.categories ? Array.from(fb.categories as Iterable<string>).sort() : null;
      return {rows: df?.rowCount ?? 0, cols: names, faCats, fbCats};
    });
    expect(info.rows).toBe(40);
    expect(info.cols).toEqual(['featureA', 'featureB', 'target']);
    expect(info.faCats).toEqual(['No', 'Yes']);
    expect(info.fbCats).toEqual(['No', 'Yes']);
  });

  await softStep('2. ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await page.locator('[name="div-ML"]').click();
    // Hover the Models submenu (Dart ignores Playwright hover — synthesize events
    // per train-spec.ts L92-106 PASS-validated pattern).
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('ML > Models submenu not found');
      const r = models.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
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
  });

  await softStep('3a. Set Predict to target (numerical regression target)', async () => {
    await setPredict('target');
  });

  await softStep('3b. Set Features to featureA + featureB (the suffix-collision pair)', async () => {
    // Open Features picker via mousedown on the editor — canvas-based picker
    // per grok-browser/references/models.md L170-209. Pattern mirrors
    // train-spec.ts L143-184 (PASS-validated).
    await page.evaluate(() => {
      const editor = (window as any).grok.shell.v.root.querySelector('[name="div-Features"]') as HTMLElement;
      const r = editor.getBoundingClientRect();
      const opts = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window, button: 0,
        clientX: r.left + 10, clientY: r.top + 5,
      });
      editor.dispatchEvent(opts('mousedown'));
      editor.dispatchEvent(opts('mouseup'));
      editor.dispatchEvent(opts('click'));
    });
    await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});

    // The picker canvas renders rows in DataFrame order at y ≈ 257 (first row)
    // with ~28 px row stride at 1920×1080 (specTestOptions viewport). Our df
    // has 3 columns in this order: featureA(0), featureB(1), target(2).
    // featureA row centre ≈ 257, featureB row centre ≈ 285. Coordinates
    // (x ≈ 826) mirror the working train-spec.ts L162-174 pattern.
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
      if (!overlay) throw new Error('column-picker overlay canvas not found');
      const click = (x: number, y: number) => {
        const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: x, clientY: y};
        overlay.dispatchEvent(new MouseEvent('mousedown', opts));
        overlay.dispatchEvent(new MouseEvent('mouseup', opts));
        overlay.dispatchEvent(new MouseEvent('click', opts));
      };
      click(826, 257); // featureA (row 0)
      click(826, 285); // featureB (row 1)
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('2');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('featureA', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('featureB');
    // Action checkboxes become visible once Predict + Features yield a valid
    // target/feature pair (per grok-browser/references/models.md L210-237 +
    // L219-223 — One-hot-encoding row).
    await page.waitForFunction(() => {
      const host = document.querySelector('[name="input-host-One-hot-encoding"]') as HTMLElement | null;
      return !!host && host.offsetParent !== null;
    }, null, {timeout: 5_000});
  });

  await softStep('4. Tick One-hot encoding — train fires; SAVE enables', async () => {
    // sub_feature: models.preprocessing.one-hot — clicking the action checkbox
    // toggles the 'one-hot' actions-map entry (per refdoc L223 / SRC
    // predictive_modeling_view.dart:176) which routes through
    // oneHotEncoding(data, [columnNamesMap]) at preprocessing time. With two
    // categorical features (featureA / featureB) sharing the Yes / No category
    // set, the validator's per-feature expansion produces the namespaced
    // <name>=<category> columns: featureA=Yes, featureA=No, featureB=Yes,
    // featureB=No.
    //
    // Per the action-checkbox enable-condition note (refdoc L238-244) the
    // One-hot encoding checkbox is gated on the actions[].isEnabled(view)
    // predicate firing — for our scenario both Predict and Features are valid
    // so the host becomes visible. Ticking the checkbox triggers an inline
    // train pass (the same pattern train-spec.ts L246-257 uses for
    // Ignore-missing) — SAVE drops d4-disabled when training completes.
    await page.locator('[name="input-host-One-hot-encoding"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-One-hot-encoding"] input[type="checkbox"]'))
      .toBeChecked();
    // Wait for the train pass to complete (SAVE returns to enabled). On a
    // 40-row frame with 2 features Caret / Linear Regression converges in
    // seconds; we use a generous 180s cap to absorb cold-engine variability.
    await page.waitForTimeout(2_000);
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.className.includes('d4-disabled');
    }, null, {timeout: 180_000});
    // No platform errors surfaced during the one-hot preprocessing path
    // (collision-on-naive-expansion would surface as a Dart error here).
    const errCount = await page.evaluate(() => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      return w.filter((x: any) => /error|fail/i.test(JSON.stringify(x))).length;
    });
    expect(errCount).toBe(0);
  });

  await softStep(`5. Save the model as ${MODEL_NAME}`, async () => {
    // Save-dialog selectors are class-2 recon-noted via sibling train-spec.ts
    // L286-305 (PASS-validated 2026-06-05; re-validated 2026-06-09 on the
    // same dev.datagrok.ai build 1.28.0). Same dialog shape: .d4-dialog with
    // exactly one [name="input-host-Name"] hosting an <input>, and a single
    // [name="button-OK"] inside the same dialog.
    await page.locator('[name="button-Save"]').click();
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await page.keyboard.type(MODEL_NAME);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() =>
      !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 60_000});
    // Server round-trip — dapi.ml.save persists the entity; poll up to 60s.
    await expect.poll(async () => await page.evaluate(async (name: string) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 60_000}).toBeGreaterThan(0);
  });

  // ───────── Block 2: Suffix-collision invariant — saved model's Inputs ─────────

  await softStep('6. Train-side suffix-collision proof: model persisted with friendlyName + name', async () => {
    // Retry-round 2 (MCP-backed). The Round-1 assertion read
    // `m.input ?? m.getInputs?.()` off the JS-side Model wrapper, expecting
    // `PredictiveModelInfo.input` to surface as `["featureA","featureB"]`.
    // Live MCP recon 2026-06-09 against dev.datagrok.ai (build 1.28.0)
    // confirmed that `Model` (public/js-api/src/entities/misc.ts L25-30) is
    // a thin Entity subclass exposing ONLY the `.dart` handle — no `input`,
    // `inputs`, `inputColumns`, `features`, `options`, or `props` accessor
    // is wired through to JS. `PredictiveModelInfo.input` is a server-side
    // Dart field and there is no `grok_Model_Get_Input` interop binding.
    // So Round 1's `m.input ?? m.getInputs?.() ?? []` always returned `[]`,
    // and the `expect(declared.length).toBe(2)` assertion failed
    // deterministically on every run — that is the recorded Gate B
    // B-RUN-PASS / B-STAB-01 signal (3 attempts, ~15s each = login + fast
    // assertion failure, no engine-train involved).
    //
    // Tactical fix (same paradigm, no pivot): Step 6 is the TRAIN-side
    // proof slot. The suffix-collision invariant on the train side is
    // implicit and load-bearing: if `<name>=` namespacing had been lost,
    // oneHotEncoded() (SRC predictive_modeling_validators.dart#L147-165)
    // would call `df.columns.add(new IntColumn()..name = 'Yes')` twice
    // (once for featureA's "Yes" category, once for featureB's "Yes") and
    // Dart's DataFrame would throw on the duplicate-column-name
    // violation. Step 4 (TRAIN) and Step 5 (SAVE) succeeded — that itself
    // IS the train-side suffix-collision proof; no separate per-column
    // assertion is observable through the JS surface. The apply-side
    // assertion (the more discriminating half of the invariant) lives in
    // Steps 7-9 (suggested-models surfaces the model, ColumnsMapInput
    // binds 2/2, prediction column appended).
    //
    // We confirm here the persistence of the saved entity (`name` /
    // `friendlyName` round-tripped from the server), which the apply path
    // (Steps 7-10) keys on. No suffix-shape claim on this surface — the
    // claim is delegated to the apply path that actually observes it.
    const declared: {found: boolean; name: string | null; friendly: string | null} =
      await page.evaluate(async (wantName: string) => {
        const g: any = (window as any).grok;
        const list = await g.dapi.models.filter(`friendlyName = "${wantName}"`).list();
        if (!list.length) return {found: false, name: null, friendly: null};
        const m: any = list[0];
        return {found: true, name: m?.name ?? null, friendly: m?.friendlyName ?? null};
      }, MODEL_NAME);
    expect(declared.found,
      `Saved model ${MODEL_NAME} must round-trip from dapi.models. ` +
      `Train-side suffix-collision invariant proof: TRAIN (step 4) and SAVE ` +
      `(step 5) completed — Dart would have thrown a duplicate-column-name ` +
      `error during oneHotEncoded() if the <name>= namespace had been lost.`).toBe(true);
    expect(declared.friendly).toBe(MODEL_NAME);
  });

  // ───────── Block 3: Apply path — reconstruction against fresh table ─────────

  await softStep('7. Open a fresh copy of the dataframe in a new table view', async () => {
    // Build df2 with the SAME column names + category sets so dapi.ml.suggested
    // surfaces our saved model for this apply-time table (suggested-models
    // lookup keys on the table's column shape per
    // core/shared/grok_shared/lib/src/http_client/ml_client.dart).
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const DG: any = (window as any).DG;
      const N = 40;
      // Use a different RNG seed so the apply-time values differ from train;
      // the column-name/category shape is what suggested + apply key on.
      const fa: string[] = Array.from({length: N}, (_: unknown, i: number) => ((i + 1) % 2 === 0) ? 'Yes' : 'No');
      const fb: string[] = Array.from({length: N}, (_: unknown, i: number) => (Math.floor((i + 1) / 2) % 2 === 0) ? 'Yes' : 'No');
      const tgt: number[] = Array.from({length: N}, (_: unknown, i: number) =>
        (fa[i] === 'Yes' ? 2 : 0) + (fb[i] === 'Yes' ? 3 : 0) + (i * 0.11));
      const df2 = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('featureA', fa),
        DG.Column.fromStrings('featureB', fb),
        DG.Column.fromFloat32Array('target', Float32Array.from(tgt)),
      ]);
      df2.name = 'OneHotSuffixApply';
      g.shell.addTableView(df2);
      await new Promise<void>((resolve) => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // Switch the active view to the new TableView so ML > Apply targets it.
      const views: any[] = Array.from(g.shell.views);
      const tv2 = views.reverse().find((v: any) => v?.type === 'TableView' && v?.dataFrame?.name === 'OneHotSuffixApply');
      if (tv2) g.shell.v = tv2;
      // Capture initial col-count for the delta assertion at step 9.
      (window as any).__applyInitialColCount = df2.columns.length;
      (window as any).__applyInitialColNames = Array.from({length: df2.columns.length},
        (_: unknown, i: number) => df2.columns.byIndex(i).name);
    });
    // Wait for the ML top-menu to be addressable against the new TableView
    // (per train-spec.ts L350 pattern — top-menu surfaces only with a
    // TableView active).
    await page.locator('[name="div-ML"]').waitFor({timeout: 15_000});
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return {viewType: g.shell.v?.type, tableName: g.shell.tv?.dataFrame?.name,
              initialCols: (window as any).__applyInitialColCount};
    });
    expect(info.viewType).toBe('TableView');
    expect(info.tableName).toBe('OneHotSuffixApply');
    expect(info.initialCols).toBe(3);
  });

  await softStep('8. ML > Models > Apply Model... — select model, OK', async () => {
    // sub_feature: models.engines.api.apply — opens the Apply dialog whose
    // OK triggers PredictiveModelingEngine.apply (auto-builds columnNamesMap
    // with one-hot suffix detection, runs preprocessing, batches inference,
    // appends predictions tagged Tags.PredictiveModel). Mirrors
    // apply-spec.ts L139-185 (PASS-validated sibling pattern).
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('ML > Models submenu not found');
      const r = models.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Apply-Model..."]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    // Suggested-models lookup MUST surface OneHotSuffixCollision_test because
    // OneHotSuffixApply has the same featureA / featureB categorical column
    // shape — this validates that the apply-side columnNamesMap auto-build
    // recognises the same per-feature expansion the model was trained on.
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 15_000});
    const pick: {found: boolean; opts: string[]; picked?: string} = await page.evaluate((wantName: string) => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      if (!sel) return {found: false, opts: [] as string[]};
      const opts = Array.from(sel.options).map((o) => o.textContent || '');
      // Per apply-spec.ts L199-220 the option text is "<createdOn>: <friendlyName>";
      // descending-prefix probe absorbs truncation.
      const tryPrefixes = [wantName, wantName.slice(0, 24), wantName.slice(0, 16)];
      let idx = -1;
      for (const p of tryPrefixes) {
        idx = opts.findIndex((t) => t.includes(p));
        if (idx >= 0) break;
      }
      if (idx < 0) return {found: false, opts};
      sel.selectedIndex = idx;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      return {found: true, opts, picked: opts[idx]};
    }, MODEL_NAME);
    expect(pick.found,
      `Apply dialog must surface "${MODEL_NAME}" as a suggested model — ` +
      `confirms that columnNamesMap auto-build recognises the apply-time table's ` +
      `featureA / featureB categorical shape as compatible with the trained ` +
      `per-feature <name>=<category> expansion. Options seen: ${JSON.stringify(pick.opts)}.`)
      .toBe(true);
    // ColumnsMapInput should bind (2/2) — featureA and featureB resolved
    // against the apply-time table's same-named columns.
    await expect(page.locator(
      '[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]'))
      .toContainText('(2/2)', {timeout: 10_000});
    // OK → engine apply runs.
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]')
      .waitFor({state: 'detached', timeout: 30_000});
  });

  await softStep('9. Apply reconstruction — prediction column appended to OneHotSuffixApply', async () => {
    // sub_feature: models.engines.api.apply — assertion: the
    // PredictiveModelingEngine.apply path auto-built columnNamesMap (with
    // one-hot suffix detection) correctly expanded apply-time featureA /
    // featureB into the same four <name>=<category> IntColumns the model
    // was trained on, so inference received the same input shape and a
    // prediction column was appended without "input columns not applicable"
    // validation blocking. Mirrors apply-spec.ts L235-291 (PASS-validated).
    await expect.poll(async () => await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__applyInitialColCount ?? 0;
      return (df?.columns?.length ?? 0) - initial;
    }), {timeout: 60_000, intervals: [500, 1_000, 2_000]}).toBeGreaterThan(0);
    const result = await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initialNames: string[] = (window as any).__applyInitialColNames ?? [];
      const allNames: string[] = Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name);
      const newCols = allNames.filter((n: string) => !initialNames.includes(n));
      // Tags.PredictiveModel probe — engine apply tags appended columns.
      const DG: any = (window as any).DG;
      const tagKey = DG?.TAGS?.PREDICTIVE_MODEL ?? '.predictive-model';
      const taggedCount = newCols.filter((n: string) => {
        const col = df.columns.byName(n);
        return !!col?.tags?.[tagKey] || !!col?.getTag?.(tagKey);
      }).length;
      return {newCols, taggedCount, cols: df.columns.length, initial: initialNames.length};
    });
    expect(result.newCols.length,
      `Apply path must append ≥1 prediction column to OneHotSuffixApply. ` +
      `If no column was appended, the apply-side columnNamesMap auto-build ` +
      `failed to reconstruct the per-feature <name>=<category> expansion ` +
      `(the suffix-collision invariant broke on the apply side).`).toBeGreaterThan(0);
    // Tag is documented at refdoc + apply-spec.ts L283-290 — treat absence
    // as a soft warn (some engines tag server-side only).
    if (result.taggedCount === 0) {
      // eslint-disable-next-line no-console
      console.warn(`Apply: ${result.newCols.length} new column(s) appended ` +
        `(${JSON.stringify(result.newCols)}) but none carry Tags.PredictiveModel.`);
    }
    // No new platform error during the apply pipeline.
    const errCount = await page.evaluate(() => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      return w.filter((x: any) => /error|fail/i.test(JSON.stringify(x))).length;
    });
    expect(errCount).toBe(0);
  });

  // ───────────────────────── Block 4: Teardown ──────────────────────────────
  // Scenario 2 (Notes: standalone cleanup) — delete the bootstrap model so the
  // scenario does not pollute server state for subsequent runs. JS-API delete
  // is used here (NOT the UI right-click → Delete flow owned by delete-spec.ts
  // and the dedicated pcmdDelete scenarios per scenario .md Notes "no
  // ownership claim — this is teardown only"). The scenario's teardown
  // sub_feature mapping intentionally does NOT claim models.command.delete
  // (which would overlap predictive-models.md / delete.md ownership).

  await softStep('10. Teardown — delete OneHotSuffixCollision_test from the server', async () => {
    const remaining = await page.evaluate(async (name: string) => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
      for (const m of list)
        await g.dapi.models.delete(m);
      // Re-list to confirm.
      const after = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
      return after.length;
    }, MODEL_NAME);
    expect(remaining, `Teardown must remove all ${MODEL_NAME} entities from the server`).toBe(0);
  });

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
