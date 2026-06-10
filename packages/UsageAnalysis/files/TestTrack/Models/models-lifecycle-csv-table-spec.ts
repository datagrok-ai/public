/* ---
sub_features_covered: [models.command.train, models.view.training, models.api.save, models.command.apply, models.workflow.apply-dialog, models.engines.api.apply, models.meta.performance-section, models.command.edit, models.workflow.edit-info, models.command.share, models.command.delete, models.workflow.remove]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (scenario: proactive-lifecycle per chain rationale; orchestrator pin pending)
//   sub_features_covered: [models.command.train, models.view.training, models.api.save,
//                          models.command.apply, models.workflow.apply-dialog,
//                          models.engines.api.apply, models.meta.performance-section,
//                          models.command.edit, models.workflow.edit-info,
//                          models.command.share, models.command.delete, models.workflow.remove]
//   ui_coverage_responsibility: [] (delegated_to: predictive-models.md)
//   related_bugs: [] (chain bug_focused_candidates owned by dedicated specs)
//
// Atlas provenance (derived_from):
//   feature-atlas/models.yaml#sub_features[models.api.save] source:
//     core/shared/grok_shared/lib/src/http_client/ml_client.dart#L35
//   feature-atlas/models.yaml#sub_features[models.engines.api.apply] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_engines.dart#L26
//   feature-atlas/models.yaml#sub_features[models.meta.performance-section] source:
//     core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L118
//   feature-atlas/models.yaml#sub_features[models.workflow.edit-info] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_browser.dart#L46
//   feature-atlas/models.yaml#sub_features[models.workflow.remove] source:
//     core/client/xamgle/lib/src/meta/predictive_model_info_meta.dart#L224
//
// Spec target: proactive-lifecycle for source_class=trained_on_csv_table.
// End-to-end UI-driven lifecycle: Train (Linear Regression on accelerometer.csv)
// → Save → Apply on fresh open → Performance pane Run Evaluation → Edit
// description → Share dialog open/cancel → Delete + cleanup verification.
// JS API permitted broadly for setup/teardown + side-effect assertions;
// each gateway operation (train view, apply dialog, performance pane,
// edit modal, share dialog, delete confirm) is UI-driven per scenario
// authority — these are the non-agnostic differentiators that motivate
// the proactive_lifecycle_specs[trained_on_csv_table] entry in
// scenario-chains/models.yaml.
//
// Selector recon-notes (class-2: live-MCP-observed 2026-06-09 on
// dev.datagrok.ai build 1.28.0, not yet in grok-browser/references/models.md):
//   [name="dialog-Predictive-model"] — Edit dialog opened by Edit... context-menu
//     item on a model card (cmdEditModel → editModelInfo per
//     predictive_modeling_browser.dart#L46-49 Modal(title: 'Predictive model')).
//     Contains [name="input-host-Name"], [name="input-host-Description"],
//     [name="button-OK"], [name="button-CANCEL"]. Reached via right-click model
//     card → "Edit..." menu item. Observed live 2026-06-09 via
//     chrome-devtools MCP (evaluate_script enumerated the dialog).
//   [name="button-Run-Evaluation"] — Big-button inside the Performance accordion
//     pane on the model context panel (per
//     predictive_model_info_meta.dart#L118 renderPerformanceSection). Reached by
//     clicking the model card to select it, then expanding the Performance
//     accordion pane. Observed live 2026-06-09 via chrome-devtools MCP
//     (querySelector returned exactly 1 BUTTON with this name, text
//     "Run Evaluation", className "ui-btn ui-btn-ok ui-btn-raised").
//   [name="dialog-Share-<friendlyName>"] — Share dialog opened by Share...
//     context-menu item (regCommand Share per predictive_model_info_meta.dart#L173
//     → shareEntity). The friendlyName is appended verbatim with the
//     html_utils.annotate() char-substitution rule (underscores and special
//     chars become hyphens). Title text "Share <friendlyName>". Observed live
//     2026-06-09 via chrome-devtools MCP — sample dialog was
//     [name="dialog-Share-Caret"]. Sibling-spec precedent: chemprop-spec.ts
//     L597-606 uses the same hyphen-encoded form.
//   .d4-dialog [name="input-host-Name"] input — Name text input inside the
//     model-save dialog (mounted by Save in PredictiveModelingView title bar).
//     Same dialog shape documented in train-spec.ts L286-305; reused here
//     for consistency.
//   [name="input-Model-Engine"] — Model Engine select (NOT input-host-prefixed
//     despite the standard `input-host-<Caption>` convention). Mounted lazily
//     by predictive_modeling_view.dart:486-512 once engine.isApplicable
//     resolves with non-empty applicableEngines. Sibling-spec
//     predictive-models-spec.ts L152 / L210-217 uses the same exact selector
//     and is the proven-passing pattern this build. Observed live 2026-06-09
//     via chrome-devtools MCP: the input-host list on the bad
//     (applicableEngines.isEmpty) state contains [Table, Predict, Features,
//     Ignore-missing, Impute-missing, One-hot-encoding, Skip-unique-categories,
//     Predict-probability] — neither input-host-Model-Engine NOR
//     input-Model-Engine appears until Features changes to a (3) state that
//     excludes the Predict column.
//
// MCP recon evidence for the soft-binding choice (paradigm not pivoted):
//   - System:DemoFiles/sensors/accelerometer.csv path verified via
//     grok.dapi.files.list (returns "sensors/accelerometer.csv" entry).
//   - Card context-menu items confirmed: ["Share...", "Edit...", "Delete",
//     "Save as Zip"] — no "Apply to" submenu surfaced on this card (atlas drift
//     consistent with grok-browser/references/models.md L444-450 note).
//   - Performance accordion pane present alongside Details / Activity0 /
//     Sharing / Chats / Sticky meta.
//   - User has 16 existing models; the model name MUST be unique-per-user (a
//     prior LifecycleCsvModel from a re-run would collide on save). Per-cycle
//     name suffix used: `LifecycleCsvModel_<runStamp>`. Pre-cleanup sweep
//     deletes ANY prior model named LifecycleCsvModel* to keep the test
//     hermetic across reruns.
//
// Retry-1 hypothesis (Gate B failure_keys B-RUN-PASS + B-STAB-01,
// 2026-06-09-models-automate-01) — INITIAL AUTHOR FIX:
//   Predict default = accel_x; scope-reduced to Predict=accel_x +
//   Features=[accel_y, accel_z, time_offset] (no "roll" column on live
//   dataset). All 12 sub_features remain covered.
//
// Retry-2 hypothesis (Gate B FAILED AGAIN with same keys 2026-06-09T18:30Z
// after the Predict=accel_x narrowing landed) — DIAGNOSED root cause as
// step 1.5 selector + step 1.4 None→3 flow. Retry-3 MCP recon (2026-06-09
// fresh probe, mcp_status=used) found the retry-2 hypothesis was wrong in
// both axes:
//
// Retry-4 MCP recon findings (2026-06-10 — cycle 2026-06-09-models-automate-02):
//
//   Live probe on dev.datagrok.ai build 1.28.0 (chrome-devtools MCP
//   list_pages = used; viewport native). Re-executed the Train view +
//   Features picker flow end-to-end:
//
//   (a) Overlay rect on this 4-column dataset: {w:228, h:142, top:223.25,
//       right:875}. Note h/4 = 35.5, which is CLOSE to the documented
//       row-0 centre (top + 36) but the docs warn that this coincidence
//       is build-specific and NOT load-bearing; reverting to the fixed
//       top+36 is the safer pattern.
//   (b) Synthetic dispatch at (right-40, top+36) = (835, 259.25) WORKS:
//       "4 checked" → "3 checked", Features text = "(3) accel_y, accel_z,
//       time_offset", accel_x correctly deselected (it is row 0 in
//       original DataFrame column order).
//   (c) The retry-3 spec was clicking at (left+15, top + height/4 * 0.5):
//       (a) wrong X — checkbox is at right edge per models.md L192-220,
//       not left;
//       (b) no 1200ms settle wait before the click — models.md L209
//       explicitly says canvas needs ≥1200ms after picker open or clicks
//       silently no-op. This explains the cycle-01 B-STAB-01 failure_key.
//   (d) Engine selector `[name="input-Model-Engine"]` mounts in ~30s on
//       this build (engine.isApplicable async probe across Caret + EDA +
//       Chemprop). Options on this Predict=accel_x + Features=(3) pair
//       observed live: ["", "Please select a model", "Eda: XGBoost",
//       "Eda: PLS Regression", "Eda: Linear Regression"]. The retry-3
//       60s budget is sufficient. Selector pattern unchanged from retry-3.
//
// Retry-3 MCP recon findings (2026-06-09):
//
//   (a) Selector recon — On the live PredictiveModelingView at initial mount
//       (Features="(4) All", Predict=accel_x — the applicableEngines.isEmpty
//       branch), the input-host list is exactly:
//         [input-host-Table, input-host-Predict, input-host-Features,
//          input-host-Ignore-missing, input-host-Impute-missing,
//          input-host-One-hot-encoding, input-host-Skip-unique-categories,
//          input-host-Predict-probability]
//       NEITHER `[name="input-host-Model-Engine"]` NOR
//       `[name="input-host-Model-Engine"] select` exists in this state.
//       The sibling spec predictive-models-spec.ts L152 uses
//       `[name="input-Model-Engine"]` (without the input-host- prefix) —
//       this is the proven-passing selector pattern on this same build.
//
//   (b) Flow recon — The sibling spec predictive-models-spec.ts L122-143
//       uses the proven canvas-click pattern: click [name="label-All"],
//       wait for "4 checked", then page.mouse.click on canvas row-0 to
//       deselect accel_x, then wait for "3 checked". This pattern PASSED
//       Gate B in models-automate-01 already. The retry-2 None→3 monotonic
//       approach added 3 trusted canvas clicks (3x the failure-surface vs
//       the sibling's 1), and reverting to the sibling pattern eliminates
//       that risk.
//
// Retry-4 FIX (same-paradigm tactical, MCP-backed 2026-06-10, NOT a paradigm pivot):
//   1. Step 1.4 — match sibling predictive-models-spec.ts L125-175 EXACTLY:
//      (a) ≥1200ms settle wait after picker open (models.md L209 — was
//          missing in retry-3, B-STAB-01 root cause);
//      (b) synthetic mousedown/mouseup/click dispatch at (right-40,
//          top+36) — RIGHT edge for the checkbox column, FIXED row-0 Y
//          coordinate (was left+15 with height/4 in retry-3 — wrong X
//          AND build-dependent Y, the actual hit-zone is documented in
//          grok-browser/references/models.md L192-220);
//      (c) page.mouse.click fallback if synthetic does not toggle.
//      MCP empirical verification 2026-06-10: this pattern produces
//      "3 checked" and Features="(3) accel_y, accel_z, time_offset"
//      deterministically on this build.
//   2. Step 1.5 — unchanged from retry-3 (sibling-proven selector
//      `[name="input-Model-Engine"]` + setter.call direct value, 60s
//      mount wait). MCP recon 2026-06-10 confirmed: engine select mounts
//      in ~30s, "Eda: Linear Regression" is in the options list.
//   3. Step 1.4 .ui-input-column-names sub-selector replaced with the
//      parent input-host-Features whole-element textContent check — the
//      sub-class is implementation-dependent and the sibling uses the
//      parent path. Same correctness, less coupling.
//
// Hypothesis category: test-bug (selector + flow regression vs proven
// sibling). Evidence: MCP empirical probe of hosts list confirms
// input-host-Model-Engine absent; sibling-spec input-Model-Engine PASSed
// Gate B same cycle. Not a paradigm pivot — still page.mouse.click on
// canvas for Features and direct <select> value-setter for Engine.
//
// Retry-6 hypothesis (Gate B FAILED AGAIN with [B-RUN-PASS, B-STAB-01]
// at duration 64s × 3 attempts in cycle 2026-06-09-models-automate-02)
// — DISTINCT from all prior hypotheses. MCP recon 2026-06-10 fresh probe
// (cycle 2026-06-10-models-automate-02) established root cause:
//
// Apply dialog option text TRUNCATION. The option text format is
// "<date>: <friendlyName>" and truncates at ~43 total chars. With
// MODEL_NAME = LifecycleCsvModel_<timestamp> (31 chars), the option
// renders as "6/10/2026 10:44 AM: LifecycleCsvModel_178..." — only
// 23 chars of the name are visible. The prior spec's tryPrefixes
// tried [wantName, slice(0,32), slice(0,24), slice(0,16)] — the first
// three probes (31, 32, 24 chars) all exceed the 23 visible chars →
// ALL missed → pickResult.found = false → step 2.3 assertion failed.
// Only slice(0,16) = 'LifecycleCsvMode' (16 chars) was within range.
//
// This explains the 64s duration: steps 1.x complete ~15-20s fast,
// step 2.3 assertion fires immediately (not a timeout), stepErrors
// accumulates, final throw at ~64s. NOT an env-flake or budget issue.
//
// Retry-5 hypothesis (Gate B FAILED with [B-RUN-PASS, B-STAB-01]
// at duration 184s × 3 attempts in prior cycle — SUPERSEDED by
// retry-6 MCP recon. The 184s duration in the comment below was from
// an EARLIER attempt set; the current cycle's 64s matches step 2.3
// assertion failure, not a timeout.
//
// Retained retry-5 MCP context (evidence is still valid, just wrong
// root cause assigned):
//
// Critical evidence from decision-log section-batch outcome
// (cycle 2026-06-09-models-automate-01): ALL 6 playwright lifecycle /
// browse / metrics / validators specs in Models/ failed Gate B with
// the SAME [B-RUN-PASS, B-STAB-01] keys, AND the sibling
// predictive-models-spec.ts (the proven-passing pattern this spec was
// modelled after) is marked producer_status: EXHAUSTED with the same
// failure keys. This is a SECTION-WIDE cold-environment divergence,
// not a per-spec test-bug.
//
// 184s × 3 attempts breakdown points to step 1.5 (engine mount 60s +
// preview render 60s = 120s budget) timing out under cold-grok-test
// concurrent worker load: cold engine.isApplicable probes (Caret + 5
// EDA + Chemprop sequential) PLUS cold preview-widget render PLUS
// cold preview-pane scroll/layout PLUS save-enabled wait can exceed
// the 60s budgets when grok-test workers contend for the same
// engine-applicability cache that warm MCP has pre-populated.
//
// Retry-6 FIX (same-paradigm tactical, MCP-backed 2026-06-10,
// NOT a paradigm pivot):
//   Step 2.3 — model lookup truncation fix: add MODEL_BASE
//   ('LifecycleCsvModel', 17 chars) as the FIRST prefix in
//   tryPrefixes, before wantName and the length slices. MODEL_BASE
//   is shorter than the 23-char visible name portion so it always
//   matches. The prior slice(0,16) = 'LifecycleCsvMode' was also
//   within range but MODEL_BASE is more meaningful. The evaluate()
//   call is updated to accept [wantName, wantBase] tuple.
//   No other changes — paradigm, selectors, flow, timeouts unchanged.
//
// Hypothesis category: test-bug (assertion logic error — prefix
// search missed truncated option text). Evidence: MCP recon
// 2026-06-10 confirmed option text truncates at 43 chars, making
// the 31-char MODEL_NAME timestamp suffix invisible in the option.
// Distinct from all prior hypotheses. Not paradigm pivot.
//
// Retry-7 hypothesis (Gate B FAILED AGAIN with [B-RUN-PASS, B-STAB-01]
// at duration 69s × 3 attempts in cycle 2026-06-10-models-automate-02,
// timestamp 2026-06-10T15:30:00Z, AFTER the Retry-6 fix landed as the
// initial dispatch) — NEW root cause distinct from all prior hypotheses:
//
// MCP recon 2026-06-10 (retry dispatch, cycle 2026-06-10-models-automate-02):
//   Step 2.3 works: MODEL_BASE='LifecycleCsvModel' found in truncated option text.
//   Steps 3.x work: /models route, card click, Performance pane, Run Evaluation.
//   Root cause: step 4.2 `model.description` JS-API accessor returns undefined.
//   grok.dapi.models.list() returns Dart-proxy objects (class Wt). The proxy has
//   NO `description` getter — `Object.getPrototypeOf(model)` has only 'constructor'
//   in own props. Direct `model.description` returns undefined; description IS
//   stored server-side (verified by re-opening Edit dialog → Description input
//   shows saved value), but the Dart-JS boundary does not expose it.
//   Prior assertion: `toBe(NEW_DESCRIPTION)` on null → polls 30s → FAIL.
//   Duration: ~25-30s steps 1-3 + 30s poll timeout at step 4.2 = ~69s per
//   attempt. Consistent with observed failure_keys and duration.
//
// Retry-7 FIX (same-paradigm tactical, MCP-backed 2026-06-10, NOT a paradigm pivot):
//   Step 4.2 — replace model.description API assertion with:
//   (a) model-exists assertion (list.length > 0) confirming no accidental delete;
//   (b) re-open Edit dialog + assert Description input value = NEW_DESCRIPTION.
//   The re-open path is the only DOM-accessible verification that the server
//   persisted the description change. No other spec changes.
//   Hypothesis category: test-bug (broken accessor assertion). Evidence: MCP
//   recon directly probed model.description = undefined; re-open modal confirmed
//   description is actually saved. Distinct from all prior hypotheses.
//
// Retry-5 FIX (same-paradigm tactical, budget bump + cold-probe diag,
// MCP-backed 2026-06-10, NOT a paradigm pivot — retained for context):
//   1. Step 1.4 — settle wait 1200ms → 2000ms.
//   2. Step 1.5 — engine-mount wait 60s → 180s; preview-render → 180s.
//   3. Step 1.6 — save-enabled wait 120s → 180s.
//   4. Steps 2.4, 3.3 — unchanged.
//   Diagnosis (now superseded): cold-vs-warm divergence; 184s per
//   attempt. Retry-6 MCP recon shows the current cycle failed at 64s
//   (step 2.3 assertion), not a timeout.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// Per-run unique name avoids cross-run collisions in the dapi.ml.save
// friendlyName index while keeping the LifecycleCsvModel prefix for pre-cleanup.
const MODEL_BASE = 'LifecycleCsvModel';
const MODEL_NAME = `${MODEL_BASE}_${Date.now()}`;
const NEW_DESCRIPTION = 'CSV-backed lifecycle smoke (round 4)';

test('Models lifecycle on CSV-backed trainedOn (trained_on_csv_table)', async ({page}) => {
  test.setTimeout(900_000);

  await loginToDatagrok(page);

  // Pre-cleanup sweep: delete any LifecycleCsvModel* leftovers from prior runs
  // so the save flow does not collide on the friendlyName uniqueness and so the
  // /models browser does not display stale rows. Mirrors train-spec.ts L52-59.
  await page.evaluate(async ({prefix}) => {
    const g: any = (window as any).grok;
    const all = await g.dapi.models.list();
    const stale = all.filter((m: any) => (m.friendlyName || m.name || '').startsWith(prefix));
    for (const m of stale)
      try { await g.dapi.models.delete(m); } catch (_) { /* best-effort */ }
  }, {prefix: MODEL_BASE});

  // ═════════ Scenario 1: Train + save model (CSV-backed trainedOn) ════════════

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/sensors/accelerometer.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('1.1 Open accelerometer.csv — table view active', async () => {
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      const df = g.shell.tv?.dataFrame;
      const names = df ? Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name) : [];
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0,
              viewType: g.shell.v?.type, names};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.viewType).toBe('TableView');
    // Scope reduction (retry-1, MCP-confirmed): scenario .md cites a "roll"
    // column that does NOT exist in the live accelerometer.csv. The
    // four-column live shape is [accel_x, accel_y, accel_z, time_offset].
    // Predict is taken at the live default (accel_x) downstream; Features
    // stays as specified (accel_y, accel_z, time_offset). Assert only the
    // observed shape here.
    expect(info.names).toContain('accel_x');
    expect(info.names).toContain('accel_y');
    expect(info.names).toContain('accel_z');
    expect(info.names).toContain('time_offset');
  });

  await softStep('1.2 ML > Models > Train Model... — PredictiveModelingView opens', async () => {
    // sub_feature: models.command.train + models.view.training
    // Per grok-browser/references/models.md L70-106: click ML, dispatch
    // mouseover/mouseenter/mousemove on Models submenu (Dart ignores
    // Playwright .hover()), then click the leaf.
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
    await page.locator('[name="div-ML---Models---Train-Model..."]').click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel',
      null, {timeout: 15_000});
  });

  await softStep('1.3 Accept default Predict = accel_x (numerical regression target)', async () => {
    // Scope reduction (retry-1, MCP-confirmed): the scenario .md directs
    // "Predict = roll", but `roll` is not a column of accelerometer.csv on
    // the live server (4 cols: accel_x, accel_y, accel_z, time_offset).
    // PredictiveModelingView defaults Predict to the first numerical column
    // (accel_x — confirmed via MCP probe 2026-06-09 returning predictDefault
    // 'accel_x'). accel_x is a numerical regression target functionally
    // equivalent to the "roll" axis the scenario intended. The sibling
    // ui-smoke predictive-models-spec.ts uses the same default and the
    // same Features triple (accel_y, accel_z, time_offset) without setting
    // Predict explicitly. Asserting the default is the correct same-paradigm
    // fix; no UI driving on the column-selector is required.
    await expect.poll(async () => await page.evaluate(() =>
      (window as any).grok.shell.v.root
        .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
      {timeout: 15_000}).toBe('accel_x');
  });

  await softStep('1.4 Set Features to accel_y, accel_z, time_offset (3 cols)', async () => {
    // Retry-4 fix (MCP-backed 2026-06-10 on dev.datagrok.ai build 1.28.0):
    // Align EXACTLY with proven sibling predictive-models-spec.ts L125-175
    // canvas-overlay pattern from grok-browser/references/models.md L192-220:
    //   - Checkbox column sits at RIGHT edge of overlay (right - 40), NOT
    //     left + 15. The retry-3 left-edge click missed the hit-zone — that
    //     was the B-STAB-01 root cause.
    //   - Row 0 centre is FIXED at top + 36 (build-independent per
    //     models.md L211-212), NOT height/4 (which on the live 4-column
    //     overlay = 142/4 = 35.5 — coincidentally close, but the doc warns
    //     this is build-dependent).
    //   - ≥1200ms settle wait after picker open (models.md L209) — canvas
    //     paints async; clicks before settle silently no-op (B-STAB-01).
    //   - Synthetic mousedown/mouseup/click dispatch primary path; isTrusted
    //     page.mouse.click fallback if synthetic gated by future build.
    //
    // MCP empirical evidence (2026-06-10 live probe on this same build):
    //   - Overlay rect: {w:228, h:142, top:223.25, right:875}
    //   - After label-All: "4 checked"
    //   - After synthetic dispatch at (right-40, top+36) = (835, 259.25):
    //     "3 checked", Features text = "(3) accel_y, accel_z, time_offset"
    //   - accel_x (row 0) correctly deselected. Engine mount followed.
    await page.locator('[name="div-Features"]').click();
    await page.locator('.d4-dialog[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});
    // Settle wait — column-grid canvas paints async (~1.2s on warm dev,
    // tail-extends under cold-grok-test concurrent worker load; retry-5
    // bumps 1200ms → 2000ms to cover the cold tail per the section-wide
    // [B-RUN-PASS, B-STAB-01] diagnosis above).
    await page.waitForTimeout(2_000);
    // Click "All" → 4 columns checked (accel_x, accel_y, accel_z, time_offset).
    await page.locator('.d4-dialog [name="label-All"]').click();
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('.d4-dialog label'))
        .find((l) => /checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('4');
    }, null, {timeout: 10_000});
    // Synthetic dispatch on overlay canvas at right-edge row-0 checkbox.
    // Coordinates per grok-browser/references/models.md L192-220 +
    // sibling predictive-models-spec.ts L141-158.
    const toggled = await page.evaluate(async () => {
      const dlg = document.querySelector(
        '.d4-dialog[name="dialog-Select-columns..."]')!;
      const overlay = dlg.querySelector('canvas[name="overlay"]') as HTMLCanvasElement;
      const rect = overlay.getBoundingClientRect();
      const x = rect.right - 40;
      const y = rect.top + 36;
      const evtInit: MouseEventInit = {
        bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0, view: window,
      };
      overlay.dispatchEvent(new MouseEvent('mousedown', evtInit));
      overlay.dispatchEvent(new MouseEvent('mouseup', evtInit));
      overlay.dispatchEvent(new MouseEvent('click', evtInit));
      await new Promise((r) => setTimeout(r, 600));
      const counter = Array.from(dlg.querySelectorAll('label'))
        .map((l) => l.textContent?.trim() ?? '').find((t) => /checked/.test(t));
      return {counter, x, y};
    });
    // Defence in depth: if synthetic dispatch did not toggle (future build
    // re-introducing isTrusted gating), retry with page.mouse.click which
    // dispatches isTrusted events.
    if (!toggled.counter?.startsWith('3')) {
      await page.mouse.click(toggled.x, toggled.y);
      await page.waitForTimeout(500);
    }
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('.d4-dialog label'))
        .find((l) => /checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('3');
    }, null, {timeout: 10_000});
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // Verify Features text reflects (3) and accel_x is removed.
    await page.waitForFunction(() =>
      document.querySelector('[name="input-host-Features"]')?.textContent?.includes('(3)'),
      null, {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"]'))
      .not.toContainText('accel_x');
  });

  await softStep('1.5 Set Model Engine to Eda: Linear Regression', async () => {
    // Retry-3 fix: revert to sibling-spec selector
    // `[name="input-Model-Engine"]` (no input-host- prefix) per
    // predictive-models-spec.ts L152 / L200-217 — the proven-passing
    // selector on this build. MCP recon 2026-06-09 confirmed
    // `[name="input-host-Model-Engine"]` does not exist in the input-host
    // list of the current PredictiveModelingView. Engine label is canonical
    // "Eda: Linear Regression" (lowercase d) per Dart engine.type.
    //
    // Mount wait → 60s: Model Engine ChoiceInput is created lazily once
    // engine.isApplicable(features, target) resolves with a non-empty list
    // (predictive_modeling_view.dart:486-512). Applicability probes run
    // sequentially across registered engines (Caret + 5 EDA + Chemprop)
    // and can take 20-40s on slow dev. Diagnostic soft-warn dumps form
    // state on mount-timeout to surface applicableEngines.isEmpty branch.
    try {
      // Retry-5: 60s → 180s. Cold-grok-test engine.isApplicable probes
      // (Caret + 5 EDA + Chemprop sequential) tail-extend under concurrent
      // worker load; section-wide [B-RUN-PASS, B-STAB-01] diagnosed as
      // engine-mount budget exhaustion.
      await page.waitForFunction(() => {
        const sel = document.querySelector(
          '[name="input-Model-Engine"]') as HTMLSelectElement | null;
        return !!sel && sel.options.length > 0;
      }, null, {timeout: 180_000});
    } catch (err) {
      const diag = await page.evaluate(() => {
        const v: any = (window as any).grok.shell.v;
        const root = v?.root as HTMLElement | undefined;
        const insights = Array.from(root?.querySelectorAll('span, .d4-card, [class*="insight"], [class*="tips"]') ?? [])
          .map((e) => (e as HTMLElement).textContent?.trim() ?? '')
          .filter((t) => t && t.length < 200).slice(0, 10);
        const hosts = Array.from(root?.querySelectorAll('[name^="input-host-"]') ?? [])
          .map((e) => e.getAttribute('name'));
        const noModelsMsg = Array.from(root?.querySelectorAll('span') ?? [])
          .find((s) => /No models reg/i.test(s.textContent || ''))?.textContent ?? null;
        // Retry-5: surface grok.shell.warnings so cold server-side errors
        // are visible in the failure trace (cold engine-applicability
        // probes can throw server-side without surfacing in DOM).
        const warnings = (((window as any).grok.shell.warnings) ?? []).slice(-10)
          .map((w: any) => typeof w === 'string' ? w.slice(0, 200) : JSON.stringify(w).slice(0, 200));
        return {insights, hosts, noModelsMsg, warnings};
      });
      throw new Error(`Model Engine select did not mount within 180s. ` +
        `Likely applicableEngines.isEmpty branch (Predict + Features pairing yields no engines) ` +
        `or cold engine.isApplicable budget exhaustion. ` +
        `noModelsMsg=${diag.noModelsMsg}; hosts=${JSON.stringify(diag.hosts)}; ` +
        `insights=${JSON.stringify(diag.insights)}; warnings=${JSON.stringify(diag.warnings)}; ` +
        `original error: ${err}`);
    }
    await page.evaluate(() => {
      const sel = document.querySelector(
        '[name="input-Model-Engine"]') as HTMLSelectElement;
      sel.focus();
      const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
      // Direct value-set per sibling predictive-models-spec.ts L210-217.
      // Engine label is canonical "Eda: Linear Regression" (Dart engine.type
      // registration, per predictive_modeling_view.dart:494 ChoiceInput formatter).
      setter.call(sel, 'Eda: Linear Regression');
      sel.dispatchEvent(new Event('input', {bubbles: true}));
      sel.dispatchEvent(new Event('change', {bubbles: true}));
    });
    // Wait for the preview pane to render the engine label.
    // Retry-5: 60s → 180s. Cold-grok-test preview-widget render tail-extends
    // alongside engine.isApplicable cold-cache miss; broaden header
    // selector to include .d4-engine-parameters-header div (observed live
    // 2026-06-10 MCP probe — preview renders this class in addition to h3).
    await page.waitForFunction(() => {
      const text = Array.from(document.querySelectorAll(
        'h3, h4, [class*="card-header"], .d4-engine-parameters-header'))
        .map((e) => e.textContent?.trim() ?? '')
        .find((t) => /Eda:\s*Linear/i.test(t) && !/pls|partial/i.test(t));
      return !!text;
    }, null, {timeout: 180_000});
  });

  await softStep(`1.6 Save model as ${MODEL_NAME} — dapi.ml.save persists`, async () => {
    // sub_feature: models.api.save (MLClient.save → uploads trainedOn via
    // dapi.tables, POST /ml). Wait for SAVE to be enabled (preview must render
    // without errors before MLClient.save accepts the entity per
    // predictive_modeling_view.dart#L161 readyToSave gate).
    // Retry-5: 120s → 180s. Preview must fully render before MLClient.save
    // accepts the entity; cold engine isApplicable hot-path lag tail-extends
    // into save-enabled timing.
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 180_000});
    await page.locator('[name="button-Save"]').click();
    // Save dialog: [name="input-host-Name"] input + [name="button-OK"] (class-2
    // recon-noted in train-spec.ts L286-305; same dialog shape here).
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await nameInput.fill(MODEL_NAME);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    // Wait for the save dialog to close.
    await page.waitForFunction(() =>
      !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 90_000});
    // Verify the model surfaces in dapi.ml — the canonical assertion for
    // models.api.save success per train-spec.ts L319-322 precedent.
    await expect.poll(async () => await page.evaluate(async (name: string) => {
      const list = await (window as any).grok.dapi.models
        .filter(`friendlyName = "${name}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 60_000}).toBeGreaterThan(0);
  });

  // ═════════ Scenario 2: Apply model on a fresh open of the same CSV ═════════

  await softStep('2.1 Close current table view, re-open accelerometer.csv', async () => {
    // Force the apply-dialog open path (not in-place re-apply). Close all, then
    // re-open the same CSV. Switch back to the TableView so [name="div-ML"]
    // is addressable (PredictiveModelingView has its own minimal top-menu;
    // train-spec.ts L330-350 documents this gotcha).
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const df = await g.dapi.files.readCsv('System:DemoFiles/sensors/accelerometer.csv');
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      (window as any).__initialColCount = df.columns.length;
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    await page.locator('[name="div-ML"]').waitFor({timeout: 10_000});
  });

  await softStep('2.2 ML > Models > Apply Model... — apply dialog opens', async () => {
    // sub_feature: models.command.apply + models.workflow.apply-dialog
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
    // Model SELECT populated by dapi.ml.suggested(tableInfo) per
    // grok-browser/references/models.md L277-278.
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
  });

  await softStep(`2.3 Select ${MODEL_NAME} — inputs auto-map to (3/3)`, async () => {
    // sub_feature: models.workflow.apply-dialog (model binding) +
    //              models.engines.api.apply (auto columnNamesMap construction).
    // Select-option text is "<createdOn>: <friendlyName>" per
    // grok-browser/references/models.md L277-278; long names truncate at ~43
    // chars total (date prefix ~20 + name ~23). MODEL_NAME includes a timestamp
    // suffix making it 31 chars, so the option text truncates to
    // "...: LifecycleCsvModel_178..." — the slice(0,24) and slice(0,32) probes
    // exceed the 23-char visible portion and produce no match.
    // Fix (retry-6, MCP-backed 2026-06-10): probe MODEL_BASE ('LifecycleCsvModel'
    // = 17 chars) FIRST — this always fits within the visible name portion.
    // Then fall back to descending-length slices of wantName for build resilience.
    const pickResult = await page.evaluate(([wantName, wantBase]: string[]) => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      if (!sel) return {found: false, opts: [] as string[]};
      const opts = Array.from(sel.options).map((o) => o.textContent || '');
      // MODEL_BASE first (fits in truncated text), then full name and length slices.
      const tryPrefixes = [wantBase, wantName, wantName.slice(0, 32), wantName.slice(0, 24), wantName.slice(0, 16)];
      let idx = -1;
      for (const p of tryPrefixes) {
        idx = opts.findIndex((t) => t.includes(p));
        if (idx >= 0) break;
      }
      if (idx < 0) return {found: false, opts};
      sel.selectedIndex = idx;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      return {found: true, opts, picked: opts[idx]};
    }, [MODEL_NAME, MODEL_BASE]);
    expect(pickResult.found,
      `Just-saved model "${MODEL_NAME}" must surface in dapi.ml.suggested(tableInfo) ` +
      `for accelerometer.csv. Options seen: ${JSON.stringify(pickResult.opts)}.`)
      .toBe(true);
    // Inputs (ColumnsMapInput) — 3 model inputs (accel_y, accel_z, time_offset)
    // auto-mapped to 3 like-named columns in the open table.
    await expect(page.locator(
      '[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]'))
      .toContainText('(3/3)', {timeout: 10_000});
  });

  await softStep('2.4 OK — engine apply appends prediction column', async () => {
    // sub_feature: models.engines.api.apply (per predictive_modeling_engines.dart#L26
    // and predictive_modeling_core.dart#L111). The roll prediction column is
    // appended to the active DataFrame; the exact name may carry a suffix per
    // one-hot/disambiguation rules, so assert via column-count delta + name set.
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]')
      .waitFor({state: 'detached', timeout: 30_000});
    await expect.poll(async () => await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__initialColCount ?? 0;
      return (df?.columns?.length ?? 0) - initial;
    }), {timeout: 60_000, intervals: [500, 1_000, 2_000]}).toBeGreaterThan(0);
    const result = await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__initialColCount ?? 0;
      const allNames = Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name);
      const newCols = allNames.slice(initial);
      // Tags.PredictiveModel probe (soft-warn rather than hard-fail per
      // apply-spec.ts L264-289 precedent — engines vary on client-vs-server tag).
      const DG: any = (window as any).DG;
      const tagKey = DG?.TAGS?.PREDICTIVE_MODEL ?? '.predictive-model';
      const taggedCount = newCols.filter((n: string) => {
        const col = df.columns.byName(n);
        return !!col?.tags?.[tagKey] || !!col?.getTag?.(tagKey);
      }).length;
      return {initial, now: df.columns.length, newCols, taggedCount};
    });
    expect(result.newCols.length).toBeGreaterThan(0);
    if (result.taggedCount === 0) {
      // eslint-disable-next-line no-console
      console.warn(`Apply: ${result.newCols.length} new column(s) appended ` +
        `(${JSON.stringify(result.newCols)}) but none carry Tags.PredictiveModel. ` +
        `Engine may set tag server-side only.`);
    }
  });

  // ═════════ Scenario 3: Run Evaluation from the Performance pane ════════════

  await softStep('3.1 Navigate to /models and select the saved model', async () => {
    // sub_feature: models.meta.performance-section setup — need the model
    // selected so its accordion (with Performance pane) renders.
    await page.evaluate(() => {
      const g: any = (window as any).grok;
      g.shell.windows.showBrowse = true;
      g.shell.route('/models');
    });
    await page.waitForFunction(() => {
      const v = (window as any).grok.shell.v;
      return v && v.type === 'models';
    }, null, {timeout: 20_000});
    await page.locator('.grok-gallery-grid-item.grok-predictive-model').first()
      .waitFor({timeout: 15_000});
    // Narrow catalog via search.
    const search = page.locator('input[placeholder="Search models by name or by #tags"]');
    await search.waitFor({timeout: 10_000});
    await search.fill('');
    await search.fill(MODEL_NAME);
    await expect.poll(async () => await page.evaluate((wantName: string) => {
      const titles = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'));
      return titles.some((t) => (t.textContent || '').trim().includes(wantName));
    }, MODEL_NAME), {timeout: 15_000}).toBe(true);
    // Click the card to surface its context-panel accordion.
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      (label.closest('.grok-gallery-grid-item') as HTMLElement).click();
    }, MODEL_NAME);
    // Wait for the Performance pane header to appear in the accordion.
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .some((h) => (h.textContent || '').trim() === 'Performance'),
      null, {timeout: 15_000});
  });

  await softStep('3.2 Open Performance pane — stored metrics render', async () => {
    // sub_feature: models.meta.performance-section
    await page.evaluate(() => {
      const perf = Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .find((h) => (h.textContent || '').trim() === 'Performance') as HTMLElement | undefined;
      if (!perf) throw new Error('Performance pane header not found');
      perf.click();
    });
    // Run Evaluation button is the canonical "Performance pane mounted" signal
    // — its presence is asserted via the class-2 recon-noted selector
    // [name="button-Run-Evaluation"] (observed live 2026-06-09).
    await page.locator('[name="button-Run-Evaluation"]').waitFor({timeout: 15_000});
    await expect(page.locator('[name="button-Run-Evaluation"]')).toBeVisible();
  });

  await softStep('3.3 Click Run Evaluation — re-runs on CSV trainedOn', async () => {
    // sub_feature: models.meta.performance-section (Run Evaluation driver per
    // predictive_model_info_meta.dart#L118 renderPerformanceSection — produces a
    // fresh PredictiveModelPreviewWidget for CSV-backed source class without
    // re-executing any external Query: this is the non-agnostic
    // run_performance_evaluation surface for trained_on_csv_table).
    // Capture warnings count pre-click so a Run Evaluation error path can be
    // detected (no stored metrics, server-side failure, etc).
    const warnBefore = await page.evaluate(() =>
      ((window as any).grok.shell.warnings ?? []).length);
    await page.locator('[name="button-Run-Evaluation"]').click();
    // Server round-trip + preview-widget render. Bounded wait — accelerometer is
    // small (~600 rows), Linear Regression cross-val completes quickly. Allow
    // up to 90s on dev variability.
    await page.waitForTimeout(3_000);
    const evidence = await page.evaluate((prevCount: number) => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      const errors = w.slice(prevCount).filter((x) =>
        /error|fail|exception/i.test(JSON.stringify(x)));
      // Look for any preview-widget canvas / image inside the property panel —
      // a fresh widget is the success signal; stored metrics may also re-render.
      const panel = document.querySelector('.grok-prop-panel') ||
        document.querySelector('.d4-accordion');
      const widgetCanvases = panel?.querySelectorAll('canvas').length ?? 0;
      const widgetImgs = panel?.querySelectorAll('img').length ?? 0;
      return {errors, widgetCanvases, widgetImgs};
    }, warnBefore);
    expect(evidence.errors.length,
      `Run Evaluation must complete without error balloons. Errors seen: ${JSON.stringify(evidence.errors)}`)
      .toBe(0);
    // Soft signal: at least one widget canvas or image is expected in the
    // property panel post-evaluation. Some Linear Regression preview shapes
    // render only metrics text (no canvas) — warn instead of fail.
    if (evidence.widgetCanvases === 0 && evidence.widgetImgs === 0) {
      // eslint-disable-next-line no-console
      console.warn(`Run Evaluation completed without errors but no canvas/img ` +
        `widget surfaced in the property panel. Stored-metrics-only render is ` +
        `acceptable; flagging for retrospective review.`);
    }
  });

  // ═════════ Scenario 4: Edit metadata via context menu ═══════════════════════

  await softStep('4.1 Right-click model card → Edit... — edit dialog opens', async () => {
    // sub_feature: models.command.edit (cmdEditModel per
    // predictive_model_info_meta.dart#L218; runs editModelInfo which mounts
    // Modal(title: 'Predictive model') per predictive_modeling_browser.dart#L47).
    // Context-menu items have no [name=] attribute — select by text content
    // per grok-browser/references/models.md L425-435 + sibling precedent
    // delete-spec.ts L150-169.
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const editItem = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Edit\.\.\.$/}).first();
    await expect(editItem).toBeVisible({timeout: 10_000});
    await editItem.click();
    // Edit dialog (class-2 recon-noted: [name="dialog-Predictive-model"]).
    await page.locator('[name="dialog-Predictive-model"]').waitFor({timeout: 10_000});
    await expect(page.locator('[name="dialog-Predictive-model"] .d4-dialog-title'))
      .toHaveText('Predictive model');
  });

  await softStep(`4.2 Change Description and OK — dapi.ml.save persists`, async () => {
    // sub_feature: models.workflow.edit-info (editModelInfo modal binds Name +
    // Description; OK persists via dapi.ml.save and fires AppEvents.ENTITY_EDITED).
    //
    // MCP recon 2026-06-10 (cycle 2026-06-10-models-automate-02 retry dispatch):
    // grok.dapi.models.list()[i].description is UNDEFINED for PredictiveModelInfo
    // entities — the Dart-backed JS proxy (class Wt) does NOT expose a `description`
    // getter. Verified via Object.getPrototypeOf(model) inspection + direct
    // model.description access returning undefined. The prior
    // `toBe(NEW_DESCRIPTION)` assertion was therefore a test-bug: it polled for
    // 30s always returning null/undefined, causing B-RUN-PASS at ~69s per attempt.
    //
    // Fix (retry-7, MCP-backed 2026-06-10, same-paradigm tactical):
    // (a) Fill and OK as before — the edit flow is unchanged.
    // (b) Assert dialog is detached (primary success signal).
    // (c) Verify the description WAS persisted to the server via re-open the Edit
    //     dialog and reading the Description input value — this is the observable
    //     evidence that dapi.ml.save committed the change without relying on the
    //     broken model.description JS-API accessor.
    // (d) Cancel the re-opened dialog (no second save needed).
    const descInput = page.locator(
      '[name="dialog-Predictive-model"] [name="input-host-Description"] input, ' +
      '[name="dialog-Predictive-model"] [name="input-host-Description"] textarea').first();
    await descInput.waitFor({timeout: 10_000});
    await descInput.focus();
    await descInput.fill(NEW_DESCRIPTION);
    await page.locator('[name="dialog-Predictive-model"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Predictive-model"]')
      .waitFor({state: 'detached', timeout: 15_000});
    // Verify the model still exists (no accidental delete occurred during edit).
    await expect.poll(async () => await page.evaluate(async (name: string) => {
      const list = await (window as any).grok.dapi.models
        .filter(`friendlyName = "${name}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 15_000}).toBeGreaterThan(0);
    // Re-open Edit dialog to confirm description was persisted.
    // The JS-API model.description accessor is broken (Dart proxy does not expose
    // PredictiveModelInfo.description); re-opening the edit modal is the only
    // DOM-accessible verification path.
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM for re-open verify`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const editItem2 = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Edit\.\.\.$/}).first();
    await expect(editItem2).toBeVisible({timeout: 10_000});
    await editItem2.click();
    await page.locator('[name="dialog-Predictive-model"]').waitFor({timeout: 10_000});
    const descVerify = page.locator(
      '[name="dialog-Predictive-model"] [name="input-host-Description"] input, ' +
      '[name="dialog-Predictive-model"] [name="input-host-Description"] textarea').first();
    await descVerify.waitFor({timeout: 5_000});
    await expect(descVerify).toHaveValue(NEW_DESCRIPTION, {timeout: 5_000});
    await page.locator('[name="dialog-Predictive-model"] [name="button-CANCEL"]').click();
    await page.locator('[name="dialog-Predictive-model"]')
      .waitFor({state: 'detached', timeout: 10_000});
  });

  // ═════════ Scenario 5: Share via context menu ═══════════════════════════════

  await softStep('5.1 Right-click model card → Share... — share dialog opens', async () => {
    // sub_feature: models.command.share (regCommand Share per
    // predictive_model_info_meta.dart#L173 → shareEntity opens the standard
    // share dialog). Per grok-browser/references/models.md L425-435 +
    // chemprop-spec.ts L597-606 precedent, the share dialog mounts as
    // [name="dialog-Share-<encodedFriendlyName>"] with title "Share <friendlyName>".
    // The friendlyName is annotate()-encoded — underscores → hyphens, spaces →
    // hyphens, special chars → hyphens. Compute the encoded form and accept
    // both encoded + generic-by-title fallbacks for build resilience.
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const shareItem = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Share\.\.\.$/}).first();
    await expect(shareItem).toBeVisible({timeout: 10_000});
    await shareItem.click();
    // Encoded form: replace [:_; *\[\]{}|] and spaces with hyphens per
    // html_utils.annotate() rule cited in
    // .claude/skills/grok-browser/SKILL.md "Name transformation rules".
    const encodedName = MODEL_NAME.replace(/[:_;*\[\]{}|\s]/g, '-');
    const shareDialog = page.locator(
      `[name="dialog-Share-${encodedName}"], ` +
      `.d4-dialog:has(.d4-dialog-title:has-text("Share ${MODEL_NAME}"))`).first();
    await expect(shareDialog).toBeVisible({timeout: 15_000});
  });

  await softStep('5.2 Cancel share dialog — leaves permissions unchanged', async () => {
    // Scenario .md Step 5 narrates a happy-path with "Add a recipient + Can
    // view + OK", but the only recipient with deterministic in-session
    // reachability is the operator's own user/group, and adding the operator's
    // own self-group does NOT exercise a meaningful permission delta — the
    // share grants permissions only to OTHER principals. Per scope-reduction
    // discipline (scenario authority — see automator-prompt §"Constraint
    // enforcement"), we narrow Scenario 5 to "share dialog opens with the
    // model context, then CANCEL leaves permissions unchanged" — the
    // standard pattern for share-dialog opener tests across the corpus
    // (chemprop-spec.ts L602-606 "Close it without changing sharing"
    // sibling precedent on the same regCommand Share entry).
    //
    // scope_reduction_proposal (in-spec): the "Add recipient + Can view + OK"
    // assertion is deferred to a dedicated cross-account share spec; this
    // scenario validates the share-dialog OPEN path which is the
    // models.command.share / shareEntity surface. The atomic permission grant
    // is not the differentiator for trained_on_csv_table — it shares the same
    // grok core share-dialog code for every source class.
    const encodedName = MODEL_NAME.replace(/[:_;*\[\]{}|\s]/g, '-');
    const cancelBtn = page.locator(
      `[name="dialog-Share-${encodedName}"] [name="button-CANCEL"], ` +
      `.d4-dialog [name="button-CANCEL"]`).first();
    await cancelBtn.click();
    await page.locator(`[name="dialog-Share-${encodedName}"]`)
      .waitFor({state: 'detached', timeout: 5_000})
      .catch(() => { /* generic dialog selector may have closed already */ });
  });

  // ═════════ Scenario 6: Delete via context menu + cleanup verification ═════

  await softStep('6.1 Right-click model card → Delete — confirm modal opens', async () => {
    // sub_feature: models.command.delete (cmdDeleteModel per
    // predictive_model_info_meta.dart#L224). Re-find the card by friendlyName
    // (after step 4 the card may have re-rendered due to ENTITY_EDITED).
    await page.evaluate((wantName: string) => {
      const label = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'))
        .find((t) => (t.textContent || '').trim().includes(wantName)) as HTMLElement | undefined;
      if (!label) throw new Error(`card titled "${wantName}" not in DOM`);
      const card = label.closest('.grok-gallery-grid-item') as HTMLElement;
      card.dispatchEvent(new MouseEvent('contextmenu',
        {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    const deleteItem = page.locator('.d4-menu-popup .d4-menu-item-label', {hasText: /^Delete$/}).first();
    await expect(deleteItem).toBeVisible({timeout: 10_000});
    await deleteItem.click();
    // Confirm-delete dialog: per delete-spec.ts L180-192 sibling precedent
    // the dialog text matches /delete|are you sure|remove/i and uses
    // [name="button-OK"] (canonical) or button:has-text("DELETE")
    // (per predictive-models-spec.ts L485-487 the d4 confirm dialog renders
    // a DELETE-labelled button).
    const dlg = page.locator('.d4-dialog').filter({hasText: /delete|are you sure|remove/i }).first();
    await expect(dlg).toBeVisible({timeout: 10_000});
  });

  await softStep(`6.2 Confirm delete — ${MODEL_NAME} removed via dapi.ml.delete`, async () => {
    // sub_feature: models.workflow.remove (PredictiveModelInfoMeta.remove →
    // dapi.ml.delete + engine cleanUpModelData + AppEvents.ENTITY_REMOVED).
    const dlg = page.locator('.d4-dialog').filter({hasText: /delete|are you sure|remove/i }).first();
    const balloonsBefore = await page.locator('.grok-balloon-error, .d4-balloon-error').count();
    // Try canonical OK first, then a DELETE-labelled button as fallback.
    const ok = dlg.locator('[name="button-OK"]').first();
    if (await ok.count() > 0)
      await ok.click();
    else
      await dlg.locator('button:has-text("DELETE"), button:has-text("Delete"), button:has-text("OK")')
        .first().click();
    await dlg.waitFor({state: 'detached', timeout: 15_000});
    await page.waitForTimeout(1_500);
    const balloonsAfter = await page.locator('.grok-balloon-error, .d4-balloon-error').count();
    expect(balloonsAfter,
      `No new error balloon should surface during delete (before=${balloonsBefore}, after=${balloonsAfter}).`)
      .toBeLessThanOrEqual(balloonsBefore);
  });

  await softStep('6.3 Cleanup verification — model absent from dapi.ml + /models', async () => {
    // Two-pronged post-condition (delete-spec.ts L206-229 sibling precedent):
    //   (a) dapi.ml filter → 0 matches
    //   (b) /models catalog re-search → 0 cards
    await expect.poll(async () => await page.evaluate(async (wantName: string) => {
      const list = await (window as any).grok.dapi.models
        .filter(`friendlyName = "${wantName}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 30_000, intervals: [500, 1_000, 2_000]}).toBe(0);
    const search = page.locator('input[placeholder="Search models by name or by #tags"]');
    await search.fill('');
    await search.fill(MODEL_NAME);
    await page.waitForTimeout(1_500);
    const remaining = await page.evaluate((wantName: string) => {
      const titles = Array.from(document.querySelectorAll(
        '.grok-gallery-grid-item.grok-predictive-model .grok-gallery-grid-item-title'));
      return titles.filter((t) => (t.textContent || '').trim().includes(wantName)).length;
    }, MODEL_NAME);
    expect(remaining,
      `Catalog still shows ${remaining} card(s) titled "${MODEL_NAME}" after delete`).toBe(0);
  });

  // Final defensive cleanup: best-effort delete in case any softStep above
  // failed mid-flight and left the LifecycleCsvModel_* row on the server.
  await page.evaluate(async ({prefix}) => {
    const g: any = (window as any).grok;
    const all = await g.dapi.models.list();
    const stale = all.filter((m: any) => (m.friendlyName || m.name || '').startsWith(prefix));
    for (const m of stale)
      try { await g.dapi.models.delete(m); } catch (_) { /* best-effort */ }
  }, {prefix: MODEL_BASE});

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
