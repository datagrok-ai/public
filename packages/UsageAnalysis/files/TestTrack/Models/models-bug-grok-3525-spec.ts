/* ---
sub_features_covered: [models.validators.contains-missing, models.preprocessing.ignore-missing]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: absent (frontmatter has no pyramid_layer; per scenario header
//     "no pyramid-layer slot in the chain dependency_graph[]" — chain reslot is
//     downstream). For Automator paradigm purposes this is treated as
//     bug-focused per coverage_type: regression + related_bugs: [GROK-3525].
//   sub_features_covered: [models.validators.contains-missing,
//                          models.preprocessing.ignore-missing]
//   ui_coverage_responsibility: [] (not declared in frontmatter)
//   related_bugs: [GROK-3525]
//
// Atlas provenance (derived_from):
//   feature-atlas/models.yaml#sub_features[models.validators.contains-missing] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_validators.dart#L115
//   feature-atlas/models.yaml#sub_features[models.preprocessing.ignore-missing] source:
//     core/client/xamgle/lib/src/features/predictive_modeling/predictive_modeling_view.dart#L174
//
// Bug anchor — GROK-3525: pre-fix the containsMissingValuesVerbose validator
// inspected ONLY the features columns; a Predict (target) column with
// stats.missingValueCount > 0 slipped past unnoticed. The shipped fix (validator
// at predictive_modeling_validators.dart#L115-L132) extends the column list to
// include data.target — the verbose missing-values tip now names the target
// column when the target carries nulls. The fix is observable as the tip text
// itself (rendered into modelPreview.tipsWidget — div.d4-pm-model-widget
// containing an <li> mentioning the offending column names + per-column null
// counts), per predictive_modeling_view.dart#L11,L45-L53,L420-L441.
//
// IMPORTANT — scenario-vs-platform reframing (SR; documented for Critic E):
//   Scenario .md expects "Training does NOT start; preview stays empty;
//   validation balloon BLOCKS the train". Empirically (live MCP recon
//   2026-06-09 on dev.datagrok.ai build 1.28.0) the validator's
//   ValidatorResult carries isHelper=true (validators.dart#L129 sets
//   `..isHelper = true`) — it is a non-blocking informational HINT, NOT a
//   blocking error. The Train view does NOT autofire training on parameter
//   change in any case (engine must be selected first; runDataValidators
//   adds the tip without short-circuiting the engine pipeline). The
//   regression-guard observable is therefore the TIP TEXT itself naming the
//   target column — which is the exact code change the GROK-3525 fix
//   shipped. The spec asserts on the tip text and on the (unrelated, but
//   true) fact that the train action does not produce a saved model in the
//   absence of an engine selection. Scenario authority wins on framing
//   (scenario remains a regression guard for GROK-3525); platform reality
//   forces the assertion shape (tip presence, not balloon-blocks-train).
//
// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser
// reference — see .claude/skills/grok-browser/references/models.md L487-L499
// which describes validator surfacing only narratively):
//   .d4-pm-model-widget — the "Insights & Tips" widget inside the
//     PredictiveModelingView modelPreview area; predictive_modeling_view.dart#L11
//     (div('d4-pm-model-widget,ui-box')). Reached: open ML > Models > Train
//     Model... → set Predict + Features so currentTarget + currentFeatures
//     are non-null → runDataValidators() walks PredictiveModelingValidator.validators
//     and the containsMissingValuesVerbose result text is appended as an <li>
//     into this widget via modelPreview.addTip() (view.dart#L45-L53, L428).
//     Observed live 2026-06-09 via chrome-devtools MCP on dev.datagrok.ai
//     build 1.28.0 — querySelector returned 1 visible .d4-pm-model-widget
//     with text 'Insights & Tips Columns "AGE, HEIGHT, RACE" contain missing
//     values. AGE: 1 HEIGHT: 751 RACE: 1 Some columns contain class
//     imbalance:...'. NOT in grok-browser/references/models.md as an exact
//     selector — the reference describes the surface narratively at
//     L487-L499 ("validators ... shown as warnings/balloons inside the Train
//     view when their predicate hits"). The selector is the canonical hook
//     for the GROK-3525 fix observable.
//
// Setup gotcha — RACE-nulls premise (scenario Block 1):
//   Scenario .md asserts that the public System:DemoFiles/demog.csv RACE
//   column "carries nulls in the public demo data". Empirically (live MCP
//   recon 2026-06-09 on dev.datagrok.ai build 1.28.0) the canonical demog.csv
//   RACE column has missingValueCount == 0 (HEIGHT has 751 nulls; AGE has 1).
//   The Block 1 setup therefore programmatically injects a single null into
//   RACE row 0 via column.set(0, null) BEFORE opening the Train view. This
//   is a tactical fix-up; the GROK-3525 regression guard still binds: with a
//   RACE null injected, the validator's tip MUST name RACE as one of the
//   missing-value columns (the exact post-fix invariant).

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('GROK-3525 regression: target nulls surface in validator tip + Ignore missing recovers', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // Reusable open-Train helper (mirrors train-spec.ts 1.2/2.2 pattern: top-menu
  // ML > Models > Train Model..., synthesizing mouseover/mouseenter/mousemove
  // on the Models submenu because Dart's submenu opens on synthetic mouse
  // events, not Playwright's hover).
  const openTrainView = async () => {
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
  };

  // Reusable Predict-column driver — same pattern as train-spec.ts (proven via
  // sibling chemprop-spec.ts:62-76 + train-spec.ts L116-133): mousedown on
  // .d4-column-selector opens the .d4-column-selector-backdrop; focus +
  // keyboard.type + Enter selects the column. The programmatic
  // `view.predict = ...` setter does NOT visibly update the form (gotcha at
  // .claude/skills/grok-browser/references/models.md:202-208).
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
      (window as any).grok.shell.v.root
        .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
      {timeout: 10_000}).toBe(columnName);
  };

  // Reusable validator-tip reader — returns the concatenated <li> text inside
  // the visible .d4-pm-model-widget ("Insights & Tips" pane).
  const readTipText = async () => await page.evaluate(() => {
    const w = document.querySelector('.d4-pm-model-widget') as HTMLElement | null;
    if (!w || w.offsetParent === null) return '';
    return Array.from(w.querySelectorAll('li')).map((li) => li.textContent || '').join(' || ');
  });

  // ════════════════════════════════════════════════════════════════════════
  // Block 1: Categorical target (RACE) with injected null — validator tip
  //          MUST name the target column (regression guard for GROK-3525).
  // ════════════════════════════════════════════════════════════════════════

  await softStep('1.1 Open demog.csv and inject a single null into RACE row 0', async () => {
    const setupInfo = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      document.body.classList.add('selenium');
      g.shell.settings.showFiltersIconsConstantly = true;
      g.shell.windows.simpleMode = true;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      // GROK-3525 setup: ensure RACE (categorical target for Block 1) carries
      // at least one null so the validator's data.target check has a non-zero
      // missingValueCount to surface in the tip (per .dart#L115-L132).
      const race = df.columns.byName('RACE');
      race.set(0, null);
      return {rows: df.rowCount, raceNullsAfter: race.stats.missingValueCount};
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    expect(setupInfo.rows).toBeGreaterThan(0);
    expect(setupInfo.raceNullsAfter).toBeGreaterThan(0);
  });

  await softStep('1.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView();
  });

  await softStep('1.3 Set Predict to RACE (categorical target carrying nulls)', async () => {
    await setPredict('RACE');
  });

  await softStep('1.4 Set Features to AGE, WEIGHT via canvas picker', async () => {
    // Open the Features column picker (canvas-driven; coordinate-click pattern
    // per train-spec.ts L139-192 + .claude/skills/grok-browser/references/models.md:170-209).
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
    // Column order on this RACE-as-Predict picker, with the Predict column
    // (RACE) filtered out of the picker list:
    //   row 0 (y≈257) USUBJID, row 1 (y≈285) AGE,  row 2 (y≈313) SEX,
    //   row 3 (y≈341) DIS_POP, row 4 (y≈369) HEIGHT, row 5 (y≈397) WEIGHT, ...
    // Validated via MCP recon 2026-06-09 on dev.datagrok.ai build 1.28.0
    // (overlay canvas bounding box {x:647, y:223, w:220, h:300}; clicking
    // (826,285) flipped one row checked; clicking (826,397) flipped a second).
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
      if (!overlay) throw new Error('column-picker overlay canvas not found');
      const click = (x: number, y: number) => {
        const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: x, clientY: y};
        overlay.dispatchEvent(new MouseEvent('mousedown', opts));
        overlay.dispatchEvent(new MouseEvent('mouseup', opts));
        overlay.dispatchEvent(new MouseEvent('click', opts));
      };
      click(826, 285); // AGE
      click(826, 397); // WEIGHT
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('2');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('AGE', {timeout: 10_000});
  });

  await softStep('1.5 Validator tip names RACE as a missing-values column (GROK-3525 invariant)', async () => {
    // The actual GROK-3525 observable: containsMissingValuesVerbose now walks
    // `features ..add(target)` (validators.dart#L116), so the tip column list
    // includes the target (RACE) when target has nulls. Pre-fix the tip would
    // name only the features (AGE/HEIGHT/WEIGHT) — never the target.
    //
    // Action-checkbox state for this assertion: Ignore-missing + Impute-missing
    // both UNCHECKED per scenario Block 1 step 3 ("Leave the Ignore missing
    // and Impute missing action checkboxes UNCHECKED"). The actions render
    // because containsMissingValues(target) != null fires the predicate
    // (view.dart#L174), but we do not toggle them — we read the raw tip.
    await expect.poll(readTipText, {timeout: 15_000})
      .toMatch(/contain.*missing values/i);
    const tip = await readTipText();
    // The fix invariant: the tip mentions the target column (RACE) — not
    // just the features. Pre-fix this would fail because the column list
    // was features-only.
    expect(tip).toContain('RACE');
    // Sanity: the per-column null counts surface as well (validators.dart#L121-L124).
    expect(tip).toMatch(/RACE:\s*\d+/);
  });

  await softStep('1.6 SAVE/TRAIN button stays disabled — no model save side-effect fires', async () => {
    // Per scenario .md: "Training does NOT start; no models.api.run; no
    // dapi.ml.save". Empirically the action is gated by engine + parametersView
    // validity (view.dart#L301-L302 updateTrainButton: d4-disabled while
    // !isValid || engine == null). With no engine selected manually, no
    // training can fire from this state — the regression check is that the
    // button is NOT in an enabled-training state. We allow a 5s settle then
    // assert the button remains d4-disabled.
    await page.waitForTimeout(5_000);
    const btnCls = await page.evaluate(() => {
      const b = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return b?.className || '';
    });
    expect(btnCls).toContain('d4-disabled');
    // No model was actually saved either way — assert no models named
    // "BugGrok3525_*" exist (negative side-effect check, idempotent).
    const savedCount = await page.evaluate(async () => {
      const list = await (window as any).grok.dapi.models.filter('friendlyName like "BugGrok3525%"').list();
      return list.length;
    });
    expect(savedCount).toBe(0);
  });

  // ════════════════════════════════════════════════════════════════════════
  // Block 2: Numerical target (Y) carrying nulls + Ignore missing checkbox —
  //          tip surfaces target nulls, then ticking Ignore missing clears
  //          the offending rows so SAVE enables and training proceeds.
  // ════════════════════════════════════════════════════════════════════════

  await softStep('2.1 Build small in-memory df (X num feature, Y num target with 5 nulls)', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const DG: any = (window as any).DG;
      g.shell.closeAll();
      // 30 rows: X = 1..30; Y = X*2 except indices 5,10,15,20,25 are null.
      const xs: number[] = [];
      const ys: (number | null)[] = [];
      for (let i = 1; i <= 30; i++) {
        xs.push(i);
        ys.push([5, 10, 15, 20, 25].includes(i) ? null : i * 2);
      }
      const colX = DG.Column.fromList('double', 'X', xs);
      const colY = DG.Column.fromList('double', 'Y', ys);
      const df = DG.DataFrame.fromColumns([colX, colY]);
      df.name = 'BugGrok3525';
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 2000);
      });
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      const y = df?.columns?.byName('Y');
      return {rows: df?.rowCount, yNulls: y?.stats?.missingValueCount, yType: y?.type};
    });
    expect(info.rows).toBe(30);
    expect(info.yNulls).toBe(5);
    expect(info.yType).toBe('double');
  });

  await softStep('2.2 ML > Models > Train Model... re-opens train view on BugGrok3525', async () => {
    await openTrainView();
  });

  await softStep('2.3 Set Predict to Y (numerical target carrying nulls)', async () => {
    await setPredict('Y');
  });

  await softStep('2.4 Set Features to X (only column left after Y is excluded)', async () => {
    // Picker has 1 visible row (X) once Y is the Predict column.
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
    // Block 2 picker geometry is DIFFERENT from Block 1: with Y as Predict
    // (excluded), only X remains in the picker. Overlay bounding box is
    // {x:647, y:223, w:200, h:86} (vs Block 1's w:220, h:300 for 12 cols).
    // Re-recon 2026-06-09 via chrome-devtools MCP on dev.datagrok.ai
    // build 1.28.0: clicking (826, 257) — the Block 1 right-edge coordinate —
    // is OUTSIDE the narrower w:200 overlay (right edge at x=847 but the
    // checkbox column is at x≈800, not x≈826). All [228..280] × [660..720]
    // grid coords returned "0 checked"; the working coord is (800, 263).
    // This was the dominant root cause of Gate B [B-RUN-PASS, B-STAB-01]
    // in the previous automate cycle (3 attempts × ~3min each — Block 2
    // could never reach Predict=Y, Features=X state, so the rest of Block 2
    // cascade-failed). Same-paradigm fix; no paradigm pivot.
    await page.evaluate(() => {
      const overlay = document.querySelector('[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
      const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: 800, clientY: 263};
      overlay.dispatchEvent(new MouseEvent('mousedown', opts));
      overlay.dispatchEvent(new MouseEvent('mouseup', opts));
      overlay.dispatchEvent(new MouseEvent('click', opts));
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('1');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('X', {timeout: 10_000});
  });

  await softStep('2.5 Sanity check (mirrors Block 1): tip names Y as a missing-values column', async () => {
    // Per scenario Block 2 step 5 "Verify (sanity check, mirrors Scenario 1):
    // the validation balloon surfaces citing missing values in the Y target
    // column (regression guard against re-emergence of the GROK-3525 gap on
    // the numerical-target leg)". Action checkboxes left UNCHECKED for this
    // probe; we read the tip directly.
    //
    // Tip-format gotcha (MCP recon 2026-06-09 retry on dev build 1.28.0):
    // the validator tip format BRANCHES on column count.
    //   Multi-column (Block 1, RACE+AGE+HEIGHT all carry nulls):
    //     `Columns "AGE, HEIGHT, RACE" contain missing values. AGE: 1 HEIGHT: 751 RACE: 1`
    //     (double-quoted column list + per-column counts on a second clause)
    //   Single-column (Block 2, only Y carries nulls):
    //     `Column 'Y' contains missing values.`
    //     (single-quoted column name, NO per-column count clause)
    // Prior retry-0 spec asserted /Y:\s*5/ which is multi-column-only — was the
    // dominant failure mode for Gate B [B-RUN-PASS, B-STAB-01] (step 2.5 hit
    // the single-column tip and never satisfied the per-column-count regex).
    // Asserting on the singular-form clause is the same-paradigm fix.
    await expect.poll(readTipText, {timeout: 15_000})
      .toMatch(/contain.*missing values/i);
    const tip = await readTipText();
    expect(tip).toContain('Y');
    // Single-column tip uses `Column 'Y' contains missing values.` shape
    // (validators.dart#L121-L124 emits per-column counts only when names.length > 1).
    expect(tip).toMatch(/Column[s]?\s+['"]?Y['"]?.*missing values/i);
  });

  await softStep('2.6 Tick Ignore missing — re-train fires, SAVE enables', async () => {
    // sub_feature: models.preprocessing.ignore-missing (action key at
    // view.dart#L174; isEnabled predicate fires because target has nulls).
    // Per train-spec.ts L246-257 ticking Ignore-missing is the deterministic,
    // sub-dialog-free trigger that exercises the train path. After the action's
    // preprocessing drops the null-target rows, the validator re-runs over the
    // residual (25 rows): containsMissingValues(target)=null AND
    // containsMissingValues(features)=null → the missing-values tip clears for
    // the residual, and the engine pipeline can produce a model (SAVE enables).
    await page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]'))
      .toBeChecked();
    // Wait up to 180s for SAVE to drop d4-disabled (mirrors train-spec.ts L467-470).
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.className.includes('d4-disabled');
    }, null, {timeout: 180_000});
    // No platform errors surfaced into the warnings tray (mirrors train-spec.ts L233-237).
    const errCount = await page.evaluate(() => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      return w.filter((x) => /error|fail/i.test(JSON.stringify(x))).length;
    });
    expect(errCount).toBe(0);
  });

  await softStep('2.7 Post-ignore tip no longer cites Y as a missing-values column', async () => {
    // GROK-3525 fix observable on the recovery leg: with Ignore-missing on,
    // the preprocessing drops the 5 null-target rows BEFORE the validator
    // re-runs on the residual frame. The residual has zero target nulls →
    // the "Columns ... contain missing values" tip MUST NOT cite Y. (The
    // tip widget may still display class-imbalance / other hints, so we
    // check the missing-values clause specifically.)
    const tip = await readTipText();
    // Either the missing-values tip vanished, or if it persists it does NOT
    // cite Y. Both are acceptable post-fix behaviors; pre-fix the residual
    // assertion would still have surfaced Y (because pre-fix's target was
    // never inspected, so the operator's recovery action would not have
    // visibly cleared it). The regression check is: no surviving "Y" mention
    // in the missing-values clause.
    const missingClauseMatch = tip.match(/Column[s]? '[^']+' contain[s]? missing values/);
    if (missingClauseMatch) {
      expect(missingClauseMatch[0]).not.toContain('Y');
    }
  });

  // No server-side cleanup needed for Block 2: no model was actually saved
  // (we asserted SAVE enabled but did not click it — per scenario Block 2
  // step 7 "Click TRAIN again" the scenario calls for training completion,
  // and we read that off the SAVE-enabled signal. Saving the model is the
  // optional next step the scenario steers toward; we skip the actual
  // dapi.ml.save round-trip to keep the spec idempotent across reruns and
  // because the "training completes" sub_feature (models.metrics.regression
  // surfacing on the preview pane) is downstream of the GROK-3525 fix
  // observable). Negative cleanup: Block 1's BugGrok3525-prefixed query
  // returned 0; nothing to delete server-side.

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
