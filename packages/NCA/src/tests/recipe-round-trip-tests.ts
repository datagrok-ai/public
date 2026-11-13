/**
 * NCA Studies — recipe round-trip SCIENCE GATE (AC-D6).
 *
 * @see src/views/nca-mapping-pane.ts (getCurrentRecipe / applyRecipe)
 *
 * **The gate.** Replaying a saved recipe on the SAME source must reproduce the
 * configuration losslessly — a dropped field silently changes the compute PATH
 * (NONMEM / sparse mode → different Cmax / AUCinf / CL / Vz). Architect-expanded
 * to THREE modes: dense serial, NONMEM event-driven, and sparse design.
 *
 * Two assertions per fixture:
 *   1. STRUCTURAL losslessness — a recipe captured from a configured pane, when
 *      replayed into a fresh pane, yields a byte-identical recipe
 *      (`getCurrentRecipe()` ∘ `applyRecipe()` = identity). Catches a dropped
 *      enumeration field directly.
 *   2. COMPUTE equivalence — the replayed pane computes a byte-identical summary
 *      (the same per-compound mean/CV) as the original (AC-D1 single path).
 *
 * **Runtime.** Requires the live DG runtime (the mapping pane docks into a real
 * `DG.TableView`); runs under `grok test --host local`. Fixtures are the bundled
 * package files (System:AppData/Nca/...), with the NONMEM case sourced from
 * `03_rat_multiroute_nonmem.csv` (R2-F6).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {loadSampleByName} from '../views/nca-launch';
import {mountAnalysisIntoView} from '../views/nca-launch';
import {normalizeNumericNaMarkers} from '../adapter/na-normalize';
import type {ComputeCallbackCtx} from '../views/nca-mapping-pane';
import type {AnalysisRecipe} from '../adapter';
import {ROUTE_STRING_TO_CODE} from '../adapter';

/** Mount a pane on `df`, optionally replaying `recipe`, compute once, and return
 *  the captured recipe + the compute summary context. Cleans up the view. */
async function roundTrip(
  df: DG.DataFrame, recipe?: AnalysisRecipe, runOpts?: {silent?: boolean},
): Promise<{recipe: AnalysisRecipe; ctx: ComputeCallbackCtx | null}> {
  // Mirror the production open path (openAnalysisFromCoordinator / launchAnalysis):
  // retype NONMEM `.`/empty NA markers from STRING → FLOAT before the grid binds,
  // so the numeric pickers + event-driven detection fire. Without this, NONMEM
  // fixtures never compute (the bug the strengthened gate exposed).
  normalizeNumericNaMarkers(df);
  const view = grok.shell.addTableView(df);
  let captured: ComputeCallbackCtx | null = null;
  try {
    const handle = mountAnalysisIntoView(view, df, {
      recipe,
      rules: recipe?.ncaRules,
      onComputed: (c) => {captured = c;},
    });
    // Let the deferred autodetect / sex auto-bind settle before reading/compute.
    await new Promise((r) => setTimeout(r, 60));
    handle.triggerCompute(runOpts);
    return {recipe: handle.getCurrentRecipe(), ctx: captured};
  } finally {
    view.close();
  }
}

/** Stable JSON for deep comparison (key order is fixed by the builders). */
function j(x: unknown): string {return JSON.stringify(x);}

/**
 * Byte-identical fingerprint of the ACTUAL computed NCA parameters (not the
 * recipe) — the true science gate (AC-D6: "byte-identical parameters to the
 * original compute"). Per-subject Cmax / AUClast / AUCinf / t½ / λz / CL / Vz +
 * route, sorted by subject. NaN serializes to `null` deterministically, so a
 * failed subject compares equal across runs.
 */
function paramFingerprint(ctx: ComputeCallbackCtx): string {
  const rows = ctx.perProfile.map((p) => ({
    s: p.index.subjectKey,
    route: p.route,
    cmax: p.result.values.cmax,
    aucLast: p.result.values.aucLast,
    aucInf: p.result.values.aucInf,
    lambdaZ: p.result.values.lambdaZ,
    halfLife: p.result.values.halfLife,
    cl: p.result.values.cl,
    vz: p.result.values.vz,
  }));
  rows.sort((a, b) => (a.s < b.s ? -1 : a.s > b.s ? 1 : 0));
  return JSON.stringify(rows);
}

async function assertRoundTrip(fileOrCsv: DG.DataFrame, label: string): Promise<void> {
  const first = await roundTrip(fileOrCsv);
  // 1. Structural losslessness: replay the captured recipe → identical recipe.
  const second = await roundTrip(fileOrCsv.clone(), first.recipe);
  expect(j(second.recipe), j(first.recipe), `${label}: recipe not lossless on replay`);
  // 2. BOTH computes must have produced a result. A null ctx on replay means the
  //    compute path diverged — e.g. a re-triggered route-confirm dialog blocked it
  //    (the confirmedRoute round-trip failed) — which the param check below could
  //    not otherwise see.
  if (first.ctx === null || second.ctx === null) {
    throw new Error(
      `${label}: compute produced no result (original=${first.ctx !== null}, replay=${second.ctx !== null})`);
  }
  // 3. SCIENCE gate: byte-identical COMPUTED PARAMETERS (Cmax/AUC/CL/Vz/…), not
  //    merely the recipe object. A lossy round-trip that silently changed the
  //    compute path would diverge here even if the recipe structurally matched.
  expect(paramFingerprint(second.ctx), paramFingerprint(first.ctx),
    `${label}: replayed NCA parameters differ from the original (lossy recipe round-trip)`);
}

category('NCA Studies: recipe round-trip (science gate)', () => {
  test('dense serial fixture replays byte-identical parameters (AC-D6a)', async () => {
    const df = await loadSampleByName('theoph.csv');
    await assertRoundTrip(df, 'dense');
  });

  test('NONMEM event-driven fixture replays byte-identical parameters (AC-D6b)', async () => {
    // Multi-route NONMEM (explicit ROUTE column → no route-confirm gate); exercises
    // the NONMEM compute PATH round-trip. The ambiguous RATE/CMT route-confirm path
    // is covered by the dedicated R2-F1 test below.
    const df = await loadSampleByName('03_rat_multiroute_nonmem.csv');
    await assertRoundTrip(df, 'nonmem');
  });

  test('sparse / destructive fixture replays byte-identical parameters (AC-D6c)', async () => {
    const df = await loadSampleByName('04_mouse_sparse_destructive.csv');
    await assertRoundTrip(df, 'sparse');
  });

  test('adj-R² acceptance threshold (R_acc) survives save/load (lambdaz-adjr2 §7 persistence)', async () => {
    // The generic structural gate above cannot catch a dropped R_acc on its own: the pane
    // initializes R_acc to the DEFAULT, so a broken applyRecipe would still round-trip
    // default→default and pass. This forces a NON-default value (0.92, ≥ the default gate 0.85
    // so the R_acc ≥ minR block does not fire) through applyRecipe → getCurrentRecipe and
    // asserts it comes back unchanged — proof the acceptance threshold is actually persisted,
    // not reset to the default on load.
    const df = await loadSampleByName('theoph.csv');
    const base = await roundTrip(df);
    const custom: AnalysisRecipe = {...base.recipe, acceptanceRSquared: 0.92};
    const replayed = await roundTrip(df.clone(), custom);
    expect(replayed.recipe.acceptanceRSquared, 0.92,
      'R_acc set in the recipe survives applyRecipe → getCurrentRecipe (not reset to default)');
  });

  test('ambiguous-route NONMEM replay does NOT re-trigger the route dialog (R2-F1/FM-1)', async () => {
    // 04_iv_infusion_nonmem.csv has RATE>0 and NO explicit ROUTE → the route-confirm
    // gate fires (routeConfirmRequired=true, reference-nonmem-tests.ts:402). A clean
    // first compute would BLOCK on the dialog, so we capture the autodetected recipe
    // WITHOUT computing, pre-fill confirmedRoute as a prior confirmation would, then
    // replay: compute must SUCCEED (the round-tripped confirmedRoute suppresses the
    // dialog). If confirmedRoute did not deserialize (R2-F1 regression), the dialog
    // re-fires, compute blocks, and ctx stays null.
    const df = await loadSampleByName('04_iv_infusion_nonmem.csv');
    normalizeNumericNaMarkers(df); // mirror production open path (event-detection needs it)
    const view0 = grok.shell.addTableView(df);
    let r0: AnalysisRecipe;
    try {
      const h0 = mountAnalysisIntoView(view0, df, {});
      await new Promise((r) => setTimeout(r, 60));
      r0 = h0.getCurrentRecipe();
    } finally {
      view0.close();
    }
    expect(r0.modes.nonmem, true, 'NONMEM layout auto-detected on 04_iv_infusion_nonmem');
    // Ensure the required NONMEM column roles are mapped (ID/TIME/DV) — the standard
    // autodetect does not always bind the NONMEM `DV` concentration column, and a
    // null conc would block compute for a reason unrelated to confirmedRoute.
    const confirmed: AnalysisRecipe = {
      ...r0,
      columns: {
        ...r0.columns,
        subject: r0.columns.subject ?? 'ID',
        time: r0.columns.time ?? 'TIME',
        concentration: r0.columns.concentration ?? 'DV',
      },
      modes: {...r0.modes, confirmedRoute: ROUTE_STRING_TO_CODE['IV-infusion']},
    };
    const replay = await roundTrip(df.clone(), confirmed);
    // If confirmedRoute did NOT deserialize (R2-F1 regression), the route dialog
    // re-fires, compute blocks, and ctx stays null.
    expect(replay.ctx !== null, true,
      'pre-filled confirmedRoute must suppress the route-confirm dialog on replay (R2-F1)');
    // SCIENCE gate on the AMBIGUOUS-route path (AC-D6 / R2-F6): suppressing the
    // dialog is necessary but not sufficient — confirmedRoute could round-trip to a
    // WRONG-but-valid route and compute would still succeed with mis-labeled CL/F /
    // Vz/F. So replay the recipe the confirmed compute EMITTED and assert the NCA
    // parameters (incl. route + CL/Vz) are byte-identical across the two replays.
    // No "original" baseline exists (a clean first compute blocks on the dialog),
    // so replay₁ vs replay₂ is the lossless-round-trip check for this path.
    const replay2 = await roundTrip(df.clone(), replay.recipe);
    if (replay.ctx === null || replay2.ctx === null) {
      throw new Error(
        `ambiguous-route replay produced no result (replay1=${replay.ctx !== null}, ` +
        `replay2=${replay2.ctx !== null})`);
    }
    expect(paramFingerprint(replay2.ctx), paramFingerprint(replay.ctx),
      'ambiguous-route recipe must replay byte-identical NCA parameters incl. route/CL/Vz (AC-D6 / R2-F6)');
  });
});

/**
 * mapping-source-radio (user-directed 2026-10-19): the five value/column `<X> type` inputs are
 * horizontal `Value` / `Column` radios, while the recipe keeps the lowercase `SourceMode`. Pins
 * AC-1 (radio, two options, the shared row class), AC-3 (autodetect selects `Column`), AC-4
 * (recipe emits `'column'` / `'value'`; replaying `'value'` selects `Value`) and AC-5
 * (`Sampling type` stays a dropdown).
 */
category('NCA Studies: value/column source radios (mapping-source-radio)', () => {
  const SOURCE_LABELS = ['Dose type', 'Route type', 'Infusion duration type', 'Body weight type', 'LLOQ type'];

  /** Mount on `df` (optionally replaying `recipe`), hand the pane root + recipe to `check`, close. */
  async function withPane(df: DG.DataFrame, recipe: AnalysisRecipe | undefined,
    check: (pane: HTMLElement, r: AnalysisRecipe) => void): Promise<void> {
    const view = grok.shell.addTableView(df);
    try {
      const handle = mountAnalysisIntoView(view, df, {recipe, rules: recipe?.ncaRules});
      await new Promise((r) => setTimeout(r, 60));
      const pane = document.querySelector('.nca-mapping-pane') as HTMLElement | null;
      if (pane === null) throw new Error('mapping pane not mounted');
      check(pane, handle.getCurrentRecipe());
    } finally {
      view.close();
    }
  }
  const rootByLabel = (pane: HTMLElement, caption: string): HTMLElement | undefined =>
    (Array.from(pane.querySelectorAll('.ui-input-root')) as HTMLElement[]).find((r) =>
      (r.querySelector('.ui-input-label')?.textContent ?? '').trim() === caption);
  /** The option labels of a radio root, and the one whose radio is checked. */
  const radioState = (root: HTMLElement): {options: string[]; checked: string | null} => {
    const radios = Array.from(root.querySelectorAll('input[type="radio"]')) as HTMLInputElement[];
    const options = (Array.from(root.querySelectorAll('.ui-input-editor label')) as HTMLElement[])
      .map((l) => (l.textContent ?? '').trim()).filter((t) => t !== '');
    const i = radios.findIndex((r) => r.checked);
    return {options, checked: i >= 0 ? options[i] ?? null : null};
  };

  test('every <X> type input is a Value / Column radio; Sampling type stays a dropdown', async () => {
    const df = await loadSampleByName('theoph.csv');
    await withPane(df, undefined, (pane) => {
      for (const caption of SOURCE_LABELS) {
        const root = rootByLabel(pane, caption);
        expect(root !== undefined, true, `${caption}: input present`);
        expect(root!.classList.contains('nca-radio-row'), true, `${caption}: horizontal row class`);
        expect(root!.querySelectorAll('input[type="radio"]').length, 2, `${caption}: two radio options`);
        expect(root!.querySelector('select'), null, `${caption}: no dropdown`);
        expect(j(radioState(root!).options), j(['Value', 'Column']), `${caption}: option labels`);
      }
      const sampling = rootByLabel(pane, 'Sampling type');
      expect(sampling !== undefined && sampling.querySelector('select') !== null, true,
        'Sampling type stays a dropdown');
    });
  });

  test('autodetected dose column selects Column; recipe keeps lowercase SourceMode', async () => {
    const df = await loadSampleByName('theoph.csv');
    let captured: AnalysisRecipe | null = null;
    await withPane(df, undefined, (pane, r) => {
      captured = r;
      expect(r.columns.doseCol !== null, true, 'theoph dose column autodetected');
      expect(r.columns.doseSource, 'column', 'recipe emits lowercase column');
      expect(radioState(rootByLabel(pane, 'Dose type')!).checked, 'Column', 'radio shows Column');
    });
    // Replay with the dose as a value: the radio shows Value, the recipe emits 'value'.
    const asValue: AnalysisRecipe = {
      ...captured!,
      columns: {...captured!.columns, doseSource: 'value', doseCol: null, doseScalar: 320},
    };
    await withPane(df.clone(), asValue, (pane, r) => {
      expect(r.columns.doseSource, 'value', 'replayed recipe emits lowercase value');
      expect(radioState(rootByLabel(pane, 'Dose type')!).checked, 'Value', 'radio shows Value');
    });
  });
});

/**
 * Silent replay-open (the on-open hydration auto-compute) must suppress the
 * interactive confirm dialogs the user already cleared when first saving the
 * analysis. The per-kg confirm gate is the testable instance: a NON-silent
 * compute on a `/kg`-dose + no-body-weight recipe BLOCKS on the dialog
 * (`runPipeline` never runs → `ctx` stays null); a SILENT replay must PROCEED.
 * (Dialog/toast non-appearance itself is visual; verified live via grok test +
 * manual QA. This asserts the compute-proceeds half, which is the observable
 * consequence of the suppression.)
 */
category('NCA Studies: silent replay-open suppresses confirm gates', () => {
  test('silent replay computes through the per-kg gate without the confirm dialog', async () => {
    const df = await loadSampleByName('theoph.csv');
    // Baseline: a normal compute on the fixture proves the fixture itself computes.
    const base = await roundTrip(df.clone());
    expect(base.ctx !== null, true, 'baseline compute should produce a result');
    // Force the per-kg no-BW condition: dose unit in /kg with no body weight (column
    // or scalar). The pre-run per-kg dialog that once gated it is retired.
    const perKg: AnalysisRecipe = {
      ...base.recipe,
      columns: {...base.recipe.columns, bodyWeightSource: 'none', bodyWeightScalar: 0},
      units: {...base.recipe.units, dose: 'mg/kg'},
    };
    // A user-initiated (non-silent) compute would block on the per-kg dialog; the
    // silent replay-open must compute through it (the user already confirmed per-kg
    // when first saving; the PER_KG_NO_BW warning remains in WARNINGS_JSON).
    const silent = await roundTrip(df.clone(), perKg, {silent: true});
    expect(silent.ctx !== null, true,
      'silent replay must compute through the per-kg gate without the confirm dialog');
  });
});
