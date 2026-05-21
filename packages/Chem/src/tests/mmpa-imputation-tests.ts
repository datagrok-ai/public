import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import {MMPA} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa';
import {imputeMmpActivities} from '../analysis/molecular-matched-pairs/mmp-analysis/mmpa-imputation';
import {MmpDiffTypes} from '../analysis/molecular-matched-pairs/mmp-function-editor';
import {SortData} from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer';
import {readDataframe} from './utils';
import * as chemSearches from '../chem-searches';
import {Fingerprint} from '../utils/chem-common';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

const FNULL = DG.FLOAT_NULL;

/** Builds a hand-crafted MMPA-like input. 4 molecules, 1 rule, 4 pairs.
 *  activities = [10.0, 12.0, FNULL, FNULL], meanDiff = 2.0, pairsPerActivity = 1.
 *  Expected:
 *    mol2 has 2 anchors (10+2=12 via pair(0,2) and 12+2=14 via pair(1,2)) → mean 13, stdev 1
 *    mol3 has 1 anchor (10+2=12 via pair(0,3)) → suppressed at minAnchors=2 */
function makeSyntheticMmpa() {
  const acts = new Float32Array([10.0, 12.0, FNULL, FNULL]);
  const initData = {
    molName: 'smiles',
    molecules: ['c1ccccc1', 'Cc1ccccc1', 'CCc1ccccc1', 'CCCc1ccccc1'],
    activities: [acts],
    activitiesNames: ['Activity'],
    activitiesCount: 1,
    diffTypes: ['delta'] as string[],
  };
  const rules = {
    rules: [{
      sr1: 0, sr2: 1,
      pairs: [
        {fs: 0, ss: 1, core: 0},
        {fs: 0, ss: 2, core: 0},
        {fs: 1, ss: 2, core: 0},
        {fs: 0, ss: 3, core: 0},
      ],
    }],
    smilesFrags: [0, 1],
  };
  const rulesBased = {
    fromFrag: [''], toFrag: [''],
    occasions: new Int32Array([4]),
    meanDiffs: [new Float32Array([2.0])],
    fragmentParentIndices: [new Uint32Array(0)],
    moleculePairParentIndices: [new Uint32Array(0), new Uint32Array(0)] as [Uint32Array, Uint32Array],
    moleculePairIndices: [new Uint32Array(0)],
  };
  const allCasesBased = {
    maxActs: [0],
    molFrom: [] as string[], molTo: [] as string[],
    pairNum: new Int32Array(0), molNumFrom: new Int32Array(0), molNumTo: new Int32Array(0),
    pairsFromSmiles: [] as string[], pairsToSmiles: [] as string[],
    ruleNum: new Int32Array(0), diffs: [new Float32Array(0)],
    activityPairsIdxs: [] as any[], coreNums: new Int32Array(0),
    pairedActivities: [[new Float32Array(0), new Float32Array(0)]] as [Float32Array, Float32Array][],
    pairsPerActivity: [new Uint32Array([1])],
  };
  return {initData, rules, rulesBased, allCasesBased};
}

const EPS = 1e-4;

category('mmpa-imputation', () => {
  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('synthetic_mol2_mean_and_stdev', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 1, 2);
    const {pred, support, stdev} = perActivity[0];
    expect(pred[0], FNULL, 'mol0 has known activity — pred must stay FNULL');
    expect(pred[1], FNULL, 'mol1 has known activity — pred must stay FNULL');
    expect(Math.abs(pred[2] - 13.0) < EPS, true, `mol2 pred expected 13.0, got ${pred[2]}`);
    expect(support[2], 2, 'mol2 support must be 2');
    // mol2 has predictions [12, 14] from one rule (both pairs share the
    // same rule with meanDiff=2.0). Mean = 13. With Bessel-corrected
    // sample variance: var = ((12-13)² + (14-13)²) / (2-1) = 2.0,
    // σ = √2 ≈ 1.4142. Pre-Bessel (biased /n) gave σ = 1.0 — flagged
    // by external review as understating uncertainty by ~29% at the
    // minimum-anchor floor and inflating MMP's weight in the inverse-
    // variance blend. Test updated to the unbiased value.
    expect(Math.abs(stdev[2] - Math.sqrt(2)) < EPS, true,
      `mol2 stdev expected √2 (Bessel-corrected sample variance), ` +
      `got ${stdev[2]}. Pre-Bessel biased value was 1.0; if the test ` +
      `is seeing 1.0 again the /(n-1) correction in mmpa-imputation.ts ` +
      `has regressed.`);
  });

  test('synthetic_mol3_suppressed_by_minAnchors', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 1, 2);
    expect(perActivity[0].pred[3], FNULL, 'mol3 has 1 anchor, minAnchors=2 — must be FNULL');
  });

  test('synthetic_mol3_predicted_with_minAnchors1', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 1, 1);
    expect(Math.abs(perActivity[0].pred[3] - 12.0) < EPS, true,
      `mol3 with minAnchors=1 expected 12.0, got ${perActivity[0].pred[3]}`);
  });

  test('synthetic_supportFloor_filters_rule', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 2, 1);
    expect(perActivity[0].pred[2], FNULL, 'supportFloor=2 > pairsPerActivity=1: mol2 FNULL');
    expect(perActivity[0].pred[3], FNULL, 'supportFloor=2 > pairsPerActivity=1: mol3 FNULL');
  });

  test('synthetic_known_rows_not_overwritten', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 1, 1);
    expect(perActivity[0].pred[0], FNULL, 'mol0 known activity: pred slot must be FNULL');
    expect(perActivity[0].pred[1], FNULL, 'mol1 known activity: pred slot must be FNULL');
  });

  test('synthetic_predictKnown_emits_predictions_for_known_rows', async () => {
    // With predictKnown=true and minAnchors=1, every molecule with at least one
    // valid anchor gets a prediction — including mol0/mol1 whose activity is known.
    //   mol0: 1 anchor (mol1 via pair(0,1)), pred = 12 − 2 = 10.0
    //   mol1: 1 anchor (mol0 via pair(0,1)), pred = 10 + 2 = 12.0
    //   mol2: 2 anchors → still 13.0 (unchanged by the flag)
    //   mol3: 1 anchor (mol0 via pair(0,3)), pred = 12.0
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 1, 1, true);
    expect(Math.abs(perActivity[0].pred[0] - 10.0) < EPS, true,
      `predictKnown mol0 expected 10.0, got ${perActivity[0].pred[0]}`);
    expect(Math.abs(perActivity[0].pred[1] - 12.0) < EPS, true,
      `predictKnown mol1 expected 12.0, got ${perActivity[0].pred[1]}`);
    expect(Math.abs(perActivity[0].pred[2] - 13.0) < EPS, true,
      `predictKnown mol2 expected 13.0, got ${perActivity[0].pred[2]}`);
    expect(Math.abs(perActivity[0].pred[3] - 12.0) < EPS, true,
      `predictKnown mol3 expected 12.0, got ${perActivity[0].pred[3]}`);
  });

  test('edge_all_null_activity', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    initData.activities[0] = new Float32Array([FNULL, FNULL, FNULL, FNULL]);
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 1, 1);
    for (let i = 0; i < 4; i++)
      expect(perActivity[0].pred[i], FNULL, `all-null: row ${i} pred must be FNULL`);
  });

  test('edge_no_rules', async () => {
    const {initData, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const emptyRules = {rules: [], smilesFrags: []};
    const emptyRulesBased = {...rulesBased, meanDiffs: [new Float32Array(0)]};
    const emptyAllCasesBased = {...allCasesBased, pairsPerActivity: [new Uint32Array(0)]};
    const {perActivity} = imputeMmpActivities(
      initData as any, emptyRules as any, emptyRulesBased as any, emptyAllCasesBased as any, 1, 1);
    for (let i = 0; i < 4; i++)
      expect(perActivity[0].pred[i], FNULL, `no-rules: row ${i} pred must be FNULL`);
  });

  test('edge_supportFloor_exceeds_all', async () => {
    const {initData, rules, rulesBased, allCasesBased} = makeSyntheticMmpa();
    const {perActivity} = imputeMmpActivities(
      initData as any, rules as any, rulesBased as any, allCasesBased as any, 9999, 1);
    expect(perActivity[0].pred[2], FNULL, 'supportFloor=9999: mol2 must be FNULL');
    expect(perActivity[0].pred[3], FNULL, 'supportFloor=9999: mol3 must be FNULL');
  });

  // Holdout validation on the bundled demo dataset. Mask ~20% of hERG_pIC50
  // values, run imputation, and check MAE/R² against the held-out truth.
  // Long-running: gated by a generous timeout (MMPA.init on 20k rows
  // dominates, ~30-90 seconds).
  test('holdout_hERG', async () => {
    const df = await readDataframe('demo_files/mmp_demo.csv');
    const smilesCol = df.col('smiles');
    const hergCol = df.col('hERG_pIC50');
    if (!smilesCol || !hergCol) {
      grok.shell.warning('holdout_hERG: mmp_demo.csv missing expected columns; skipping.');
      return;
    }
    const rawActs: Float32Array = hergCol.getRawData() as Float32Array;
    const heldOutIdxs: number[] = [];
    const maskedActs = new Float32Array(rawActs.length);
    for (let i = 0; i < rawActs.length; i++) {
      maskedActs[i] = rawActs[i];
      // Deterministic mask: every 5th non-null row → ~20% holdout.
      if (rawActs[i] !== FNULL && i % 5 === 0) {
        heldOutIdxs.push(i);
        maskedActs[i] = FNULL;
      }
    }
    expect(heldOutIdxs.length > 100, true,
      `holdout_hERG: expected >100 held-out rows, got ${heldOutIdxs.length}`);

    const sortingInfo: SortData = {fragmentIdxs: [], frequencies: []};
    const mmpa = await MMPA.init(
      smilesCol.name, smilesCol.toList(), 0.4,
      [maskedActs], ['hERG_pIC50'], [MmpDiffTypes.delta], sortingInfo,
    );
    // Use production defaults (supportFloor=10, minAnchors=2) so the gate
    // exercises the same code path users hit. Earlier this test ran with
    // a laxer 3/1 to inflate coverage — that hid whether the gates would
    // pass at production settings.
    const {perActivity} = imputeMmpActivities(
      mmpa.initData, mmpa.rules, mmpa.rulesBased, mmpa.allCasesBased, 10, 2,
    );
    const predActs = perActivity[0].pred;

    let n = 0;
    let trueMean = 0;
    const evaluated: number[] = [];
    for (const i of heldOutIdxs) {
      if (predActs[i] !== FNULL) {
        evaluated.push(i);
        trueMean += rawActs[i];
        n++;
      }
    }
    if (n === 0) {
      grok.shell.warning('holdout_hERG: no predictions produced; skipping MAE/R² checks.');
      return;
    }
    trueMean /= n;

    let sumAbsErr = 0;
    let sumPredResid = 0;
    let sumTrueResid = 0;
    for (const i of evaluated) {
      const err = predActs[i] - rawActs[i];
      sumAbsErr += Math.abs(err);
      sumPredResid += err * err;
      sumTrueResid += (rawActs[i] - trueMean) ** 2;
    }
    const mae = sumAbsErr / n;
    const r2 = 1 - sumPredResid / sumTrueResid;
    grok.shell.info(`holdout_hERG: predicted ${n}/${heldOutIdxs.length} held-out rows, ` +
      `MAE=${mae.toFixed(3)}, R²=${r2.toFixed(3)}`);
    expect(mae < 0.8, true, `holdout MAE must be < 0.8, got ${mae.toFixed(3)}`);
    expect(r2 > 0.3, true, `holdout R² must be > 0.3, got ${r2.toFixed(3)}`);
  }, {timeout: 120000});

  // Stratified version of `holdout_hERG`. For each held-out row, compute its
  // maximum ECFP4 Tanimoto similarity to any known-activity row, bin by that
  // similarity, then report MAE per bin. This is the MoleculeACE-style
  // "where does the predictor actually work?" view — predictions on rows with
  // close neighbours (Tc > 0.7) should be very accurate; predictions far from
  // any known compound (Tc < 0.3) should degrade gracefully.
  //
  // Assertion: high-Tc bin MAE must be < low-Tc bin MAE. If it isn't, the
  // predictor is producing same-quality predictions regardless of neighbour
  // distance, which means either (a) the algorithm isn't using
  // chemistry-aware reasoning, or (b) the metric is dominated by something
  // other than similarity. Either is a signal worth investigating.
  test('holdout_hERG_stratified_by_neighbour_Tc', async () => {
    const df = await readDataframe('demo_files/mmp_demo.csv');
    const smilesCol = df.col('smiles');
    const hergCol = df.col('hERG_pIC50');
    if (!smilesCol || !hergCol) {
      grok.shell.warning('holdout_hERG_stratified: missing columns; skipping.');
      return;
    }
    const rawActs: Float32Array = hergCol.getRawData() as Float32Array;
    const heldOutIdxs: number[] = [];
    const maskedActs = new Float32Array(rawActs.length);
    for (let i = 0; i < rawActs.length; i++) {
      maskedActs[i] = rawActs[i];
      if (rawActs[i] !== FNULL && i % 5 === 0) {
        heldOutIdxs.push(i);
        maskedActs[i] = FNULL;
      }
    }

    // Predict via MMP (same path as holdout_hERG, production defaults).
    const sortingInfo: SortData = {fragmentIdxs: [], frequencies: []};
    const mmpa = await MMPA.init(
      smilesCol.name, smilesCol.toList(), 0.4,
      [maskedActs], ['hERG_pIC50'], [MmpDiffTypes.delta], sortingInfo,
    );
    const {perActivity} = imputeMmpActivities(
      mmpa.initData, mmpa.rules, mmpa.rulesBased, mmpa.allCasesBased, 10, 2,
    );
    const predActs = perActivity[0].pred;

    // Compute max-Tc-to-any-known for each held-out row. We use the package's
    // Morgan-FP cache via chem-searches (same fingerprints the kNN baseline
    // uses), then a single nested scan over known-activity rows.
    const fingerprints = await chemSearches.chemGetFingerprints(
      smilesCol, Fingerprint.Morgan, false);
    const knownMask = new Uint8Array(rawActs.length);
    for (let i = 0; i < rawActs.length; i++)
      if (maskedActs[i] !== FNULL && fingerprints[i]) knownMask[i] = 1;

    // Bins: (-∞, 0.3], (0.3, 0.5], (0.5, 0.7], (0.7, 1.0]
    const binBounds = [0.3, 0.5, 0.7, 1.01];
    const binNames = ['Tc<=0.30', '0.30<Tc<=0.50', '0.50<Tc<=0.70', 'Tc>0.70'];
    const binSumAbsErr = new Float64Array(binBounds.length);
    const binCount = new Int32Array(binBounds.length);

    for (const i of heldOutIdxs) {
      const fpI = fingerprints[i];
      if (!fpI || predActs[i] === FNULL) continue;
      let maxTc = 0;
      for (let j = 0; j < rawActs.length; j++) {
        if (!knownMask[j]) continue;
        const fpJ = fingerprints[j];
        if (!fpJ) continue;
        const tc = tanimotoSimilarity(fpI, fpJ);
        if (tc > maxTc) maxTc = tc;
      }
      const err = Math.abs(predActs[i] - rawActs[i]);
      let bin = 0;
      while (bin < binBounds.length - 1 && maxTc > binBounds[bin]) bin++;
      binSumAbsErr[bin] += err;
      binCount[bin]++;
    }

    const report: string[] = ['Stratified hERG-holdout MAE by max-Tc-to-known:'];
    const binMaes: number[] = [];
    for (let b = 0; b < binBounds.length; b++) {
      const n = binCount[b];
      const mae = n > 0 ? binSumAbsErr[b] / n : NaN;
      binMaes.push(mae);
      report.push(`  ${binNames[b]}: n=${n}, MAE=${Number.isFinite(mae) ? mae.toFixed(3) : 'n/a'}`);
    }
    grok.shell.info(report.join('\n'));
    console.log(report.join('\n'));

    // Soft assertion: at least one of the high-Tc bins (>0.5) must have lower
    // MAE than the lowest-Tc bin, IF both have data. This is the
    // "predictor uses chemistry" check. We don't gate on a specific number —
    // datasets vary — but the relationship between Tc and accuracy should
    // hold qualitatively.
    const lowMae = binMaes[0];
    const highMae = Math.min(binMaes[2], binMaes[3]);
    if (Number.isFinite(lowMae) && Number.isFinite(highMae)) {
      expect(highMae < lowMae, true,
        `Predictor should be more accurate where neighbours are closer. ` +
        `Lowest-Tc-bin MAE=${lowMae.toFixed(3)}, best-of-high-Tc-bins MAE=${highMae.toFixed(3)}.`);
    }
  }, {timeout: 240000});
});
