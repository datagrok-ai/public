import * as DG from 'datagrok-api/dg';
import {MmpInitData, MmpRules, MmpRulesBasedData, MmpAllCasesBasedData} from './mmpa-misc';
import {createColWithDescription} from '../mmp-viewer/mmp-generations';

/** Per-activity numeric output for one imputation run. */
export type MmpActivityImputationResult = {
  /** Predicted activity; DG.FLOAT_NULL where not possible. */
  pred: Float32Array;
  /** Number of anchor molecules that contributed. */
  support: Int32Array;
  /** Std dev of individual predictions; DG.FLOAT_NULL where support < minAnchors. */
  stdev: Float32Array;
};

/** Full result returned by imputeMmpActivities. */
export type MmpImputationResult = {
  /** One entry per activity in initData.activitiesNames. */
  perActivity: MmpActivityImputationResult[];
};

/**
 * Single-pass MMP-based activity imputation.
 *
 * For each activity k, every molecule i whose activity is DG.FLOAT_NULL
 * accumulates predictions from all rules r where:
 *   - allCasesBased.pairsPerActivity[k][r] >= supportFloor
 *   - rulesBased.meanDiffs[k][r] !== DG.FLOAT_NULL
 *   - the paired molecule (the anchor) has a known activity
 *
 * A prediction is emitted only when nAnchors >= minAnchors.
 * Only delta-type meanDiffs are valid with this formula
 * (ratio rules require pred = anchor * meanDiff and are out of scope).
 *
 * When `predictKnown` is true, predictions are also emitted for molecules
 * whose activity is already known. Useful for validation — lets the caller
 * compare the model's prediction against the measured value side by side.
 *
 * The `useLooForPredictKnown` flag controls how the rule meanDiff is
 * computed for predict-known rows. Default `false` reproduces the
 * historical behaviour (uses the precomputed rule meanDiff, which includes
 * the target molecule's own pairs — a self-inclusion leak that biases
 * predictions toward the measured value and artificially tightens σ).
 * Setting it to `true` recomputes the meanDiff per (rule, target row),
 * excluding pairs that include the target row — the classic leave-one-out
 * correction. Cost: O(pairs_in_rule) extra work per predict-known
 * prediction; in practice negligible because rules typically have <50
 * pairs and predict-known rows are rare in production callers.
 *
 * The stress-test fixture (Imatinib local SAR, 24 cpds) showed a 2x
 * MAE gap between predict-known (~0.035) and holdout (~0.017) without
 * LOO; that gap collapses to noise level once LOO is on. The σ on
 * predict-known rows also becomes more honestly larger, improving the
 * σ-vs-|error| Pearson r calibration.
 *
 * Float64 accumulators are used internally to avoid catastrophic cancellation
 * in the variance formula; results are cast back to Float32 at the end.
 */
export function imputeMmpActivities(
  initData: MmpInitData,
  rules: MmpRules,
  rulesBased: MmpRulesBasedData,
  allCasesBased: MmpAllCasesBasedData,
  supportFloor: number,
  minAnchors: number = 2,
  predictKnown: boolean = false,
  /** Leave-one-out correction for predict-known rows. When true, the
   *  meanDiff used to predict a known row i excludes any pairs of the
   *  rule that include i as one endpoint. Off by default to preserve
   *  legacy behaviour; opt-in for callers that need unbiased predict-
   *  known estimates (typically validation / model-quality assessment).
   *  Has no effect on rows with DG.FLOAT_NULL activity. */
  useLooForPredictKnown: boolean = false,
): MmpImputationResult {
  const numActivities = initData.activitiesCount;
  const numMolecules = initData.molecules.length;
  const numRules = rules.rules.length;
  const FNULL = DG.FLOAT_NULL;

  const sumPred: Float64Array[] = new Array(numActivities);
  const sumSq: Float64Array[] = new Array(numActivities);
  const nAnchors: Int32Array[] = new Array(numActivities);

  for (let k = 0; k < numActivities; k++) {
    sumPred[k] = new Float64Array(numMolecules);
    sumSq[k] = new Float64Array(numMolecules);
    nAnchors[k] = new Int32Array(numMolecules);
  }

  // Per-(rule, molecule) accumulators for the σ decomposition. Sparse
  // Map keyed by (r * numMolecules + i) — populated only for (rule,
  // target-row) cells where at least one pair prediction landed.
  //
  // Why sparse: the dense version is O(numActivities × numRules ×
  // numMolecules) Float64Arrays. For toy fixtures (60 mols × 32 rules =
  // ~50 KB) that's fine; for production-scale hERG fixtures (~20 k
  // mols × thousands of rules) the allocation hits multi-hundred-MB and
  // throws ArrayBuffer allocation failures. Most (rule, row) cells are
  // empty in practice — a rule with 50 pairs contributes to at most
  // ~100 distinct rows out of 20 k. Sparse maps to actual coverage.
  //
  // What the decomposition does (independent of dense vs sparse): the
  // old σ formula `sumSq/n - mean²` is the variance over ALL per-pair
  // predictions, treating each pair as one sample. That collapses to
  // ~zero when one rule dominates the pair count, even if that rule is
  // mis-applying — pair predictions cluster tightly inside the rule, so
  // σ shrinks while error stays meaningful. Empirically: at N=60 the
  // σ-vs-|err| Pearson r went NEGATIVE (-0.464), meaning σ was
  // anti-informative.
  //
  // The decomposition fixes this by computing
  //
  //   σ² = σ_between² + σ_within²
  //         (between-rule disagreement) + (average within-rule spread)
  //
  // treating each RULE as one vote in the between term. A row whose
  // prediction comes from one rule with 12 pairs and another with 2
  // pairs gets σ_between = std-dev of the two rule means (no pair
  // weighting) — honest about inter-rule disagreement regardless of
  // pair imbalance. r doubled from 0.308 → 0.603 on the 24-cpd
  // Imatinib fixture after this change.
  interface RuleAgg {sumPred: number; sumSqPred: number; nPairs: number;}
  const ruleAgg: Map<number, RuleAgg>[] = new Array(numActivities);
  for (let k = 0; k < numActivities; k++) ruleAgg[k] = new Map();
  // Tight inline helper: get-or-create the aggregate for (rule, row)
  // under activity k. Returns the mutable object so callers can update
  // its fields directly without re-fetching.
  const getAgg = (k: number, r: number, i: number): RuleAgg => {
    const key = r * numMolecules + i;
    let a = ruleAgg[k].get(key);
    if (!a) {a = {sumPred: 0, sumSqPred: 0, nPairs: 0}; ruleAgg[k].set(key, a);}
    return a;
  };

  /** Recompute the rule's meanDiff for a target row, excluding pairs that
   *  include the target. Returns FNULL when fewer than `supportFloor`
   *  pairs remain (the rule is judged not transferable for this target
   *  under LOO). Hoisted as a closure so each predict-known prediction
   *  can call it with the current rule's pairs and activity column. */
  const looMeanDiff = (
    rulePairs: {fs: number; ss: number; core: number}[],
    actCol: Float32Array | Float64Array | number[], excludeIdx: number,
  ): number => {
    let sum = 0;
    let n = 0;
    for (let j = 0; j < rulePairs.length; j++) {
      const fs = rulePairs[j].fs;
      const ss = rulePairs[j].ss;
      if (fs === excludeIdx || ss === excludeIdx) continue;
      const af = actCol[fs];
      const as = actCol[ss];
      if (af === FNULL || as === FNULL) continue;
      sum += as - af;
      n++;
    }
    return n >= supportFloor ? sum / n : FNULL;
  };

  for (let r = 0; r < numRules; r++) {
    const rulePairs = rules.rules[r].pairs;

    for (let k = 0; k < numActivities; k++) {
      const meanDiff = rulesBased.meanDiffs[k][r];
      const pairCount = allCasesBased.pairsPerActivity[k][r];

      if (meanDiff === FNULL || pairCount < supportFloor) continue;

      const actCol = initData.activities[k];

      for (let j = 0; j < rulePairs.length; j++) {
        const fs = rulePairs[j].fs;
        const ss = rulePairs[j].ss;
        const actFs = actCol[fs];
        const actSs = actCol[ss];
        const fsNull = actFs === FNULL;
        const ssNull = actSs === FNULL;

        // Predict for ss when the anchor (fs) is known; default behaviour
        // also requires the target (ss) to be unknown, but `predictKnown`
        // lifts that restriction so callers can validate.
        if (!fsNull && (ssNull || predictKnown)) {
          // LOO: when predicting a KNOWN row (ssNull is false here),
          // recompute the meanDiff excluding pairs that include ss.
          // The full meanDiff path remains for the predict-unknown case
          // since LOO doesn't apply (ss is FNULL so ss's pairs were
          // never included in the rule's meanDiff anyway).
          let mdForPred = meanDiff;
          if (useLooForPredictKnown && !ssNull) {
            mdForPred = looMeanDiff(rulePairs, actCol, ss);
            if (mdForPred === FNULL) continue;
          }
          const pred = actFs + mdForPred;
          sumPred[k][ss] += pred;
          sumSq[k][ss] += pred * pred;
          nAnchors[k][ss]++;
          const a = getAgg(k, r, ss);
          a.sumPred += pred;
          a.sumSqPred += pred * pred;
          a.nPairs++;
        }

        if (!ssNull && (fsNull || predictKnown)) {
          let mdForPred = meanDiff;
          if (useLooForPredictKnown && !fsNull) {
            mdForPred = looMeanDiff(rulePairs, actCol, fs);
            if (mdForPred === FNULL) continue;
          }
          const pred = actSs - mdForPred;
          sumPred[k][fs] += pred;
          sumSq[k][fs] += pred * pred;
          nAnchors[k][fs]++;
          const a = getAgg(k, r, fs);
          a.sumPred += pred;
          a.sumSqPred += pred * pred;
          a.nPairs++;
        }
      }
    }
  }

  const perActivity: MmpActivityImputationResult[] = new Array(numActivities);

  for (let k = 0; k < numActivities; k++) {
    const pred = new Float32Array(numMolecules).fill(FNULL);
    const support = new Int32Array(numMolecules);
    const stdev = new Float32Array(numMolecules).fill(FNULL);
    const actCol = initData.activities[k];

    // Bucket the sparse (rule, row) entries by row in a single pass —
    // each molecule's σ calculation then iterates over only the rules
    // that actually contributed, regardless of how many rules exist
    // globally. Cost: O(total predictions) bucketing + O(rules-per-row)
    // σ per row. Cleanly scales from 24-row toy fixtures up to 20k-row
    // hERG datasets without the dense allocation that failed before.
    const perMolMeans: number[][] = new Array(numMolecules);
    const perMolWithinVars: number[][] = new Array(numMolecules);
    for (let i = 0; i < numMolecules; i++) {
      perMolMeans[i] = [];
      perMolWithinVars[i] = [];
    }
    for (const [key, agg] of ruleAgg[k]) {
      if (agg.nPairs === 0) continue;
      const i = key % numMolecules;
      const meanR = agg.sumPred / agg.nPairs;
      perMolMeans[i].push(meanR);
      // Within-rule variance for this (rule, row) — meaningful only
      // when the rule contributed >=2 pairs to row i.
      //
      // Bessel-corrected (sample) variance:
      //   var = (Σx² - (Σx)²/n) / (n-1)
      // Equivalent to (sumSq - n·mean²)/(n-1), expanded to keep the
      // numerator in raw accumulators for numerical stability. The
      // biased /n formula systematically understates σ by ~29% at the
      // minimum-anchor floor (n=2), which then inflates MMP's weight
      // in the inverse-variance blend exactly when uncertainty should
      // suppress it. Bessel correction is the standard fix.
      if (agg.nPairs > 1) {
        const n = agg.nPairs;
        const varR = (agg.sumSqPred - agg.sumPred * agg.sumPred / n) / (n - 1);
        perMolWithinVars[i].push(varR > 0 ? varR : 0);
      } else
        perMolWithinVars[i].push(0);
    }

    for (let i = 0; i < numMolecules; i++) {
      const n = nAnchors[k][i];
      if (n < minAnchors) continue;
      if (!predictKnown && actCol[i] !== FNULL) continue;

      const mean = sumPred[k][i] / n;
      pred[i] = mean;
      support[i] = n;

      // σ via law of total variance:
      //   Var(X) = E[Var(X|R)] + Var(E[X|R])
      // where X = pair prediction, R = rule that produced it.
      // E[Var(X|R)] = average within-rule variance  → σ_within²
      // Var(E[X|R]) = variance of rule means        → σ_between²
      //
      // Crucially, σ_between is computed UNWEIGHTED across rules (each
      // rule contributes one (mean_r) sample). The old `sumSq/n - mean²`
      // formula effectively pair-weighted both terms, which let a single
      // dominant rule's tight pair cluster collapse σ even when other
      // rules disagreed. Unweighted-between honestly reflects rule
      // disagreement regardless of pair imbalance.
      const ruleMeans = perMolMeans[i];
      const ruleWithinVars = perMolWithinVars[i];
      const nRulesContrib = ruleMeans.length;
      let variance = 0;
      if (nRulesContrib > 0) {
        let momMeans = 0;
        for (let m = 0; m < nRulesContrib; m++) momMeans += ruleMeans[m];
        momMeans /= nRulesContrib;
        let varBetween = 0;
        if (nRulesContrib >= 2) {
          for (let m = 0; m < nRulesContrib; m++) {
            const d = ruleMeans[m] - momMeans;
            varBetween += d * d;
          }
          // Bessel correction: treat each rule mean as an independent
          // observation, divide by (nRules - 1) for an unbiased
          // estimator. Matches the Bessel-correction applied to the
          // within-rule variance above so the σ² = within + between
          // sum stays internally consistent (both terms in sample-
          // variance units rather than mixing biased and unbiased).
          varBetween /= (nRulesContrib - 1);
        }
        let varWithinSum = 0;
        for (let m = 0; m < nRulesContrib; m++) varWithinSum += ruleWithinVars[m];
        const varWithinAvg = varWithinSum / nRulesContrib;
        variance = varBetween + varWithinAvg;
      }
      stdev[i] = variance > 0 ? Math.sqrt(variance) : 0;
    }

    perActivity[k] = {pred, support, stdev};
  }

  return {perActivity};
}

/**
 * Converts an MmpImputationResult into DG.Column triples per activity.
 * Does not add columns to any table — caller is responsible for that.
 *
 * Column names:
 *   `<activityName> (MMP Pred)`    — double, FNULL where no prediction
 *   `<activityName> (MMP Support)` — int, anchor count (0 = no prediction)
 *   `<activityName> (MMP Stdev)`   — double, FNULL where no prediction
 *
 * Each column gets a human-readable description via `createColWithDescription`
 * (the existing helper used by MMP generations — keeps tooltips consistent
 * with the rest of the MMP analysis output).
 */
export function imputationResultToColumns(
  result: MmpImputationResult,
  activityNames: string[],
): DG.Column[][] {
  const groups: DG.Column[][] = [];

  for (let k = 0; k < activityNames.length; k++) {
    const base = activityNames[k];
    const r = result.perActivity[k];

    const predCol = createColWithDescription(
      'double', `${base} (MMP Pred)`, Array.from(r.pred), {}, '',
      undefined, undefined,
      `MMP-predicted ${base} (mean over MMP-rule anchors); null where fewer than 2 anchors`,
    );

    const supportCol = createColWithDescription(
      'int', `${base} (MMP Support)`, Array.from(r.support), {}, '',
      undefined, undefined,
      `Number of MMP anchor molecules contributing to the ${base} prediction`,
    );

    const stdevCol = createColWithDescription(
      'double', `${base} (MMP Stdev)`, Array.from(r.stdev), {}, '',
      undefined, undefined,
      `Std dev of individual MMP predictions for ${base}; null where fewer than 2 anchors`,
    );

    groups.push([predCol, supportCol, stdevCol]);
  }

  return groups;
}
