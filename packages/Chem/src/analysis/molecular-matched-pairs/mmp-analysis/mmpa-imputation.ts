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
 * Note: this is NOT a leave-one-out estimate; the rule's meanDiff still
 * includes the target molecule's own pairs, so the leak is mild but not zero.
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
          const pred = actFs + meanDiff;
          sumPred[k][ss] += pred;
          sumSq[k][ss] += pred * pred;
          nAnchors[k][ss]++;
        }

        if (!ssNull && (fsNull || predictKnown)) {
          const pred = actSs - meanDiff;
          sumPred[k][fs] += pred;
          sumSq[k][fs] += pred * pred;
          nAnchors[k][fs]++;
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

    for (let i = 0; i < numMolecules; i++) {
      const n = nAnchors[k][i];
      if (n < minAnchors) continue;
      if (!predictKnown && actCol[i] !== FNULL) continue;

      const mean = sumPred[k][i] / n;
      pred[i] = mean;
      support[i] = n;

      const variance = sumSq[k][i] / n - mean * mean;
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
