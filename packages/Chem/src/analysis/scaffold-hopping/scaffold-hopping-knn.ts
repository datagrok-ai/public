/** kNN-on-Morgan-FP activity-prediction baseline. For each row in
 *  `predictTargets`, take the k=10 nearest known-activity neighbours by
 *  ECFP4 Tanimoto similarity (computed live, no caching), and predict the
 *  row's activity as the Tc-weighted mean of those neighbours' activities.
 *  The Tc-weighted std-dev of the same 10 activities is emitted as a
 *  per-row confidence proxy (see comment inside).
 *
 *  Reference baseline: Sheridan 2004 (DOI 10.1021/ci049782w) established
 *  k=10 with Tc weighting as the standard non-parametric chem-activity
 *  predictor against which other models are compared.
 *
 *  Extracted from `scaffold-hopping.ts` as part of the multi-module split;
 *  no behaviour change. */

import * as DG from 'datagrok-api/dg';

import * as chemSearches from '../../chem-searches';
import {Fingerprint} from '../../utils/chem-common';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

/** Per-row kNN prediction output. Three Float32Arrays of length N (the input
 *  table's row count). Rows that didn't get a prediction (no fingerprint, no
 *  neighbours, or excluded from `predictTargets`) carry `DG.FLOAT_NULL` /
 *  zero — caller checks against FNULL to know whether a value is real. */
export interface KnnPredictionResult {
  pred: Float32Array;
  k: Int32Array;
  stdev: Float32Array;
}

/** Runs the kNN-on-Morgan-FP prediction.
 *
 *  @param molecules        molecule column (SMILES).
 *  @param acts             length-N Float32Array of measured activities;
 *                          slots with `DG.FLOAT_NULL` are treated as
 *                          unknown (the prediction targets / the rows
 *                          where measured value is missing) and excluded
 *                          from the neighbour pool.
 *  @param N                table row count.
 *  @param predictTargets   set of row indices to PREDICT for. Restricting
 *                          to survivors avoids O(N * |known|) work where
 *                          most output is discarded (the proximity loop
 *                          only reads predictions for survivors).
 *  @param progress         progress indicator for cancel checks.
 *
 *  Throws if cancelled mid-loop. Catch upstream and degrade gracefully. */
export async function computeKnnFpPrediction(
  molecules: DG.Column, acts: Float32Array, N: number,
  predictTargets: Set<number>,
  progress: DG.TaskBarProgressIndicator,
): Promise<KnnPredictionResult> {
  const FNULL = DG.FLOAT_NULL;
  const KNN_K = 10;
  const fingerprints = await chemSearches.chemGetFingerprints(
    molecules, Fingerprint.Morgan, false);
  const knownIdxs: number[] = [];
  for (let i = 0; i < N; i++) {
    if (acts[i] !== FNULL && fingerprints[i]) knownIdxs.push(i);
  }
  const pred = new Float32Array(N).fill(FNULL);
  const k = new Int32Array(N);
  const stdev = new Float32Array(N).fill(FNULL);

  // Restrict prediction TARGETS to the caller-provided set — non-survivors
  // are filtered out of the displayed result anyway, and the proximity loop
  // only reads predictions for survivors. Predicting for all N rows was
  // O(N * |known|) work where most of it was discarded; restricting to
  // survivors brings worst-case from ~25 min (50k×5k) down to a few
  // seconds. The neighbour POOL stays full (all known-activity rows in
  // the table) so each prediction is unbiased — only the set we compute
  // for shrinks.
  if (knownIdxs.length > 0) {
    const topTc = new Float32Array(KNN_K);
    const topAct = new Float32Array(KNN_K);
    let kIter = 0;
    for (const i of predictTargets) {
      // Cooperative cancel + yield to event loop every 128 rows
      // so the browser stays responsive on big survivor pools.
      if ((kIter++ & 127) === 0) {
        if (progress.canceled) throw new Error('Scaffold hopping cancelled by user.');
        await new Promise((r) => setTimeout(r, 0));
      }
      const fpI = fingerprints[i];
      if (!fpI) continue;
      let topLen = 0;
      let minTc = Infinity;
      let minPos = 0;
      for (const j of knownIdxs) {
        if (i === j) continue;
        const fpJ = fingerprints[j];
        if (!fpJ) continue;
        const tc = tanimotoSimilarity(fpI, fpJ);
        if (topLen < KNN_K) {
          topTc[topLen] = tc;
          topAct[topLen] = acts[j];
          if (tc < minTc) {minTc = tc; minPos = topLen;}
          topLen++;
        } else if (tc > minTc) {
          topTc[minPos] = tc;
          topAct[minPos] = acts[j];
          minTc = Infinity;
          for (let m = 0; m < KNN_K; m++) {
            if (topTc[m] < minTc) {minTc = topTc[m]; minPos = m;}
          }
        }
      }
      if (topLen > 0) {
        let sumW = 0;
        let sumWV = 0;
        for (let m = 0; m < topLen; m++) {
          sumW += topTc[m];
          sumWV += topTc[m] * topAct[m];
        }
        if (sumW > 0) {
          const mean = sumWV / sumW;
          pred[i] = mean;
          k[i] = topLen;
          // Tc-weighted std-dev of the k neighbour activities —
          // "how much do my nearest neighbours disagree?". A
          // cheap, instance-local confidence proxy in the same
          // units as the prediction, directly comparable to the
          // MMP-anchor stdev so the combined Stdev column stays
          // semantically homogeneous. Not a calibrated posterior
          // (a Tanimoto-kernel GP would be) — just an honest
          // signal of neighbour agreement.
          let sumWDD = 0;
          for (let m = 0; m < topLen; m++) {
            const d = topAct[m] - mean;
            sumWDD += topTc[m] * d * d;
          }
          const variance = sumWDD / sumW;
          stdev[i] = variance > 0 ? Math.sqrt(variance) : 0;
        }
      }
    }
  }
  return {pred, k, stdev};
}
