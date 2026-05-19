import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {addSubstructProvider, ISubstruct, ISubstructProvider}
  from '@datagrok-libraries/chem-meta/src/types';

import * as chemCommonRdKit from '../../utils/chem-common-rdkit';
import {removeWaterAndSaltsSingle} from '../../utils/reactions/reactions';
import {MMPA} from '../molecular-matched-pairs/mmp-analysis/mmpa';
import {imputeMmpActivities} from '../molecular-matched-pairs/mmp-analysis/mmpa-imputation';
import {MmpDiffTypes} from '../molecular-matched-pairs/mmp-function-editor';
import {SortData} from '../molecular-matched-pairs/mmp-viewer/mmp-viewer';
import {
  extractFragmentMolblockWithRGroups, buildCleanFragmentMolblock,
  buildR1AnchorMolblock,
} from './scaffold-hopping-molblock';
import {computeCatsCosine} from './scaffold-hopping-cats';
import {
  computeTanimoto, computeTopMcsRatio, TOP_N_FOR_MCS,
} from './scaffold-hopping-mcs';
import {computeKnnFpPrediction} from './scaffold-hopping-knn';

// sh-color moved to ./scaffold-hopping-mcs.ts

/** Composite-score weights. Sum to 1. Two complementary low-cost descriptors —
 *  ECFP4 Tanimoto (substructure / scaffold "shape" signal) and CATS2D cosine
 *  (Schneider 1999 topological pharmacophore-pair float vector with feature-
 *  count normalisation, purpose-built for low-Tc scaffold-hop retrieval).
 *
 *  Why no Pharm2D in the blend: empirically Pharm2D Tanimoto correlates ~0.97
 *  with ECFP4 on drug-like sets (Pearson on the BCR-ABL test rows), so adding
 *  it as a third descriptor inflates a signal we already have rather than
 *  contributing orthogonal information. CATS, by contrast, captures the
 *  pharmacophore-distance distribution that survives scaffold swaps — exactly
 *  the signal Schneider 1999 designed it for. ErG (Stiefl 2006) is a better
 *  candidate than Pharm2D for a future bit-vector pharmacophore signal.
 *
 *  Empirical basis for the 0.4/0.6 split: tuned by inspection on the 5-row
 *  BCR-ABL test set (Imatinib reference vs. Nilotinib/Dasatinib/Bosutinib/
 *  Ponatinib), NOT on held-out validation data. CATS values were uniformly
 *  0.92-0.97 across all four candidates (correctly identifying them as
 *  pharmacophore-similar to imatinib), while Tc varied 0.20-0.52, so giving
 *  CATS the heavier weight produced more stable hop-discovery ranking than
 *  the inverse split. The 0.6 also preserves the total pharmacophore weight
 *  from the previous Pharm2D-included blend (Tc 0.4 + Pharm 0.3 + CATS 0.3),
 *  which had been empirically functional. Users with different priors are
 *  expected to edit these constants directly; we deliberately don't surface
 *  them as user inputs because composite-weight tuning is a research task,
 *  not a per-run knob. */
const W_TANIMOTO = 0.4;
const W_CATS = 0.6;

// top-n moved to ./scaffold-hopping-mcs.ts

// timeout moved to ./scaffold-hopping-mcs.ts

// exact-atoms moved to ./scaffold-hopping-mcs.ts

// Pharmacophore-features.csv cache moved to ./scaffold-hopping-cats.ts
// exact-bonds moved to ./scaffold-hopping-mcs.ts

const COL_TANIMOTO = 'Scaffold Hop Tanimoto';
const COL_CATS = 'Scaffold Hop CATS Sim';
const COL_MCS_RATIO = 'Scaffold Hop MCS Ratio';
const COL_SCORE = 'Scaffold Hop Score';
const COL_FLAG = 'Scaffold Hop';
/** Per-row "where in the candidate is the conserved region" molecule
 *  column. Renders the FULL candidate SMILES with the shared-region atoms
 *  (atoms IN the MCS, restricted to the marked-region image when atoms
 *  were marked) highlighted in green, via a per-row substruct provider
 *  attached to `col.temp` — same per-row highlight mechanism MMP uses for
 *  from/to pair columns. Sister column to `Scaffold Hop Shared Region`,
 *  which shows the same shared atoms but as a standalone fragment SMILES;
 *  together the two columns answer "where in this candidate?" and "what
 *  does it look like alone?" from the same underlying MCS data. */
const COL_REASON = 'Scaffold Hop Reason';
/** Per-row activity-proximity factor in [0, 1], populated only when the
 *  caller passed an `activityColumnName`. 1.0 = same activity as the
 *  reference; values decay smoothly toward 0 as the candidate's activity
 *  drifts further from the reference, with the decay scale set to the
 *  activity column's standard deviation. Surfaced as a column so the user
 *  can sort/filter by activity proximity independently of the structural
 *  Score, and so the score-multiplier behaviour is transparent. */
const COL_ACTIVITY_PROX = 'Scaffold Hop Activity Proximity';
/** Per-row predicted activity, emitted when `imputeActivity` is on. Waterfall
 *  source: MMP-rule-based prediction first (mechanistic, sparse); falls back
 *  to Morgan-FP kNN (k=10, Tc-weighted mean of nearest known-activity
 *  neighbours; always-filled baseline). `Scaffold Hop Predicted Activity
 *  Source` says which engine fired per row, plus a "Δ=X vs kNN" suffix
 *  on MMP rows where the two engines disagree above the threshold.
 *  `Scaffold Hop Predicted Activity Stdev` is in the activity column's
 *  units regardless of which engine fired: MMP rows carry the anchor
 *  spread; kNN rows carry the Tc-weighted std-dev of the k neighbour
 *  activities. The column is suppressed when every value would be
 *  null-or-zero (degenerate / binary activity) to dodge a Datagrok
 *  filter-panel histogram NaN.floor() crash. */
const COL_PRED_ACTIVITY = 'Scaffold Hop Predicted Activity';
const COL_PRED_SOURCE = 'Scaffold Hop Predicted Activity Source';
const COL_PRED_STDEV = 'Scaffold Hop Predicted Activity Stdev';
/** Visible integer rank column. Drives the "reference → hops by score desc
 *  → survivors → dropped" view ordering via the standard grid sort (`grid.
 *  sort([COL_RANK], [true])` = ascending), and doubles as a "this is hop
 *  #N" readout for the user. Reference = 0; hops = 1..nHops; survivors =
 *  nHops+1..; dropped at the end. We tried `grid.setRowOrder` first but it
 *  silently no-ops on some Datagrok builds — column-sort is reliable. */
const COL_RANK = 'Scaffold Hop Rank';

const COL_REPLACEMENT = 'Scaffold Hop Replacement';
/** Per-row "what corresponds to your marked region" fragment column.
 *  Extracted as a standalone fragment SMILES. Two modes:
 *  - **Local preset (R-group decomposition):** the unmarked region is
 *    treated as an exact core; for each candidate that contains the
 *    core, the R-group(s) at the marked-region attachment positions
 *    are surfaced here. Candidates that don't contain the exact core
 *    get an empty cell — those aren't local hops by this definition.
 *  - **Other presets (ErG pharmacophore matching):** the candidate
 *    atoms matched to the user's marked reference atoms via ErG
 *    (Stiefl 2006 reduced-graph) — pharmacophore-equivalent rings can
 *    match where exact-MCS refuses, and path-completion fills the
 *    linker atoms in between.
 *  Renamed from the earlier "Shared Region" label, which suggested
 *  literal atom-level sharing — "Replaced Region" more accurately
 *  conveys that the column shows what the candidate has *in place of*
 *  the user's marked region. */
const COL_REPLACED = 'Scaffold Hop Replaced Region';
/** Legacy column name; cleaned up on re-run so users with older result
 *  tables don't end up with both the old and new columns side-by-side. */
const LEGACY_SHARED_COLS = ['Scaffold Hop Shared Region'];

/** Legacy column names from earlier iterations of this feature, cleaned up
 *  on re-run so users with old result tables don't end up with multiple
 *  pharmacophore columns sitting next to each other:
 *  - `Scaffold Hop Pharm Sim`     — v2 Pharm2D Tanimoto (dropped in favour of CATS-only)
 *  - `Scaffold Hop Pharm Overlap` — v1 7-bit family Jaccard
 *  - `Scaffold Hop Pharm Cosine`  — v1.5 brief: ErG cosine, never shipped
 *  - `Scaffold Hop TcMCS`         — earlier mis-applied bond-Tanimoto formula */
const LEGACY_PHARM_COLS = [
  'Scaffold Hop Pharm Sim',
  'Scaffold Hop Pharm Overlap',
  'Scaffold Hop Pharm Cosine',
];
const LEGACY_MCS_COLS = ['Scaffold Hop TcMCS'];

/** Orchestrates the 2D scaffold-hopping pipeline.
 *
 *  1. Standardize SMILES via `removeWaterAndSaltsSingle`.
 *  2. ECFP4 Tanimoto pre-filter — retain rows with `Tc ∈ [tcMin, tcMax]`.
 *  3. CATS2D cosine for survivors — Schneider 1999 topological-pharmacophore-
 *     pair descriptor (7×7×10 = 490-dim float vector with `count(A)+count(B)`
 *     normalisation per pair-and-distance bin), the canonical low-Tc
 *     scaffold-hop descriptor.
 *  4. Composite score = `0.4 × Tc + 0.6 × catsCosine`.
 *  5. MCS atom-ratio top-200 — Maeda 2024 SH classifier:
 *     `ratio_atom = atoms(MCS) / atoms(query)`.
 *  6. Scaffold-hop flag: `ratio_atom ≤ mcsRatioMax (default 0.4 per Maeda)`,
 *     plus optional Tc-window and CATS-Sim conditions toggleable by the user.
 *     If user marked atoms, additionally require the candidate's MCS does NOT
 *     cover all marked atoms.
 *  7. Output 7 columns (Tanimoto, CATS, MCS ratio, Score, Hop flag,
 *     Replacement, Shared Region) + apply pass/fail text color coding +
 *     apply hops-first row order (reference → flagged hops by score desc →
 *     other survivors by score desc → dropped rows by Tc desc).
 *
 *     The Replacement column shows the FULL candidate molecule with the
 *     shared atoms highlighted in green via a per-row substruct provider
 *     ("where in this candidate is the conserved region"). The Shared
 *     Region column shows the same shared atoms extracted as a standalone
 *     fragment SMILES ("what does the conserved region look like alone").
 *     With marked atoms both columns restrict to the marked-region image;
 *     without marks both columns show the entire candidate-side MCS.
 *     Fragment extraction reuses each candidate's V2000 molblock —
 *     boundary atoms get implicit hydrogens filled in by RDKit at parse
 *     time. */
export async function runScaffoldHopping(
  table: DG.DataFrame,
  molecules: DG.Column,
  referenceRowIdx: number,
  tanimotoMin: number,
  tanimotoMax: number,
  mcsRatioMax: number,
  minCatsSim: number,
  replaceableAtomsJson: string,
  useTcInFlag: boolean = true,
  useCatsInFlag: boolean = true,
  /** Activity-aware re-rank — empty disables. When set, must name a
   *  numeric column in `table`; each candidate's composite score is
   *  multiplied by `exp(-|activity - refActivity| / scale)` where `scale`
   *  is the column's std-dev, so candidates near the reference's activity
   *  rank higher. The unblended proximity factor is also surfaced as a
   *  separate `Scaffold Hop Activity Proximity` column for transparency
   *  and independent sorting. */
  activityColumnName: string = '',
  /** Blend strength of the activity proximity in the score multiplier,
   *  in [0, 1]. 0 leaves the score unchanged (proximity column still
   *  populated for inspection); 1 lets activity proximity dominate. */
  activityWeight: number = 0.3,
  /** When true, write the `Scaffold Hop Reason` string column listing the
   *  per-row pass/fail breakdown of each flag condition. Defaults to
   *  false — the reason column adds visual noise to the result table for
   *  users who just want the boolean flag + scores, so it's opt-in. */
  showReason: boolean = false,
  /** When true, the per-row Replacement / Replaced Region columns are
   *  computed via R-group decomposition (exact-core matching, clean
   *  attachment-point semantics). When false, the columns use ErG
   *  pharmacophore matching (allows ring chemotype swaps, includes
   *  shortest-path linker atoms). Local preset turns this on; Easy /
   *  Middle / Hard presets leave it off. */
  useRGroupReplacement: boolean = false,
  /** When true, run MMP-based activity imputation on `activityColumnName`
   *  before the proximity step. Adds prediction columns and lets the
   *  proximity loop use the filled (measured ∪ predicted) activity values.
   *  Ignored when no activity column is selected. */
  imputeActivity: boolean = false,
  /** Min number of supporting pairs an MMP rule must have to contribute
   *  to a prediction. mmpdb's preferred threshold is 10 (Dalke 2018,
   *  DOI 10.1021/acs.jcim.8b00173); 5 is the documented fallback for
   *  sparse datasets. */
  mmpSupportFloor: number = 10,
  /** When true, MMP also emits predictions for rows with measured
   *  activity (validation mode). Off in production to avoid the mild
   *  self-inclusion bias from the rule meanDiff including the target's
   *  own pairs. */
  predictKnown: boolean = false,
): Promise<void> {
  if (molecules.semType !== DG.SEMTYPE.MOLECULE) {
    grok.shell.error(`Column ${molecules.name} is not of Molecule semantic type`);
    return;
  }
  if (referenceRowIdx < 0 || referenceRowIdx >= table.rowCount) {
    grok.shell.error(`Reference row index ${referenceRowIdx} is out of range`);
    return;
  }
  const refSmilesRaw: string = molecules.get(referenceRowIdx);
  if (!refSmilesRaw) {
    grok.shell.error(`Reference row ${referenceRowIdx} has empty SMILES`);
    return;
  }

  const replaceableAtoms = parseReplaceableAtoms(replaceableAtomsJson);
  const refSmiles = removeWaterAndSaltsSingle(refSmilesRaw);
  const N = table.rowCount;

  const progress = DG.TaskBarProgressIndicator.create('Scaffold Hopping...', {cancelable: true});
  try {
    progress.update(5, 'Computing ECFP4 Tanimoto...');
    const tanimoto = await computeTanimoto(molecules, refSmiles, N);

    const survivorIdxs: number[] = [];
    for (let i = 0; i < N; i++) {
      if (i === referenceRowIdx) continue;
      const tc = tanimoto[i];
      if (tc >= tanimotoMin && tc <= tanimotoMax) survivorIdxs.push(i);
    }

    // NaN defaults so users can distinguish "wasn't computed" from "= 0".
    const catsCosine = new Float32Array(N).fill(NaN);
    const mcsRatio = new Float32Array(N).fill(NaN);
    const score = new Float32Array(N).fill(NaN);
    const isHop = new Uint8Array(N);

    // Reference-row pinning. Hoisted BEFORE the no-survivors early-return so
    // the reference always shows 1.0 across every metric — without this the
    // narrow / inverted-Tc-range case left the reference row's CATS / Score /
    // MCS as NaN, which looked like a per-metric failure rather than "the
    // reference is by definition perfectly similar to itself".
    score[referenceRowIdx] = 1.0;
    tanimoto[referenceRowIdx] = 1.0;
    catsCosine[referenceRowIdx] = 1.0;
    mcsRatio[referenceRowIdx] = 1.0;

    if (survivorIdxs.length === 0) {
      grok.shell.warning(
        `No candidates passed the Tanimoto pre-filter [${tanimotoMin}, ${tanimotoMax}]. ` +
        `Try widening the bounds.`);
      // No MCS step ran, so no per-row shared substructs / fragments to surface.
      // Reason column still populated so the user can see WHY each row was
      // dropped (pre-filtered + Tc value).
      const reasonEmpty = new Array<string>(N).fill('');
      const fmt = (v: number) => Number.isFinite(v) ? v.toFixed(2) : '—';
      for (let i = 0; i < N; i++) {
        if (i === referenceRowIdx) {reasonEmpty[i] = '(reference)'; continue;}
        const tc = tanimoto[i];
        reasonEmpty[i] = Number.isNaN(tc) ?
          'Tc — (could not compute)' :
          `Tc ${fmt(tc)} ✗ (pre-filtered; outside [${fmt(tanimotoMin)}, ${fmt(tanimotoMax)}])`;
      }
      const activityProxEmpty = new Float32Array(N).fill(NaN);
      writeOutputColumns(table, molecules, tanimoto, catsCosine, mcsRatio, score, isHop,
        reasonEmpty, activityProxEmpty, '',
        new Map(), new Map(), replaceableAtoms.length > 0, showReason,
        refSmiles, new Set(replaceableAtoms), referenceRowIdx);
      applyColorCoding(table, tanimotoMin, tanimotoMax, mcsRatioMax, minCatsSim);
      applyHopsFirstOrder(table, referenceRowIdx, isHop, score, tanimoto);
      return;
    }

    progress.update(25, 'Computing CATS2D pharmacophore-pair descriptors...');
    try {
      await computeCatsCosine(molecules, refSmiles, survivorIdxs, catsCosine);
    } catch (e: any) {
      grok.shell.warning(
        `CATS2D pharmacophore computation failure (${e?.message ?? e}). ` +
        `CATS similarity will be NaN; score will fall back to ECFP4 Tc only and ` +
        `the CATS-in-flag condition (if enabled) will be skipped. ` +
        `Check that the Chem Python environment has rdkit and numpy available.`);
    }

    progress.update(50, 'Computing composite score...');
    // CATS-aware blend: when CATS is NaN (Python failure) fall back to ECFP4
    // Tc only, so the run still produces a usable ranking instead of NaN.
    for (const i of survivorIdxs) {
      const tc = tanimoto[i];
      const ca = catsCosine[i];
      score[i] = Number.isNaN(ca) ? tc : (W_TANIMOTO * tc + W_CATS * ca);
    }

    // Activity-aware re-rank. Empty `activityColumnName` disables. Per-row
    // proximity factor lives in `activityProx`; when blending is enabled
    // (`activityWeight > 0`), the composite score is multiplied by
    // `(1 - activityWeight) + activityWeight * proximity` so a fully-far
    // candidate keeps `(1 - activityWeight)` of its original score and a
    // same-activity candidate keeps it all. The proximity values are
    // surfaced as `Scaffold Hop Activity Proximity` regardless of weight
    // so the user can sort by it independently.
    //
    // When `imputeActivity` is on, we first run MMP-rule-based prediction
    // on the activity column (same machinery the MMP analysis viewer's
    // Generations tab uses: rule meanDiffs from `MMPA.init` aggregated over
    // anchor pairs). Three columns get added to the table —
    // `<actName> (MMP Pred / Support / Stdev)` — and the proximity loop
    // uses `measured ∪ predicted` so masked / unknown rows still get
    // ranked.
    const activityProx = new Float32Array(N).fill(NaN);
    let activityColUsed = '';
    if (activityColumnName) {
      const activityCol = table.col(activityColumnName);
      if (activityCol && activityCol.type !== DG.COLUMN_TYPE.STRING) {
        const FNULL = DG.FLOAT_NULL;
        const acts = new Float32Array(N);
        for (let i = 0; i < N; i++) {
          const v = activityCol.get(i);
          acts[i] = (typeof v === 'number' && Number.isFinite(v)) ? v : FNULL;
        }

        // Filled = measured where known, predicted where missing.
        // Without imputation, filled == acts.
        let filled = acts;
        if (imputeActivity) {
          // Compute MMP + kNN predictions, then collapse into ONE waterfall
          // column (MMP > kNN > FNULL) plus a Source string column and a
          // Stdev column. Avoids the five-column visual clutter of emitting
          // MMP and kNN separately; the Source string keeps per-row
          // attribution so you can still tell which engine fired.
          let mmpPred: Float32Array | null = null;
          let mmpSupport: Int32Array | null = null;
          let mmpStdev: Float32Array | null = null;
          let knnPred: Float32Array | null = null;
          let knnK: Int32Array | null = null;
          let knnStdev: Float32Array | null = null;

          progress.update(48, 'MMP activity prediction...');
          if (progress.canceled) throw new Error('Scaffold hopping cancelled by user.');
          console.log(`[scaffold-hopping] MMP imputation: activity=${activityColumnName}, N=${N}`);
          try {
            const sortingInfo: SortData = {fragmentIdxs: [], frequencies: []};
            const mmpa = await MMPA.init(
              molecules.name, molecules.toList(), 0.4,
              [acts], [activityColumnName], [MmpDiffTypes.delta], sortingInfo,
            );
            console.log(`[scaffold-hopping] MMPA built: ${mmpa.rules.rules.length} rules`);
            const result = imputeMmpActivities(
              mmpa.initData, mmpa.rules, mmpa.rulesBased, mmpa.allCasesBased,
              mmpSupportFloor, 2, predictKnown,
            );
            mmpPred = result.perActivity[0].pred;
            mmpSupport = result.perActivity[0].support;
            mmpStdev = result.perActivity[0].stdev;
            const nAllPred = Array.from(mmpPred).filter((p) => p !== FNULL).length;
            console.log(`[scaffold-hopping] MMP predictions: ${nAllPred} total`);
          } catch (e: any) {
            console.error('[scaffold-hopping] MMP imputation failed', e);
            grok.shell.warning(
              `MMP activity prediction skipped - falling back to kNN only. ${e?.message ?? e}`);
          }

          // kNN-on-Morgan-FP baseline. Always-filled fallback for rows where
          // MMP couldn't find enough anchors: take the k nearest
          // known-activity neighbours by ECFP4 Tanimoto and average their
          // activity weighted by similarity. Same-row excluded.
          // k = 10 with Tc weighting is a standard kNN-FP regression
          // baseline for chem activity (Sheridan 2004,
          // DOI 10.1021/ci049782w).
          progress.update(56, 'kNN-Morgan-FP prediction...');
          console.log('[scaffold-hopping] kNN-FP prediction starting...');
          try {
            // Restrict prediction TARGETS to the survivor pool (plus the
            // reference) — see scaffold-hopping-knn.ts header for rationale.
            const knnTargets = new Set<number>(survivorIdxs);
            knnTargets.add(referenceRowIdx);
            const r = await computeKnnFpPrediction(molecules, acts, N, knnTargets, progress);
            knnPred = r.pred;
            knnK = r.k;
            knnStdev = r.stdev;
            const nKnnTotal = Array.from(knnPred).filter((p) => p !== FNULL).length;
            console.log(`[scaffold-hopping] kNN-FP predictions: ${nKnnTotal} total`);
          } catch (e: any) {
            console.error('[scaffold-hopping] kNN-FP prediction failed', e);
            grok.shell.warning(`kNN-FP prediction skipped. ${e?.message ?? e}`);
          }

          // Combine into 3 columns: Predicted Activity, Source, Stdev.
          // Source string keeps the per-row attribution without the column
          // clutter of separate MMP-Pred and kNN-Pred columns.
          const predValue = new Float32Array(N).fill(FNULL);
          const predSource = new Array<string>(N).fill('');
          const predStdev = new Float32Array(N).fill(FNULL);
          // Disagreement threshold for surfacing the MMP/kNN delta in the
          // Source string. 1.0 activity units aligns with the activity-
          // cliff literature's definition of a "meaningful SAR shift"
          // (Auer 2014, PMC3448951) and stays above typical inter-assay
          // noise (~0.2-0.5 pIC50 units, Kramer 2014). A lower threshold
          // would fire on routine measurement noise and train users to
          // ignore the flag. Binary-label disagreements (|Δ| = 1) also
          // pass this check exactly at the boundary.
          const DISAGREEMENT_THRESHOLD = 1.0;
          // σ floor for the inverse-variance blend. Without it, a zero-σ
          // engine (e.g. MMP on discrete binary activity, where all anchors
          // agree exactly) gets infinite weight and collapses the blend
          // to itself. 0.05 pIC50 is the order of inter-assay precision on
          // a single replicate — below this the data isn't telling us
          // anything meaningful about uncertainty anyway.
          const SIGMA_FLOOR = 0.05;
          let nMmpFired = 0;
          let nKnnFired = 0;
          let nBlended = 0;
          let nDisagreement = 0;
          for (let i = 0; i < N; i++) {
            const mmpHas = mmpPred && mmpPred[i] !== FNULL;
            const knnHas = knnPred && knnPred[i] !== FNULL;
            if (mmpHas && knnHas) {
              // Both engines fired — combine by inverse-variance blending.
              // Each prediction's weight is 1/σ², so the more confident
              // engine dominates. On Local-mode rows where MMP-σ is small
              // (anchors agree tightly) the blend stays ~95 % MMP — no
              // behavioural change vs the old waterfall. On non-Local rows
              // where MMP fires with few anchors and large σ, kNN can
              // genuinely contribute and pull the prediction toward its
              // (typically tighter) chemistry-local estimate.
              const mmpSigma = Math.max(mmpStdev![i], SIGMA_FLOOR);
              const knnSigma = Math.max(knnStdev![i], SIGMA_FLOOR);
              const wMmp = 1 / (mmpSigma * mmpSigma);
              const wKnn = 1 / (knnSigma * knnSigma);
              const wSum = wMmp + wKnn;
              const blendedPred = (mmpPred![i] * wMmp + knnPred![i] * wKnn) / wSum;
              const blendedVar = 1 / wSum;
              const mmpWeightPct = Math.round((wMmp / wSum) * 100);
              predValue[i] = blendedPred;
              predStdev[i] = Math.sqrt(blendedVar);
              // Source string surfaces who carried the blend and (when the
              // two engines genuinely disagreed) the magnitude of the
              // resolved disagreement.
              let sourceStr = `MMP+kNN (n=${mmpSupport![i]}, k=${knnK![i]}, ${mmpWeightPct}% MMP)`;
              const delta = Math.abs(mmpPred![i] - knnPred![i]);
              if (delta >= DISAGREEMENT_THRESHOLD) {
                sourceStr += `, Δ=${delta.toFixed(2)}`;
                nDisagreement++;
              }
              predSource[i] = sourceStr;
              nBlended++;
            } else if (mmpHas) {
              predValue[i] = mmpPred![i];
              predSource[i] = `MMP (n=${mmpSupport![i]})`;
              predStdev[i] = mmpStdev![i];
              nMmpFired++;
            } else if (knnHas) {
              predValue[i] = knnPred![i];
              predSource[i] = `kNN (k=${knnK![i]})`;
              if (knnStdev) predStdev[i] = knnStdev[i];
              nKnnFired++;
            }
          }

          if (table.col(COL_PRED_ACTIVITY)) table.columns.remove(COL_PRED_ACTIVITY);
          if (table.col(COL_PRED_SOURCE)) table.columns.remove(COL_PRED_SOURCE);
          if (table.col(COL_PRED_STDEV)) table.columns.remove(COL_PRED_STDEV);
          const predCol = DG.Column.fromFloat32Array(COL_PRED_ACTIVITY, predValue);
          predCol.setTag('description',
            `Predicted ${activityColumnName} per row. When both engines fire, ` +
            `predictions are combined by inverse-variance weighting (each ` +
            `engine's contribution weighted by 1/σ², so the more confident ` +
            `one dominates). When only one engine fires, that engine's ` +
            `prediction is used directly. Engines: MMP-rule-based (mechanistic, ` +
            `needs >=2 anchors from rules with >=5 pairs) and Morgan-FP kNN ` +
            `(k=10 ECFP4-Tanimoto nearest known-activity neighbours, ` +
            `Tc-weighted mean). The Source column shows which engines fired ` +
            `and how much each contributed.`);
          const sourceCol = DG.Column.fromStrings(COL_PRED_SOURCE, predSource);
          sourceCol.setTag('description',
            `Which engine produced the prediction for this row: ` +
            `"MMP+kNN (n=N, k=K, X% MMP)" = inverse-variance blend with X% ` +
            `weight on MMP; "MMP (n=N)" = MMP only (kNN didn't fire); ` +
            `"kNN (k=K)" = kNN only (MMP didn't have enough anchors); ` +
            `"" = no prediction. Trailing "Δ=X.XX" on blended rows means ` +
            `the two engines disagreed by at least ${DISAGREEMENT_THRESHOLD} ` +
            `activity units before being combined — a reliability flag worth ` +
            `filtering on.`);
          table.columns.add(predCol);
          table.columns.add(sourceCol);

          // Only emit the Stdev column when it carries information — i.e.
          // at least one MMP prediction actually had spread across its
          // anchors. On binary / discrete activity columns every MMP row
          // ends up with stdev == 0 (anchors agree exactly), and a column
          // that is all-FNULL-or-zero crashes Datagrok's filter-panel
          // histogram with `Unsupported operation: NaN.floor()` because
          // the histogram bin width collapses to 0. Skipping the column
          // in the degenerate case avoids the renderer crash without
          // losing real information (degenerate stdev = always 0 anyway).
          let hasNonZeroStdev = false;
          for (let i = 0; i < N; i++) {
            if (predStdev[i] !== FNULL && predStdev[i] > 0) {
              hasNonZeroStdev = true;
              break;
            }
          }
          if (hasNonZeroStdev) {
            const stdevCol = DG.Column.fromFloat32Array(COL_PRED_STDEV, predStdev);
            stdevCol.setTag('description',
              `Std-dev of the prediction. For MMP rows: spread across the ` +
              `rule anchors. For kNN rows: Tc-weighted std-dev of the k ` +
              `neighbour activities. Both are in the activity column's ` +
              `units, so they are directly comparable as a per-row ` +
              `confidence proxy (smaller = more agreement).`);
            table.columns.add(stdevCol);
          }

          // Fill measured-gaps in `filled[]` from predValue (MMP > kNN).
          let nFilled = 0;
          filled = new Float32Array(N);
          for (let i = 0; i < N; i++) {
            if (acts[i] !== FNULL) filled[i] = acts[i];
            else if (predValue[i] !== FNULL) {filled[i] = predValue[i]; nFilled++;}
            else filled[i] = FNULL;
          }
          console.log(`[scaffold-hopping] Pred fired: ${nBlended} blended (MMP+kNN), ` +
            `${nMmpFired} MMP-only, ${nKnnFired} kNN-only, ` +
            `filled ${nFilled} missing rows, ${nDisagreement} blended disagreements |Δ|>=${DISAGREEMENT_THRESHOLD}`);
          const disagreementSuffix = nDisagreement > 0 ?
            ` ${nDisagreement} blended row${nDisagreement === 1 ? '' : 's'} ` +
            `had |Δ|>=${DISAGREEMENT_THRESHOLD} between MMP and kNN (see Source column).` : '';
          grok.shell.info(
            `Predicted Activity: ${nBlended} rows blended (MMP+kNN), ` +
            `${nMmpFired} MMP-only, ${nKnnFired} kNN-only. ` +
            `${nFilled} previously-missing rows filled.${disagreementSuffix}`);
        }

        const refRawF = filled[referenceRowIdx];
        const refA = (refRawF !== FNULL && Number.isFinite(refRawF)) ? refRawF : NaN;
        if (!Number.isNaN(refA)) {
          // Use the column's std-dev as the decay scale; fall back to 1.0
          // when std-dev is 0 (single-value column) so the factor stays
          // numerically defined.
          let mean = 0;
          let n = 0;
          for (let i = 0; i < N; i++) {
            const v = filled[i];
            if (v !== FNULL && Number.isFinite(v)) {mean += v; n++;}
          }
          if (n > 0) mean /= n;
          let variance = 0;
          for (let i = 0; i < N; i++) {
            const v = filled[i];
            if (v !== FNULL && Number.isFinite(v)) {
              const d = v - mean;
              variance += d * d;
            }
          }
          const scale = n > 1 ? Math.sqrt(variance / (n - 1)) : 0;
          if (scale === 0) {
            grok.shell.warning(
              `Activity column "${activityColumnName}" has zero std-dev — ` +
              `every proximity factor will collapse to 1.0 and the score ` +
              `multiplier will be inert. Pick a column with real activity ` +
              `variance to get a meaningful re-rank.`);
          }
          const safeScale = scale > 0 ? scale : 1;
          activityColUsed = activityColumnName;
          activityProx[referenceRowIdx] = 1.0;
          for (const i of survivorIdxs) {
            const v = filled[i];
            if (v !== FNULL && Number.isFinite(v)) {
              const p = Math.exp(-Math.abs(v - refA) / safeScale);
              activityProx[i] = p;
              if (activityWeight > 0 && !Number.isNaN(score[i]))
                score[i] = score[i] * ((1 - activityWeight) + activityWeight * p);
            }
          }
        }
      }
    }

    progress.update(60, `MCS atom-ratio for top ${Math.min(survivorIdxs.length, TOP_N_FOR_MCS)} candidates...`);
    if (progress.canceled) throw new Error('Scaffold hopping cancelled by user.');
    const topByScore = survivorIdxs.slice().sort((a, b) => score[b] - score[a]);
    const topMcsRows = topByScore.slice(0, TOP_N_FOR_MCS);
    const markedRefSet = new Set(replaceableAtoms);
    const {refMcsAtoms: mcsResults, sharedSubstructs, sharedFragments} =
      await computeTopMcsRatio(
        molecules, refSmiles, topMcsRows, mcsRatio, markedRefSet,
        useRGroupReplacement, progress);

    // Reference-row populating. With marks: the reference row's Shared
    // Region cell shows the marked-region fragment (so the user can read
    // the column straight down: "what I marked → each candidate's
    // version of it"), and the Replacement cell shows the full reference
    // with the marked atoms highlighted in green for visual continuity
    // with the candidate-row highlights. Without marks: both cells stay
    // empty (would just duplicate the smiles column).
    if (replaceableAtoms.length > 0) {
      try {
        const rdKitModule = chemCommonRdKit.getRdKitModule();
        const refMol = rdKitModule.get_mol(refSmiles);
        if (refMol) {
          try {
            const refMolblock = refMol.get_molblock?.() ?? '';
            if (refMolblock) {
              // Use the with-RGroups variant so the reference's marked
              // region carries the same `R1` / `R2` labels every
              // candidate row gets. Then the column reads top-to-bottom
              // as "this is what I marked (with R1 here) → here's each
              // candidate's version of it (with R1 at the equivalent
              // position)."
              const refFrag = extractFragmentMolblockWithRGroups(
                refMolblock, markedRefSet, rdKitModule);
              if (refFrag) sharedFragments.set(referenceRowIdx, refFrag);
            }
            // Replacement-column highlight payload for the reference row.
            const refHighlightColor = DG.Color.hexToPercentRgb('#5BC85B');
            if (refHighlightColor && markedRefSet.size > 0) {
              const highlightAtomColors: {[k: number]: number[]} = {};
              for (const a of markedRefSet) highlightAtomColors[a] = [...refHighlightColor];
              sharedSubstructs.set(referenceRowIdx, {
                atoms: [...markedRefSet],
                bonds: [],
                highlightAtomColors,
                highlightBondColors: {},
              });
            }
          } finally {
            refMol.delete();
          }
        }
      } catch {/* leave reference cell empty on failure */}
    }

    // Flag composition. Maeda 2024 atom-ratio is ALWAYS applied (the paper's
    // canonical SH classifier). The two optional conditions — Tc window and
    // CATS Sim — are user-toggleable via `useTcInFlag` / `useCatsInFlag`,
    // letting users drop into "strict-Maeda" mode by turning both off. The
    // marked-atoms refinement is implicit (only active when atoms were
    // marked) and is independent of the toggles.
    //
    // Per-row Reason string is produced in parallel — same conditions,
    // formatted for human-readable "why was this flagged / not flagged"
    // diagnosis. Format: dot-separated condition tokens, each tagged with
    // its numeric value and a ✓/✗ glyph. Survivor rows that failed Maeda
    // get a fully-populated reason; pre-filtered rows get a short
    // "(pre-filtered)" string with the Tc value; the reference row gets
    // "(reference)"; rows whose MCS timed out get "MCS — (timed out)".
    const reason = new Array<string>(N).fill('');
    const survivorSet = new Set(survivorIdxs);
    const fmtNum = (v: number) => Number.isFinite(v) ? v.toFixed(2) : '—';

    for (let i = 0; i < N; i++) {
      if (i === referenceRowIdx) {reason[i] = '(reference)'; continue;}
      if (!survivorSet.has(i)) {
        const tc = tanimoto[i];
        if (Number.isNaN(tc))
          reason[i] = 'Tc — (could not compute)';
        else
          reason[i] =
            `Tc ${fmtNum(tc)} ✗ (pre-filtered; outside [${fmtNum(tanimotoMin)}, ${fmtNum(tanimotoMax)}])`;
        continue;
      }
      const tc = tanimoto[i];
      const ca = catsCosine[i];
      const m = mcsRatio[i];
      const parts: string[] = [];
      let pass = true;

      const tcInWindow = tc >= tanimotoMin && tc <= tanimotoMax;
      parts.push(`Tc ${fmtNum(tc)} ${tcInWindow ? '✓' : '✗'}`);
      if (useTcInFlag && !tcInWindow) pass = false;

      if (Number.isNaN(ca)) {
        parts.push('CATS — (Python failure)');
        // CATS-in-flag: NaN degrades to "Tc + Maeda" rather than blocking.
      } else {
        const catsOk = ca >= minCatsSim;
        parts.push(`CATS ${fmtNum(ca)} ${catsOk ? '✓' : '✗'}`);
        if (useCatsInFlag && !catsOk) pass = false;
      }

      if (Number.isNaN(m)) {
        parts.push('MCS — (timed out)');
        pass = false;
      } else {
        const mcsOk = m <= mcsRatioMax;
        parts.push(`MCS ${fmtNum(m)} ${mcsOk ? '✓' : '✗'}`);
        if (!mcsOk) pass = false;
      }

      if (replaceableAtoms.length > 0) {
        const refMcsAtoms = mcsResults.get(i);
        const markedPreserved = refMcsAtoms &&
          replaceableAtoms.every((a) => refMcsAtoms.has(a));
        if (markedPreserved) {
          parts.push('marked region preserved ✗');
          pass = false;
        } else if (refMcsAtoms) {
          parts.push('marked region changed ✓');
        } else {
          // MCS skipped — already covered by the "MCS —" token above.
        }
      }

      if (pass) isHop[i] = 1;
      reason[i] = parts.join(' · ');
    }

    progress.update(95, 'Writing result columns...');
    writeOutputColumns(table, molecules, tanimoto, catsCosine, mcsRatio, score, isHop,
      reason, activityProx, activityColUsed,
      sharedSubstructs, sharedFragments, replaceableAtoms.length > 0, showReason,
      refSmiles, markedRefSet, referenceRowIdx);
    applyColorCoding(table, tanimotoMin, tanimotoMax, mcsRatioMax, minCatsSim);
    applyHopsFirstOrder(table, referenceRowIdx, isHop, score, tanimoto);

    const flaggedCount = (isHop as any).reduce((s: number, x: number) => s + x, 0);
    grok.shell.info(
      `Scaffold Hopping done.\n` +
      `Reference: row ${referenceRowIdx} (${refSmiles})\n` +
      `Pre-filter survivors: ${survivorIdxs.length} / ${N - 1}\n` +
      `Top-200 MCS atom-ratio: ${mcsResults.size} pairs computed (${topMcsRows.length - mcsResults.size} timed out)\n` +
      `Flagged as scaffold hops: ${flaggedCount}\n` +
      `Replaceable region: ${replaceableAtoms.length === 0 ? 'none (any global hop)' :
        `${replaceableAtoms.length} marked atoms`}`);
  } finally {
    progress.close();
  }
}

// ---------------------------------------------------------------------------

function parseReplaceableAtoms(json: string): number[] {
  try {
    const parsed = JSON.parse(json || '[]');
    return Array.isArray(parsed) ? parsed.filter((n) => Number.isInteger(n)) : [];
  } catch {
    return [];
  }
}

// computeTanimoto + computeTopMcsRatio moved to ./scaffold-hopping-mcs.ts


function writeOutputColumns(
  table: DG.DataFrame, molecules: DG.Column,
  tanimoto: Float32Array, catsCosine: Float32Array,
  mcsRatio: Float32Array, score: Float32Array, isHop: Uint8Array,
  reason: string[],
  activityProx: Float32Array,
  activityColUsed: string,
  sharedSubstructs: Map<number, ISubstruct>,
  sharedFragments: Map<number, string>,
  hasMarkedRegion: boolean,
  showReason: boolean,
  refSmiles: string,
  markedRefAtoms: Set<number>,
  referenceRowIdx: number,
): void {
  // Drop legacy columns from earlier versions so the table doesn't accumulate
  // orphan columns when re-running across upgrades.
  for (const legacy of LEGACY_PHARM_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);
  for (const legacy of LEGACY_MCS_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);
  for (const legacy of LEGACY_SHARED_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);

  replaceFloatColumn(table, COL_TANIMOTO, tanimoto);
  replaceFloatColumn(table, COL_CATS, catsCosine);
  replaceFloatColumn(table, COL_MCS_RATIO, mcsRatio);
  replaceFloatColumn(table, COL_SCORE, score);
  replaceBoolColumn(table, COL_FLAG, isHop);
  // Opt-in Reason column. Defaults to off so the result table stays
  // compact for users who only want the boolean flag + numeric scores;
  // power users who need the per-row pass/fail breakdown can toggle it
  // on in the dialog. On re-run with the checkbox off, drop any column
  // left over from a prior run that had it on.
  if (showReason) {
    replaceStringColumn(table, COL_REASON, reason);
  } else {
    if (table.col(COL_REASON)) table.columns.remove(COL_REASON);
  }
  // Activity proximity column is only meaningful when an activity column
  // was selected. Without it the array is all-NaN; drop the column so the
  // table doesn't carry a row of empty cells.
  if (activityColUsed) {
    replaceFloatColumn(table, COL_ACTIVITY_PROX, activityProx);
  } else {
    if (table.col(COL_ACTIVITY_PROX)) table.columns.remove(COL_ACTIVITY_PROX);
  }
  // The Replacement and Shared Region columns answer "what corresponds to
  // the region I marked?" — without marks the answer is "the entire shared
  // scaffold" which mostly duplicates the input molecule + adds visual
  // noise (~80% of every candidate would be highlighted green). Show
  // the columns only when the user actually asked a region-specific
  // question. On re-run without marks, drop any columns left over from
  // a prior marked-region run so the table stays clean.
  if (hasMarkedRegion) {
    // Build the alignment scaffolds once for both columns:
    // - Replacement column: aligns full candidates to the reference's
    //   UNMARKED region (the conserved core). Every candidate shares
    //   this substructure (that's what makes them local hops), so the
    //   alignment succeeds across the column.
    // - Replaced Region column: aligns fragments to a MINIMAL anchor
    //   (R1 + its bonded atom). Tiny scaffold so substructure search
    //   succeeds on every variant ring (piperidine, morpholine, …)
    //   regardless of the atom-type mismatch with the reference's
    //   piperazine.
    const rdkit = chemCommonRdKit.getRdKitModule();
    let unmarkedRegionMolblock: string | null = null;
    let r1AnchorMolblock: string | null = null;
    let refMol: any = null;
    try {
      refMol = rdkit.get_mol(refSmiles);
      if (refMol) {
        const refMolblock = refMol.get_molblock?.() ?? '';
        const refNumAtoms = refMol.get_num_atoms?.() ?? 0;
        if (refMolblock && refNumAtoms > 0) {
          const unmarkedAtoms = new Set<number>();
          for (let i = 0; i < refNumAtoms; i++)
            if (!markedRefAtoms.has(i)) unmarkedAtoms.add(i);
          if (unmarkedAtoms.size > 0)
            unmarkedRegionMolblock = buildCleanFragmentMolblock(
              refMolblock, unmarkedAtoms, rdkit);
        }
      }
    } catch {/* alignment is best-effort */} finally {refMol?.delete();}

    const refMarkedFragmentMolblock = sharedFragments.get(referenceRowIdx);
    if (refMarkedFragmentMolblock)
      r1AnchorMolblock = buildR1AnchorMolblock(refMarkedFragmentMolblock);

    replaceMoleculeColumnWithProvider(table, molecules, COL_REPLACEMENT,
      sharedSubstructs, unmarkedRegionMolblock ?? undefined);
    replaceFragmentColumn(table, COL_REPLACED, sharedFragments,
      r1AnchorMolblock ?? undefined);
  } else {
    if (table.col(COL_REPLACEMENT)) table.columns.remove(COL_REPLACEMENT);
    if (table.col(COL_REPLACED)) table.columns.remove(COL_REPLACED);
  }
  applyColumnDescriptions(table, activityColUsed);
}

/** Attaches per-column descriptions so hovering the column header in the
 *  Datagrok grid surfaces a plain-English explanation of what the column
 *  means — the user doesn't have to remember that CATS is Schneider 1999
 *  or that the Maeda atom-ratio is one-sided over the reference. Each
 *  description is written to the `description` tag, the same surface
 *  Datagrok's column tooltip reads. Idempotent — re-running overwrites
 *  with the current text. */
function applyColumnDescriptions(table: DG.DataFrame, activityColUsed: string): void {
  const desc = (name: string, text: string) => {
    const col = table.col(name);
    if (col) col.setTag('description', text);
  };
  desc(COL_TANIMOTO,
    'ECFP4 Tanimoto similarity to the reference row (radius 2, 2048 bits). ' +
    'Used as the pre-filter shortlist — only rows inside the [Tc min, Tc max] ' +
    'window proceed to CATS / MCS scoring. Higher = more similar. ' +
    'Reference: Rogers & Hahn 2010 JCIM 50:742.');
  desc(COL_CATS,
    'CATS2D pharmacophore-pair cosine similarity to the reference (Schneider ' +
    '1999, Angew. Chem. 38:2894). CATS is a 490-dim float vector (7 pharmacophore ' +
    'families × 7 × 10 topological-distance bins) with count-based normalisation. ' +
    'Designed for low-Tc scaffold-hop retrieval — pharmacophore-equivalent ' +
    'molecules with different scaffolds score high here even when ECFP4 says ' +
    'they\'re far apart. Cosine values typically run 0.7+ for related chemotypes.');
  desc(COL_MCS_RATIO,
    'Scaffold-hop atom-ratio: ratio_atom = atoms(MCS) / atoms(reference). ' +
    'A candidate is flagged as a hop iff ratio_atom ≤ the user-set threshold ' +
    '(0.4 in the Hard preset, looser in Middle / Easy). Computed only for ' +
    'the top-200 candidates by composite score — rows further down may show NaN.');
  desc(COL_SCORE,
    'Composite ranking score = 0.4 × ECFP4 Tc + 0.6 × CATS2D cosine. CATS gets ' +
    'the heavier weight because it survives scaffold swaps where ECFP4 doesn\'t. ' +
    'When an activity column is selected, the score is also multiplied by ' +
    '(1 - activityWeight) + activityWeight × proximity so candidates near the ' +
    'reference\'s activity rank higher. When CATS fails (Python script ' +
    'unavailable), score falls back to ECFP4 Tc only.');
  desc(COL_FLAG,
    'Scaffold-hop flag. TRUE iff the candidate passes the always-applied ' +
    'atom-ratio criterion (ratio_atom = atoms(MCS) / atoms(reference) ≤ ' +
    'threshold), plus the optionally-enabled Tc-window and CATS-Sim ' +
    'conditions, and (if the user marked atoms) the candidate\'s MCS does ' +
    'NOT preserve the entire marked region. Enable the Reason column in ' +
    'the dialog for the per-row pass/fail breakdown.');
  desc(COL_REASON,
    'Per-row pass/fail breakdown of every flag condition. Format: each ' +
    'condition shows its numeric value and a ✓ or ✗ glyph. "Tc 0.18 ✓ · ' +
    'CATS 0.81 ✓ · MCS 0.52 ✗" reads as: passed Tc-window, passed CATS Sim, ' +
    'failed atom-ratio. Pre-filtered rows show "(pre-filtered)" + the Tc ' +
    'value; the reference shows "(reference)". Filter on this column to ' +
    'find specific failure modes ("show rows where MCS = ✗ but everything ' +
    'else = ✓").');
  if (activityColUsed) {
    desc(COL_ACTIVITY_PROX,
      `Activity-proximity factor in [0, 1] based on the "${activityColUsed}" ` +
      'column. 1.0 = exact match with the reference row\'s activity; values ' +
      'decay as exp(-|delta| / sigma) where sigma is the column\'s standard ' +
      'deviation. When the score is activity-weighted, this is the unblended ' +
      'proximity that went into the multiplier — sort by this column to find ' +
      'hops near the reference\'s potency independent of structure.');
  }
  desc(COL_REPLACEMENT,
    'Full candidate molecule with the shared region highlighted in green. ' +
    'The matching primitive depends on the preset: in Local mode it\'s strict ' +
    'MCS via R-group decomposition (exact-core matching, clean attachment-' +
    'point semantics); in Easy / Middle / Hard modes it\'s ErG (Stiefl 2006 ' +
    'JCIM 46:208) reduced-graph matching, which collapses each ring system to ' +
    'a pharmacophore-labelled node so chemically equivalent heterocycles ' +
    '(pyrimidine ↔ thiazole) match where strict-MCS refuses to. Path-completion ' +
    'extends the highlight along the shortest linker between matched rings so ' +
    'the green region reads as a connected scaffold, not isolated atom dots. ' +
    'Visible only when the user marked atoms on the reference.');
  desc(COL_REPLACED,
    'Standalone fragment SMILES of the candidate\'s region that corresponds ' +
    'to the user\'s marked region on the reference. In Local mode, this is ' +
    'the R-group decomposition output: the candidate\'s R-group fragment at ' +
    'the marked attachment positions (exact-core matching, clean attachment-' +
    'point semantics). In Easy / Middle / Hard modes, this is the ErG-matched ' +
    'region (pharmacophore equivalence, allows ring chemotype swaps). ' +
    'Filterable: use Datagrok\'s substructure filter on this column to find ' +
    '"candidates whose replaced region contains X". Visible only when the ' +
    'user marked atoms on the reference.');
}

/** Removes every Scaffold-Hop result column from `table` (current + legacy ' +
 *  versions) and returns the count removed. Exposed as a top-menu action so
 *  users can undo a run without re-running. Safe to call on a table that
 *  never ran scaffold-hopping (no-op, returns 0). */
export function cleanUpScaffoldHopColumns(table: DG.DataFrame): number {
  const names = [
    COL_TANIMOTO, COL_CATS, COL_MCS_RATIO, COL_SCORE, COL_FLAG,
    COL_REASON, COL_ACTIVITY_PROX, COL_REPLACEMENT, COL_REPLACED,
    COL_PRED_ACTIVITY, COL_PRED_SOURCE, COL_PRED_STDEV, COL_RANK,
    ...LEGACY_PHARM_COLS, ...LEGACY_MCS_COLS, ...LEGACY_SHARED_COLS,
  ];
  let removed = 0;
  for (const name of names) {
    if (table.col(name)) {
      table.columns.remove(name);
      removed++;
    }
  }
  return removed;
}

/** Adds (or replaces) a Molecule-semType column whose cell strings duplicate
 *  the input molecule SMILES, with a per-row substruct provider attached to
 *  `col.temp` so the cell renderer paints highlights on the atoms in the
 *  per-row `ISubstruct`. Same per-row highlight mechanism MMP uses for
 *  from/to pair columns. Used by the `Scaffold Hop Replacement` column to
 *  render the full candidate molecule with the shared region highlighted. */
function replaceMoleculeColumnWithProvider(
  table: DG.DataFrame, molecules: DG.Column, columnName: string,
  byRow: Map<number, ISubstruct>,
  alignByScaffold?: string,
): void {
  if (table.col(columnName)) table.columns.remove(columnName);
  const N = table.rowCount;
  const values: string[] = new Array<string>(N);
  for (let i = 0; i < N; i++) values[i] = molecules.get(i) ?? '';
  const col = DG.Column.fromStrings(columnName, values);
  col.semType = DG.SEMTYPE.MOLECULE;
  table.columns.add(col);
  addSubstructProvider(col.temp, new ScaffoldHopSubstructProvider(byRow, alignByScaffold));
}

/** Per-row substruct provider for the Replacement column. Returns the
 *  shared-region `ISubstruct` for rows where MCS produced one; for other
 *  rows returns an alignment-only substruct (`alignByScaffold` set, no
 *  highlight info) so the cell renderer still re-orients them.
 *
 *  When the optional `alignByScaffold` molblock is supplied, every row
 *  gets that scaffold attached — the cell renderer reads it and calls
 *  `mol.generate_aligned_coords` so the whole column orients to a
 *  common reference. The substruct-provider path explicitly tags
 *  `color: NO_SCAFFOLD_COLOR` and `highlight: false`, so alignment
 *  does NOT trigger the red-highlight side-effect that the column-tag
 *  `.%chem-scaffold-align` path does. */
class ScaffoldHopSubstructProvider implements ISubstructProvider {
  constructor(
    private readonly _byRow: Map<number, ISubstruct>,
    private readonly _alignByScaffold?: string,
  ) {}

  getSubstruct(rowIdx: number | null): ISubstruct | undefined {
    if (rowIdx == null || rowIdx < 0) return undefined;
    const base = this._byRow.get(rowIdx);
    if (this._alignByScaffold) {
      // Compose: keep the per-row highlight (replaced atoms tinted
      // green) and ALSO attach the column-wide alignment scaffold.
      if (base) return {...base, alignByScaffold: this._alignByScaffold};
      // Rows without a highlight still get aligned — `alignByScaffold`
      // alone is enough for the renderer to re-orient the cell.
      return {alignByScaffold: this._alignByScaffold};
    }
    return base;
  }
}

/** Alignment-only substruct provider for the Replaced Region column.
 *  The column's cells are fragment molblocks (with R1 markers); they
 *  don't carry per-row highlights — there's nothing to draw on a
 *  fragment that's already the "thing of interest." All this provider
 *  does is hand the renderer an `alignByScaffold` so every fragment
 *  cell gets oriented to a common anchor (typically a minimal "R1 +
 *  neighbour atom" molblock so the substructure search succeeds on
 *  every piperidine / morpholine / piperazine variant). */
class ScaffoldHopAlignmentOnlyProvider implements ISubstructProvider {
  constructor(private readonly _alignByScaffold: string) {}
  getSubstruct(rowIdx: number | null): ISubstruct | undefined {
    if (rowIdx == null || rowIdx < 0) return undefined;
    return {alignByScaffold: this._alignByScaffold};
  }
}

/** Adds (or replaces) a Molecule-semType column whose cell content is the
 *  per-row fragment representation captured during MCS. Cell values can
 *  be either SMILES (legacy / no-marks path) or V2000 molblocks with
 *  M RGP records (Local-mode path — render attachments as R1, R2, …).
 *  The Chem cell renderer auto-detects the format from the cell content
 *  (multi-line = molblock), so a single column can carry both. Rows
 *  without an entry in `byRow` get an empty cell. */
function replaceFragmentColumn(
  table: DG.DataFrame, columnName: string, byRow: Map<number, string>,
  alignByScaffold?: string,
): void {
  if (table.col(columnName)) table.columns.remove(columnName);
  const N = table.rowCount;
  const values: string[] = new Array<string>(N);
  let anyMolblock = false;
  for (let i = 0; i < N; i++) {
    const v = byRow.get(i) ?? '';
    values[i] = v;
    if (v.length > 0 && v.includes('\n')) anyMolblock = true;
  }
  const col = DG.Column.fromStrings(columnName, values);
  col.semType = DG.SEMTYPE.MOLECULE;
  // Alignment-only substruct provider for the fragment column. Uses the
  // per-row provider mechanism (not the column-tag mechanism) so the
  // renderer sets `color: NO_SCAFFOLD_COLOR` and `highlight: false` on
  // the alignment scaffold — pure alignment, no red highlights on the
  // matched atoms. Skipped when no scaffold supplied (e.g. no marks).
  if (alignByScaffold)
    addSubstructProvider(col.temp, new ScaffoldHopAlignmentOnlyProvider(alignByScaffold));
  // When at least one cell is a molblock (newlines present), tag the
  // column units as MolBlock so the renderer picks up the V2000 M RGP
  // records and shows attachment points as R1/R2 instead of `*:1`.
  // Empty / SMILES-only columns leave the default rendering alone.
  if (anyMolblock) col.meta.units = DG.chem.Notation.MolBlock;
  table.columns.add(col);
}

function replaceFloatColumn(table: DG.DataFrame, name: string, data: Float32Array) {
  if (table.col(name)) table.columns.remove(name);
  table.columns.add(DG.Column.fromFloat32Array(name, data));
}

function replaceBoolColumn(table: DG.DataFrame, name: string, data: Uint8Array) {
  if (table.col(name)) table.columns.remove(name);
  const col = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, name,
    Array.from(data, (x) => x !== 0));
  table.columns.add(col);
}

function replaceStringColumn(table: DG.DataFrame, name: string, data: string[]) {
  if (table.col(name)) table.columns.remove(name);
  table.columns.add(DG.Column.fromStrings(name, data));
}

/** Applies pass/fail color coding using the user-picked thresholds.
 *  Colors are applied to the **text only** (not the cell background) — full
 *  cell shading was visually overwhelming on tables with many result columns.
 *  The `.color-coding-text` column tag (`ddt.api.g.ts:32`) is what flips
 *  Datagrok's grid renderer between text-tint and cell-fill modes. */
function applyColorCoding(
  table: DG.DataFrame, tanimotoMin: number, tanimotoMax: number,
  mcsRatioMax: number, minCatsSim: number,
): void {
  const PASS = DG.Color.fromHtml('#1B9E3F');     // saturated green — readable as foreground
  const MID = DG.Color.fromHtml('#B8860B');      // dark yellow — readable as foreground
  const FAIL = DG.Color.fromHtml('#C92020');     // deep red — readable as foreground
  const fmt = (n: number) => Number.isInteger(n) ? `${n}` : n.toFixed(2);

  const setTextColor = (col: DG.Column | null) => {
    if (col) col.tags['.color-coding-text'] = 'true';
  };

  const tanCol = table.col(COL_TANIMOTO);
  tanCol?.meta.colors.setConditional({
    [`< ${fmt(tanimotoMin)}`]: FAIL,
    [`${fmt(tanimotoMin)}-${fmt(tanimotoMax)}`]: PASS,
    [`> ${fmt(tanimotoMax)}`]: FAIL,
  });
  setTextColor(tanCol);

  const catsCol = table.col(COL_CATS);
  catsCol?.meta.colors.setConditional({
    [`< ${fmt(minCatsSim)}`]: FAIL,
    [`${fmt(minCatsSim)}-1`]: PASS,
  });
  setTextColor(catsCol);

  const mcsCol = table.col(COL_MCS_RATIO);
  mcsCol?.meta.colors.setConditional({
    [`0-${fmt(mcsRatioMax)}`]: PASS,
    [`> ${fmt(mcsRatioMax)}`]: FAIL,
  });
  setTextColor(mcsCol);

  const scoreCol = table.col(COL_SCORE);
  scoreCol?.meta.colors.setLinear([FAIL, MID, PASS]);
  setTextColor(scoreCol);

  // Flag column intentionally has no color coding. The boolean cell renderer
  // ignores the `.color-coding-text` tag and applies categorical colors as
  // cell backgrounds, which dominated the visual once the rest of the table
  // had only text-tinted columns. The checkbox icon already conveys the
  // true/false state on its own.
}

/** Applies the hops-first row order on the active grid (no-op otherwise):
 *
 *    [reference]
 *    [hops by Score desc]
 *    [non-hop survivors by Score desc]
 *    [dropped rows (NaN score) by Tc desc]
 *
 *  Uses `grid.setRowOrder` with an explicit index array — same primitive
 *  `sortBySimilarity` uses (`packages/Chem/src/package.ts:1953`). */
function applyHopsFirstOrder(
  table: DG.DataFrame, referenceRowIdx: number,
  isHop: Uint8Array, score: Float32Array, tanimoto: Float32Array,
): void {
  // Build the explicit display order, then assign integer ranks 0..N-1.
  // Reference = 0; hops sorted by score desc get 1..nHops; non-hop
  // survivors get nHops+1..nHops+nSurv; dropped rows go last sorted by
  // raw Tanimoto desc. Sort the grid ASC on this integer column — clean
  // numeric ranking that doubles as a "this is hop #N" readout for the
  // user.
  const N = table.rowCount;
  const hops: number[] = [];
  const survivors: number[] = [];
  const dropped: number[] = [];
  for (let i = 0; i < N; i++) {
    if (i === referenceRowIdx) continue;
    if (isHop[i]) hops.push(i);
    else if (!Number.isNaN(score[i])) survivors.push(i);
    else dropped.push(i);
  }
  const byScoreDesc = (a: number, b: number) => score[b] - score[a];
  const byTcDesc = (a: number, b: number) => tanimoto[b] - tanimoto[a];
  hops.sort(byScoreDesc);
  survivors.sort(byScoreDesc);
  dropped.sort(byTcDesc);

  const rank = new Int32Array(N);
  rank[referenceRowIdx] = 0;
  let r = 1;
  for (const i of hops) rank[i] = r++;
  for (const i of survivors) rank[i] = r++;
  for (const i of dropped) rank[i] = r++;

  if (table.col(COL_RANK)) table.columns.remove(COL_RANK);
  const rankCol = DG.Column.fromInt32Array(COL_RANK, rank);
  rankCol.setTag('description',
    `Display rank for this scaffold-hopping run. ` +
    `0 = reference row. 1..N = scaffold hops sorted by Scaffold Hop Score ` +
    `descending (so 1 is the top-scoring hop). Then non-hop survivors by ` +
    `Score, then dropped rows by raw Tanimoto.`);
  table.columns.add(rankCol);

  console.log(`[scaffold-hopping] reorder via Rank column: ref=${referenceRowIdx}, ` +
    `hops=${hops.length}, non-hop survivors=${survivors.length}, ` +
    `dropped=${dropped.length}`);

  const tv = grok.shell.tv;
  if (!tv || tv.dataFrame !== table) {
    console.log('[scaffold-hopping] skip grid sort - no active table view for this dataframe');
    return;
  }
  try {
    tv.grid.sort([COL_RANK], [true]);
  } catch (e) {
    console.error('[scaffold-hopping] grid sort failed', e);
  }
}
