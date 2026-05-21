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

/** Local-mode score weights. Sum to 1. Introduced because the default
 *  `0.4·Tc + 0.6·CATS` formula doesn't penalise candidates that drop
 *  unmarked-side substituents — ECFP4 and CATS2D degrade gracefully when
 *  one ring goes missing, so a "truncated" candidate (e.g. compound 7 in
 *  the GSK650394 series, which keeps the pyrrolopyridine + central
 *  benzamide + cyclopentyl but loses the second phenyl) scores nearly as
 *  high as a fully-preserved candidate whose marked region has been
 *  rewired (the actual intent of a Local hop).
 *
 *  Adding the MCS atom-ratio with a heavy weight (0.5) makes "scaffold
 *  preserved" the DOMINANT scoring criterion for Local mode. Tc and
 *  CATS still contribute (mostly as tie-breakers among similarly-
 *  preserved candidates), but the primary ranking signal is "does this
 *  candidate keep the reference's scaffold?". Combined with the "marked
 *  region changed" check, this produces the ordering the user wants:
 *  candidates that keep the WHOLE scaffold AND rewire the marked region
 *  rank decisively above candidates that drop unmarked substituents.
 *
 *  Weight tuning history: started at 0.3 for a balanced blend. Promoted
 *  to 0.5 after the GSK650394 CAMKK2 case showed that the gap between
 *  truncated (MCS≈0.79) and full-scaffold (MCS≈1.0) candidates needed a
 *  more decisive numerical separation — a 0.5 weight means MCS-difference
 *  of 0.21 contributes 0.105 score difference, large enough to dominate
 *  the small Tc/CATS perturbations between same-family analogs.
 *
 *  Non-Local presets (Easy / Middle / Hard) keep the original
 *  `0.4·Tc + 0.6·CATS` formula intentionally — Maeda 2024's atom-ratio
 *  is already a gate for those presets, and inverting the MCS direction
 *  (`1 - MCS` as bonus) for "real hops" would compete with the gate
 *  rather than complement it. Direction-by-mode added complexity for no
 *  clear win on those presets; trade-off settled in favour of fewer
 *  knobs. */
const W_LOCAL_TANIMOTO = 0.2;
const W_LOCAL_CATS = 0.3;
const W_LOCAL_MCS = 0.5;

/** Maximum MCS atom-ratio value used in the Local score when the marked-
 *  region detector verdict is "changed". Applied as `min(mcsRatio, cap)`,
 *  so candidates whose ratio is BELOW the cap are unaffected.
 *
 *  Rationale: the Maeda atom-ratio metric is composition-only. A
 *  connectivity isomer (same atom inventory as the reference, only the
 *  bonds within the marked region differ — e.g. GSK650394 vs compound
 *  29: same pyrrolopyridine atoms, swapped 2,4 ↔ 3,5 substitution
 *  pattern) saturates at `ratio_atom = 1.0`, because FMCS finds all of
 *  the reference's atoms in the candidate (potentially as disconnected
 *  pieces with mismatched boundary bonds). Bond-Tanimoto would
 *  correctly score this lower (~0.7-0.8), but bond-Tc is a bigger
 *  refactor and an at-rest measure for the next pass.
 *
 *  Capping at 0.95 trims only the very top of the artifact range:
 *    - Connectivity isomer (MCS=1.0, changed) → contribution = 0.5·0.95.
 *    - Same-scaffold + R-swap (MCS≈0.97, preserved) → NOT capped (not
 *      flagged as a hop anyway; `markedOverride === false`).
 *    - Real scaffold hop (MCS<0.95, changed) → NOT capped (below cap).
 *  Keeps the chemically correct ranking (real hop > truncated analog)
 *  while removing the artifactual ceiling that pushed connectivity
 *  isomers to rank #1 by virtue of `ratio_atom = 1.0`. Below ~0.9 the
 *  cap would start inverting the real-hop > truncation ordering — that
 *  trade-off is documented at the call site. */
const W_LOCAL_MCS_CAP_WHEN_CHANGED = 0.95;

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
/** Applicability-domain Tc column. ECFP4 Tanimoto from each predicted
 *  row to its NEAREST measured-activity row in the table (self-excluded).
 *
 *  Why this is the right AD signal for MMP imputation specifically:
 *
 *    1. MMP rules are atom-level substructure transformations learned
 *       on training molecules. For a rule to transfer, the query
 *       molecule must share enough atom-level structure with the
 *       molecules the rule was trained on. ECFP4 (Morgan, radius 2)
 *       captures exactly that — circular atom environments out to
 *       2 bonds — at the same granularity MMP rules operate.
 *    2. Tanimoto is size-invariant. A small fragment with one matching
 *       feature shouldn't be falsely judged "close" to a giant molecule
 *       sharing that one feature; Tanimoto normalises by union.
 *    3. The ECFP4 fingerprints are already computed for the kNN step,
 *       so this column is essentially free (one extra max-Tc pass).
 *
 *  Why not other fingerprints:
 *    - MACCS: too coarse (166 fixed pharmacophore-like bits), misses
 *      the atom-level structural detail MMP rules depend on.
 *    - Atom-pairs / topological torsions: capture longer-range
 *      structure, but largely correlated with ECFP4 on drug-like sets;
 *      not worth the extra fingerprint computation cost.
 *    - 3D / pharmacophore: answers a different question ("similar
 *      binding mode"), not "similar substructure transferability".
 *    - Learned (GNN, ChemBERTa, etc.): could be better but requires
 *      a trained model, validation, and GPU; massive overkill for an
 *      AD signal that's just a "is this row in or out of the training
 *      neighbourhood" check.
 *
 *  Interpretation guidance for users:
 *    - Tc >= 0.6: same chemotype as known compounds → predictions
 *      typically reliable within the model's stress-tested MAE.
 *    - 0.3 <= Tc < 0.6: structurally adjacent but distinct → take
 *      predictions as ballpark estimates with σ-level error bars.
 *    - Tc < 0.3: extrapolation / out-of-domain → predictions are
 *      best-guess only; σ underestimates true uncertainty here.
 *
 *  Reference row gets Tc = 1.0 (it's identical to itself). Rows with
 *  no fingerprint (malformed SMILES) get FNULL. */
const COL_PRED_IN_DOMAIN_TC = 'Scaffold Hop Predicted Activity In-Domain Tc';
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
        refSmiles, new Set(replaceableAtoms), referenceRowIdx, useRGroupReplacement);
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
        // Datagrok column null-sentinels:
        //   FLOAT / DOUBLE → DG.FLOAT_NULL (the JS sentinel; not finite)
        //   INT / BIGINT   → DG.INT_NULL (= -2_147_483_648, a finite int)
        // `Number.isFinite(-2147483648)` returns TRUE, so the old guard
        // `typeof v === 'number' && Number.isFinite(v)` let INT_NULL pass
        // through unchanged. An integer activity column with even one null
        // row therefore silently corrupted MMP's pair deltas, kNN's
        // weighted means, and the proximity factor — all without any
        // visible warning. Branch on column type to use the right sentinel.
        const isIntCol = activityCol.type === DG.COLUMN_TYPE.INT ||
                         activityCol.type === DG.COLUMN_TYPE.BIG_INT;
        const intNullSentinel = DG.INT_NULL;
        for (let i = 0; i < N; i++) {
          const v = activityCol.get(i);
          if (typeof v !== 'number' || !Number.isFinite(v)) {acts[i] = FNULL; continue;}
          if (isIntCol && v === intNullSentinel) {acts[i] = FNULL; continue;}
          acts[i] = v;
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
          // Tc-to-nearest-known-activity-row, per predicted row. Filled
          // by the kNN pass as a free byproduct of the top-K Tanimoto
          // search. Used to populate the applicability-domain column
          // (Scaffold Hop Predicted Activity In-Domain Tc).
          let knnNearestTc: Float32Array | null = null;

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
            // Opt into LOO meanDiff for predict-known rows. Reasoning:
            //
            // Without LOO, the rule's meanDiff includes the target row's
            // own pairs. When we then predict the target row using that
            // rule, we're effectively letting the row contribute to its
            // own prediction (a "leak"). This makes predictions look
            // closer to the row's MEASURED value, but the model is
            // really just fitting the row's own noise rather than
            // learning what other compounds say about it.
            //
            // With LOO, the meanDiff excludes the target row's pairs,
            // so the prediction is what other compounds in the rule
            // collectively predict — the right "validation" estimate
            // when a user has predictKnown=true. Counter-intuitive
            // consequence: on noisy data, LOO predictions look slightly
            // WORSE against measured values, because they predict the
            // signal (noise-free) rather than the signal+noise the
            // assay actually reports. That's the correct behaviour for
            // a validation column; the user asked "what would the model
            // predict if it didn't know this row?" — they didn't ask
            // for the row's own measurement back.
            //
            // The Imatinib stress fixture confirms this: LOO leaves
            // holdout MAE unchanged (0.017), pushes predict-known MAE
            // from 0.035 to 0.042 (signal-vs-noise gap widens), and σ
            // on high-anchor rows becomes a touch larger (more honest).
            // Cost is negligible (O(pairs_in_rule) per predict-known
            // row) and only paid when predictKnown=true.
            const result = imputeMmpActivities(
              mmpa.initData, mmpa.rules, mmpa.rulesBased, mmpa.allCasesBased,
              mmpSupportFloor, 2, predictKnown, predictKnown,
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
          // baseline; k=10 chosen as the conservative end of the k=5-10
          // community range. No single paper establishes 10 specifically —
          // Sheridan 2004 (DOI 10.1021/ci049782w) is the canonical
          // applicability-domain reference (similarity predicts QSAR
          // confidence), which is a sibling concept, not this exact choice.
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
            knnNearestTc = r.nearestKnownTc;
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
              `Estimated spread of the prediction. For MMP rows: ` +
              `decomposes into σ² = σ_between² + σ_within² (between-rule ` +
              `disagreement + average within-rule pair spread), with ` +
              `each rule contributing equally to the between term — so ` +
              `a single dominant rule's tight cluster doesn't ` +
              `artificially collapse σ. Bessel-corrected sample variance ` +
              `(/(n-1)) so the absolute σ value is unbiased rather than ` +
              `systematically understated. For kNN rows: Tc-weighted ` +
              `std-dev of the k neighbour activities. Both in the ` +
              `activity column's units.\n\n` +
              `On the 24-compound Imatinib local-SAR stress fixture σ ` +
              `correlates with the actual |prediction - measured| error ` +
              `at Pearson r ≈ 0.3-0.4 — weakly positive, useful for ` +
              `flagging the highest-σ rows for follow-up but NOT for ` +
              `ranking individual predictions. On larger / more ` +
              `uniformly-easy fixtures (e.g. 60 compounds where every ` +
              `prediction is near the noise floor) r can drop further ` +
              `because |error| is dominated by random measurement noise ` +
              `that σ can't predict. Treat as a directional reliability ` +
              `indicator, not a calibrated uncertainty.\n\n` +
              `For applicability-domain filtering, prefer the ` +
              `companion "In-Domain Tc" column: it answers "is this ` +
              `prediction within the chemical-space coverage of the ` +
              `training set?" — a complementary question that σ ` +
              `cannot answer.`);
            table.columns.add(stdevCol);
          }

          // Applicability-domain Tc column. Emitted whenever kNN ran
          // (regardless of whether kNN actually produced a regression
          // prediction — the Tc-to-nearest is meaningful even on rows
          // where kNN abstained). Skipped when kNN failed entirely
          // (knnNearestTc null) since we have no fingerprint data to
          // populate from.
          if (knnNearestTc !== null) {
            // Reference row's nearest-known-Tc is itself (skipped in the
            // kNN loop via `i === j` exclusion), so it currently sits at
            // FNULL. Pin it to 1.0 — the reference is by definition
            // in-domain for itself, and an "FNULL on the reference"
            // confuses users looking at the column.
            knnNearestTc[referenceRowIdx] = 1.0;
            if (table.col(COL_PRED_IN_DOMAIN_TC))
              table.columns.remove(COL_PRED_IN_DOMAIN_TC);
            const inDomainCol = DG.Column.fromFloat32Array(
              COL_PRED_IN_DOMAIN_TC, knnNearestTc);
            inDomainCol.setTag('description',
              `Applicability-domain signal for the activity prediction. ` +
              `ECFP4 Tanimoto similarity (radius 2, 2048 bits) between ` +
              `this row and its nearest measured-activity neighbour ` +
              `(self-excluded). Higher = the prediction sits within a ` +
              `well-covered area of chemical space, where MMP rules and ` +
              `kNN have anchor support; lower = the prediction is an ` +
              `extrapolation to chemistry unlike anything in the ` +
              `training set.\n\n` +
              `Suggested interpretation:\n` +
              `  • Tc ≥ 0.6 — same chemotype, predictions typically ` +
              `reliable within the stress-tested MAE.\n` +
              `  • 0.3 ≤ Tc < 0.6 — structurally adjacent but distinct; ` +
              `treat predictions as ballpark with σ-level error bars.\n` +
              `  • Tc < 0.3 — extrapolation; predictions are best-guess ` +
              `only, and σ underestimates true uncertainty here.\n\n` +
              `Why ECFP4 + Tanimoto: MMP rules are atom-level ` +
              `substructure transformations, and ECFP4 captures atom ` +
              `environments at the same granularity (radius 2). ` +
              `Tanimoto is size-invariant so analogs of different size ` +
              `are compared fairly. The fingerprints are reused from ` +
              `the kNN step, so this column is essentially free.`);
            table.columns.add(inDomainCol);
          }

          // Fill measured-gaps in `filled[]` from predValue (MMP > kNN).
          let nFilled = 0;
          filled = new Float32Array(N);
          for (let i = 0; i < N; i++) {
            if (acts[i] !== FNULL) filled[i] = acts[i];
            else if (predValue[i] !== FNULL) {filled[i] = predValue[i]; nFilled++;} else filled[i] = FNULL;
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
    const {
      refMcsAtoms: mcsResults, sharedSubstructs, sharedFragments,
      refToCandByRow, candAdjByRow, refAdj,
    } = await computeTopMcsRatio(
      molecules, refSmiles, topMcsRows, mcsRatio, markedRefSet,
      useRGroupReplacement, progress);

    // Pre-compute the marked-region-changed verdict per row. Used by BOTH
    // the Local re-scoring (below) — to cap the MCS contribution when the
    // verdict is "changed" so connectivity isomers can't ride a saturated
    // 1.0 atom-ratio to the top — AND the flag composition loop further
    // down (where it triggers the MCS-gate override). Centralising the
    // computation here avoids re-doing it twice. Empty when no marks.
    const markedRegionChanged = new Set<number>();
    if (replaceableAtoms.length > 0) {
      for (const i of topMcsRows) {
        const refMcsAtoms = mcsResults.get(i);
        if (!refMcsAtoms) continue;
        const compositionPreserved =
          replaceableAtoms.every((a) => refMcsAtoms.has(a));
        let preserved = compositionPreserved;
        if (compositionPreserved) {
          const refToCand = refToCandByRow.get(i);
          const candAdj = candAdjByRow.get(i);
          if (refToCand && candAdj && refAdj.size > 0) {
            preserved = isMarkedRegionTrulyPreserved(
              markedRefSet, refToCand, refAdj, candAdj);
          }
        }
        if (!preserved) markedRegionChanged.add(i);
      }
    }

    // Local-mode MCS-aware re-score. The initial score (0.4·Tc + 0.6·CATS)
    // is used ONLY to pick the top-N for MCS; once MCS is known, candidates
    // that preserve more of the reference scaffold should rank higher in
    // Local mode. Without this, truncated candidates — same marked region
    // but missing an unmarked substituent (e.g. compound 7 in the
    // GSK650394 series) — rank as high as candidates that fully reproduce
    // the scaffold with the marked region rewired (the actual Local-hop
    // intent). Only rows that GOT MCS computed get re-scored; the rest
    // keep the original score and will rank below re-scored survivors.
    //
    // Connectivity-isomer cap: when the marked-region verdict is "changed",
    // the MCS contribution is capped at `W_LOCAL_MCS_CAP_WHEN_CHANGED`
    // (0.95). Strips the artifactual ceiling at `ratio_atom = 1.0` that
    // FMCS produces for connectivity isomers (same atom set, different
    // bonds in the marked region — compound 29 of the GSK650394 series
    // is the canonical example). MCS values below the cap pass through
    // unchanged.
    //
    // Restricted to Local (`useRGroupReplacement=true`). Non-Local presets
    // use Maeda's atom-ratio as a gate, not a continuous signal — adding
    // `(1 - MCS)` as a positive bonus would compete with that gate.
    //
    // Activity-multiplier re-application: the activity-proximity blend
    // (`score *= (1 - aw) + aw * proximity`) was applied to the initial
    // score above. We re-apply it to the MCS-aware score below so the
    // activity-aware re-rank survives the re-scoring step.
    if (useRGroupReplacement) {
      for (const i of topMcsRows) {
        const m = mcsRatio[i];
        if (Number.isNaN(m)) continue; // MCS timed out — keep initial score
        const tc = tanimoto[i];
        const ca = catsCosine[i];
        // Apply cap if marked region changed (likely connectivity-isomer
        // artifact where atom ratio saturates at 1.0 despite real
        // structural difference).
        const mForScore = markedRegionChanged.has(i) ?
          Math.min(m, W_LOCAL_MCS_CAP_WHEN_CHANGED) : m;
        // CATS-NaN fallback: when CATS unavailable, drop the CATS term and
        // proportionally redistribute its weight to Tc + MCS so the score
        // still ranks on real signals (matches the spirit of the original
        // CATS-NaN fallback). Tc and MCS split CATS's weight equally.
        let s: number;
        if (Number.isNaN(ca)) {
          const halfCatsW = W_LOCAL_CATS / 2;
          s = (W_LOCAL_TANIMOTO + halfCatsW) * tc +
              (W_LOCAL_MCS + halfCatsW) * mForScore;
        } else
          s = W_LOCAL_TANIMOTO * tc + W_LOCAL_CATS * ca + W_LOCAL_MCS * mForScore;

        // Re-apply activity-proximity multiplier when it was active.
        const p = activityProx[i];
        if (activityColUsed && activityWeight > 0 && Number.isFinite(p))
          s = s * ((1 - activityWeight) + activityWeight * p);
        score[i] = s;
      }
    }

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
        else {
          reason[i] =
            `Tc ${fmtNum(tc)} ✗ (pre-filtered; outside [${fmtNum(tanimotoMin)}, ${fmtNum(tanimotoMax)}])`;
        }
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

      if (Number.isNaN(ca))
        parts.push('CATS — (Python failure)');
        // CATS-in-flag: NaN degrades to "Tc + Maeda" rather than blocking.
      else {
        const catsOk = ca >= minCatsSim;
        parts.push(`CATS ${fmtNum(ca)} ${catsOk ? '✓' : '✗'}`);
        if (useCatsInFlag && !catsOk) pass = false;
      }

      // MCS gate. Failure tracked separately from `pass` so the marked-
      // region check below can override the gate when the user marked
      // atoms and the connectivity check verdict is "changed" — the
      // marked-region signal is more specific than the global MCS atom-
      // ratio and should win when they disagree. Without the override,
      // the GSK650394 / compound 29 case (same atoms, rewired
      // substituent positions) lands at mcsRatio ≈ 1.0 because FMCS
      // finds all 29 reference atoms in the candidate (potentially as
      // disconnected pieces), and the gate then vetoes a verdict the
      // marked-region detector has already correctly flagged as a hop.
      let mcsFailed = false;
      if (Number.isNaN(m)) {
        parts.push('MCS — (timed out)');
        mcsFailed = true;
      } else {
        const mcsOk = m <= mcsRatioMax;
        parts.push(`MCS ${fmtNum(m)} ${mcsOk ? '✓' : '✗'}`);
        if (!mcsOk) mcsFailed = true;
      }

      // Marked-region check. Three outcomes:
      //  - No marks given: skip; `mcsFailed` directly drives `pass`.
      //  - Marks given AND `preserved`: hard veto (user wanted region
      //    changed; it wasn't, so this isn't the candidate they want).
      //  - Marks given AND `changed`: override MCS gate (mark-region
      //    signal is more specific than global MCS).
      //
      // Verdict comes from the `markedRegionChanged` Set pre-computed
      // above; that set was populated only for rows where the MCS-extract
      // produced `refMcsAtoms`. So `refMcsAtoms == null` is equivalent to
      // "the row is absent from the set" (MCS was skipped on it).
      let markedOverride = false;
      if (replaceableAtoms.length > 0) {
        const refMcsAtoms = mcsResults.get(i);
        if (refMcsAtoms) {
          if (markedRegionChanged.has(i)) {
            // "changed ✓" → override the MCS gate (see comment block
            // above the MCS check). Tc / CATS gates still apply.
            parts.push('marked region changed ✓');
            markedOverride = true;
          } else {
            parts.push('marked region preserved ✗');
            pass = false;
          }
        }
        // If `refMcsAtoms == null`, the "MCS —" token already explains
        // the missing data; we don't add a redundant marked-region tag.
      }

      // Apply MCS gate UNLESS the marked-region check overrode it.
      if (mcsFailed && !markedOverride) pass = false;

      if (pass) isHop[i] = 1;
      reason[i] = parts.join(' · ');
    }

    // Suppress degenerate Replaced Region outputs — when the algorithm
    // finds no real analog (very different chemotypes), it falls back
    // to a small acyclic fragment (e.g. just an NH from FMCS pairing
    // on the marked aniline-N). These mislead the user by suggesting an
    // analog where none exists. Rule: clear the Replaced Region cell
    // when the output has 0 ring atoms AND < 5 heavy atoms AND the
    // ECFP4 Tc to the reference is < 0.2. The combined gate keeps
    // small-but-valid outputs (e.g. an isolated chlorofluoroaniline
    // analog at Tc=0.4 stays visible) while clearing the
    // truly-degenerate ones.
    const rdkitForSuppress = chemCommonRdKit.getRdKitModule();
    for (const idx of [...sharedFragments.keys()]) {
      if (idx === referenceRowIdx) continue;
      const tc = tanimoto[idx];
      if (tc >= 0.2) continue;
      const molblock = sharedFragments.get(idx);
      if (!molblock) continue;
      let mol: any = null;
      try {
        mol = rdkitForSuppress.get_mol(molblock);
        if (!mol) continue;
        const numAtoms = mol.get_num_atoms?.() ?? 0;
        if (numAtoms >= 5) continue;
        const smi = mol.get_smiles?.() ?? '';
        const strippedSmi = smi.replace(/\[\d*\*\]/g, '');
        if (/\d/.test(strippedSmi)) continue; // contains ring closure → keep
        // 0 rings AND < 5 atoms AND Tc < 0.2 → suppress.
        sharedFragments.delete(idx);
        sharedSubstructs.delete(idx);
      } catch {/* parse failure — leave untouched */} finally {mol?.delete?.();}
    }

    progress.update(95, 'Writing result columns...');
    writeOutputColumns(table, molecules, tanimoto, catsCosine, mcsRatio, score, isHop,
      reason, activityProx, activityColUsed,
      sharedSubstructs, sharedFragments, replaceableAtoms.length > 0, showReason,
      refSmiles, markedRefSet, referenceRowIdx, useRGroupReplacement);
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

/** Topology-aware "marked region preserved" check. Returns `true` only if
 *  ALL three of the following hold:
 *    (1) Composition: every marked reference atom pairs to a candidate
 *        atom via the MCS — i.e. the candidate contains the marked atoms.
 *    (2) Boundary degree: each marked ref atom has the same number of
 *        cross-boundary bonds (= bonds to UNMARKED atoms) as its
 *        candidate-side correspondent has bonds to atoms outside the
 *        cand-image of the marked region. Catches the case where the
 *        same marked atoms appear in the candidate but with substituents
 *        at different positions on the marked region (the bicycle).
 *    (3) Matched-neighbour: when the MCS DID extend through a boundary
 *        bond (so we know how the unmarked-side atom corresponds), the
 *        pair must align — the cand correspondent of the boundary
 *        neighbour must be a graph-neighbour of the cand correspondent
 *        of the marked atom.
 *
 *  Why this matters: composition-only (the previous behaviour) miscalls
 *  "preserved" when a candidate has the same marked atoms wired at
 *  different attachment positions on the unmarked core. The canonical
 *  example is compound 29 of the GSK650394 series (Asquith et al.,
 *  J. Med. Chem. 2020, 63, 13750, DOI 10.1021/acs.jmedchem.0c00200):
 *  the pyrrolopyridine bicycle of the reference is "preserved" by
 *  composition (same N, C, and aromatic-ring atom set), but the
 *  benzamide attaches at the 2-position and the phenyl at the
 *  4-position — swapped relative to the reference's 3,5-substitution
 *  pattern. The paper classifies this as a scaffold hop; composition-
 *  only blocked it from being flagged. (3) alone catches some cases but
 *  fails when FMCS doesn't extend into the substituent — it pairs the
 *  bicycle and stops at the boundary, so there's no matched neighbour
 *  to check against. (2) is the load-bearing fix: it doesn't need the
 *  MCS to extend through the boundary; the boundary-degree count alone
 *  is enough to detect attachment-position swaps and N-H ↔ N-Me changes.
 *
 *  Algorithm: O(|marked| × avg-degree) per row.
 *
 *  Safe-default: returns `true` (preserved) on data-missing edge cases —
 *  if the cand mapping or adjacency couldn't be built, we don't flip
 *  the user's call from "preserved" to "changed" on incomplete data. */
function isMarkedRegionTrulyPreserved(
  markedRefAtoms: Set<number>,
  refToCand: Map<number, number>,
  refAdj: Map<number, number[]>,
  candAdj: Map<number, number[]>,
): boolean {
  // (1) Composition: every marked ref atom is paired by the MCS.
  for (const a of markedRefAtoms)
    if (!refToCand.has(a)) return false;
  // Cand-side image of the marked region — used by (2) to compute the
  // candidate's boundary-degree without needing the MCS to extend
  // beyond the marked atoms.
  const markedCandAtoms = new Set<number>();
  for (const a of markedRefAtoms) markedCandAtoms.add(refToCand.get(a)!);
  for (const refMarked of markedRefAtoms) {
    const candCorresp = refToCand.get(refMarked)!;
    const candNbrs = candAdj.get(candCorresp) ?? [];
    // (2) Boundary-degree check. Count bonds out of the marked region on
    // each side. Mismatch = substitution-pattern change (e.g. compound 29).
    let refBoundary = 0;
    for (const nbr of refAdj.get(refMarked) ?? [])
      if (!markedRefAtoms.has(nbr)) refBoundary++;
    let candBoundary = 0;
    for (const nbr of candNbrs)
      if (!markedCandAtoms.has(nbr)) candBoundary++;
    if (refBoundary !== candBoundary) return false;
    // (3) Matched-neighbour check. When the MCS extended through a boundary
    // bond on the ref side, the corresponding bond must exist on the cand
    // side at the paired atoms. Catches different-but-same-degree cases
    // (e.g. bicycle rearrangement where boundary degrees coincidentally
    // match).
    for (const refNbr of refAdj.get(refMarked) ?? []) {
      if (markedRefAtoms.has(refNbr)) continue;
      if (!refToCand.has(refNbr)) continue;
      const expectedCandNbr = refToCand.get(refNbr)!;
      if (!candNbrs.includes(expectedCandNbr)) return false;
    }
  }
  return true;
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
  /** Local-preset flag — when true, the Replacement / Replaced Region
   *  columns were populated via strict MCS + R-group decomposition;
   *  when false, via ErG pharmacophore matching. The column
   *  descriptions are gated on this so the user reads the rules that
   *  actually fired on THIS run rather than a hedged both-mode summary
   *  that mentions ErG even when the actual primitive was strict MCS. */
  useRGroupReplacement: boolean = false,
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
  if (showReason)
    replaceStringColumn(table, COL_REASON, reason);
  else
    if (table.col(COL_REASON)) table.columns.remove(COL_REASON);

  // Activity proximity column is only meaningful when an activity column
  // was selected. Without it the array is all-NaN; drop the column so the
  // table doesn't carry a row of empty cells.
  if (activityColUsed)
    replaceFloatColumn(table, COL_ACTIVITY_PROX, activityProx);
  else
    if (table.col(COL_ACTIVITY_PROX)) table.columns.remove(COL_ACTIVITY_PROX);

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
          if (unmarkedAtoms.size > 0) {
            unmarkedRegionMolblock = buildCleanFragmentMolblock(
              refMolblock, unmarkedAtoms, rdkit);
          }
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
  applyColumnDescriptions(table, activityColUsed, useRGroupReplacement);
}

/** Attaches per-column descriptions so hovering the column header in the
 *  Datagrok grid surfaces a plain-English explanation of what the column
 *  means — the user doesn't have to remember that CATS is Schneider 1999
 *  or that the Maeda atom-ratio is one-sided over the reference. Each
 *  description is written to the `description` tag, the same surface
 *  Datagrok's column tooltip reads. Idempotent — re-running overwrites
 *  with the current text. */
function applyColumnDescriptions(
  table: DG.DataFrame, activityColUsed: string, useRGroupReplacement: boolean,
): void {
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
    'unavailable), score falls back to ECFP4 Tc only.\n\n' +
    'Calibration provenance: the 0.4 / 0.6 split was tuned by inspection ' +
    'on a 5-row BCR-ABL test set (Imatinib + Nilotinib + Dasatinib + ' +
    'Bosutinib + Ponatinib), NOT on a held-out validation set. The ' +
    'direction (CATS-heavier) is consistent with the original CATS ' +
    'methodology (Schneider 1999) and the Maeda 2024 scaffold-hop ' +
    'classifier, but the exact weights have NOT been validated on a ' +
    'large diverse dataset. Treat the Score as a useful ranking signal ' +
    'within a single run, not as an absolute calibrated metric. To ' +
    're-rank by your own priors, sort independently by Tanimoto, CATS ' +
    'Sim, or MCS Ratio columns.');
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
      'hops near the reference\'s potency independent of structure. ' +
      'EXPECTS LOG-SCALED ACTIVITY (pIC50 / pKi / pEC50): raw IC50/Ki/Kd ' +
      'columns are log-normal-distributed and blow sigma out, so the ' +
      'proximity factor collapses to ~1 for every pair (rank no longer ' +
      'shifts). Pre-compute -log10(activity) into a new column and re-run.');
  }
  // Mode-aware descriptions. Pre-fix, both columns emitted a hedged
  // both-mode summary that mentioned ErG even when the actual primitive
  // was strict MCS (in Local mode). Now we describe ONLY the primitive
  // that fired on this run, which is what the user sees in the column.
  if (useRGroupReplacement) {
    desc(COL_REPLACEMENT,
      'Full candidate molecule with the shared region highlighted in green. ' +
      'Local preset: matching is strict MCS via R-group decomposition — ' +
      'the unmarked region of the reference is treated as an EXACT core, ' +
      'and only candidates that contain that core get a populated cell. ' +
      'Same-core analogs (Murcko-equal candidates) get the Murcko scaffold ' +
      'as a placeholder fragment, signalling "scaffold preserved." ' +
      'Visible only when the user marked atoms on the reference.');
    desc(COL_REPLACED,
      'Standalone fragment SMILES of the candidate\'s R-group at the ' +
      'marked attachment positions on the reference. Local-preset R-group ' +
      'decomposition output — exact-core matching with clean attachment-' +
      'point semantics. Empty cell = candidate doesn\'t contain the ' +
      'reference\'s exact core. Filterable: use Datagrok\'s substructure ' +
      'filter on this column to find "candidates whose R-group contains ' +
      'X". Visible only when the user marked atoms on the reference.');
  } else {
    desc(COL_REPLACEMENT,
      'Full candidate molecule with the shared region highlighted in green. ' +
      'Easy / Middle / Hard preset: matching is ErG (Stiefl 2006 JCIM ' +
      '46:208) reduced-graph — each ring system collapses to a ' +
      'pharmacophore-labelled node, so chemically equivalent heterocycles ' +
      '(pyrimidine ↔ thiazole) match where strict-MCS refuses to. ' +
      'Path-completion extends the highlight along the shortest linker ' +
      'between matched rings so the green region reads as a connected ' +
      'scaffold, not isolated atom dots. Visible only when the user ' +
      'marked atoms on the reference.');
    desc(COL_REPLACED,
      'Standalone fragment SMILES of the candidate\'s region that ' +
      'corresponds to the user\'s marked region on the reference. ' +
      'Easy / Middle / Hard preset: ErG-matched region — pharmacophore ' +
      'equivalence allows ring chemotype swaps (pyrimidine ↔ thiazole, ' +
      'phenyl ↔ pyridine). Filterable: use Datagrok\'s substructure ' +
      'filter on this column to find "candidates whose replaced region ' +
      'contains X". Visible only when the user marked atoms on the reference.');
  }
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
  const PASS = DG.Color.fromHtml('#1B9E3F'); // saturated green — readable as foreground
  const MID = DG.Color.fromHtml('#B8860B'); // dark yellow — readable as foreground
  const FAIL = DG.Color.fromHtml('#C92020'); // deep red — readable as foreground
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
