import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as chemSearches from '../../chem-searches';
import * as chemCommonRdKit from '../../utils/chem-common-rdkit';
import {removeWaterAndSaltsSingle} from '../../utils/reactions/reactions';

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

/** Top-N rows by composite score that get the (expensive) per-pair MCS step. */
const TOP_N_FOR_MCS = 200;

/** Per-pair MCS timeout (seconds). Forwarded to RDKit's FMCS via the
 *  `Timeout` option, so the search bounds itself in C++ rather than relying
 *  on a JS-side `Promise.race` (which can't free a synchronous WASM call —
 *  the worker stays stuck and every subsequent pair queues behind it).
 *  8s gives FMCS room to finish on ~40-atom drugs while keeping the
 *  worst-case top-200 sweep bounded at ~27 minutes instead of unbounded. */
const MCS_TIMEOUT_SEC = 8;

/** Atoms must match by element type — `AtomCompare: 'Elements'`. The default
 *  (`'Any'`) is exponential on drug-sized molecules and exceeds the timeout. */
const MCS_EXACT_ATOMS = true;
/** Bonds must match by exact order — `BondCompare: 'OrderExact'`. Prevents
 *  single↔double bond mis-matches in MCS extraction. */
const MCS_EXACT_BONDS = true;

const COL_TANIMOTO = 'Scaffold Hop Tanimoto';
const COL_CATS = 'Scaffold Hop CATS Sim';
const COL_MCS_RATIO = 'Scaffold Hop MCS Ratio';
const COL_SCORE = 'Scaffold Hop Score';
const COL_FLAG = 'Scaffold Hop';

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
 *  7. Output 5 columns + apply pass/fail color coding + apply hops-first row
 *     order (reference → flagged hops by score desc → other survivors by score
 *     desc → dropped rows by Tc desc). */
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

  const progress = DG.TaskBarProgressIndicator.create('Scaffold Hopping...');
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
      writeOutputColumns(table, tanimoto, catsCosine, mcsRatio, score, isHop);
      applyColorCoding(table, tanimotoMin, tanimotoMax, mcsRatioMax, minCatsSim);
      applyHopsFirstOrder(table, referenceRowIdx, isHop, score, tanimoto);
      return;
    }

    progress.update(25, 'Computing CATS2D pharmacophore-pair descriptors...');
    try {
      await computeCatsCosine(molecules, refSmiles, survivorIdxs, catsCosine);
    } catch (e: any) {
      grok.shell.warning(
        `CATS2D pharmacophore computation failed (${e?.message ?? e}). ` +
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

    progress.update(60, `MCS atom-ratio for top ${Math.min(survivorIdxs.length, TOP_N_FOR_MCS)} candidates...`);
    const topByScore = survivorIdxs.slice().sort((a, b) => score[b] - score[a]);
    const topMcsRows = topByScore.slice(0, TOP_N_FOR_MCS);
    const mcsResults = await computeTopMcsRatio(molecules, refSmiles, topMcsRows, mcsRatio, progress);

    // Flag composition. Maeda 2024 atom-ratio is ALWAYS applied (the paper's
    // canonical SH classifier). The two optional conditions — Tc window and
    // CATS Sim — are user-toggleable via `useTcInFlag` / `useCatsInFlag`,
    // letting users drop into "strict-Maeda" mode by turning both off. The
    // marked-atoms refinement is implicit (only active when atoms were
    // marked) and is independent of the toggles.
    for (const i of survivorIdxs) {
      const tc = tanimoto[i];
      const ca = catsCosine[i];
      const m = mcsRatio[i];
      if (Number.isNaN(m)) continue;
      if (m > mcsRatioMax) continue;                                  // Maeda — always applied
      if (useTcInFlag && (tc < tanimotoMin || tc > tanimotoMax)) continue;
      // CATS-in-flag: skip if CATS is NaN (Python failure) — degrades to
      // "Tc + Maeda" rather than blocking the whole run.
      if (useCatsInFlag && !Number.isNaN(ca) && ca < minCatsSim) continue;
      if (replaceableAtoms.length > 0) {
        const refMcsAtoms = mcsResults.get(i);
        if (refMcsAtoms && replaceableAtoms.every((a) => refMcsAtoms.has(a)))
          continue;
      }
      isHop[i] = 1;
    }

    progress.update(95, 'Writing result columns...');
    writeOutputColumns(table, tanimoto, catsCosine, mcsRatio, score, isHop);
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

async function computeTanimoto(
  molecules: DG.Column, refSmiles: string, N: number,
): Promise<Float32Array> {
  const out = new Float32Array(N);
  const tcCol = await chemSearches.chemGetSimilarities(molecules, refSmiles);
  if (!tcCol) {
    grok.shell.error('Failed to compute ECFP4 Tanimoto similarities');
    out.fill(NaN);
    return out;
  }
  for (let i = 0; i < N; i++) out[i] = tcCol.get(i) ?? NaN;
  return out;
}

/** Computes Schneider 1999 CATS2D pharmacophore-pair descriptors for the
 *  reference + every survivor in one server round-trip, then writes the
 *  cosine similarity vs. the reference into `catsCosine` for each survivor.
 *
 *  CATS2D is a *float* vector (not a bit vector) over (familyA, familyB,
 *  topological_distance) triples — 7 families × 7 families × 10 distance
 *  bins = 490 dims — with each (A,B) slice divided by `count(A) + count(B)`
 *  per the Schneider scaling so the descriptor is scale-invariant.
 *
 *  Cosine (not Tanimoto) is the canonical CATS similarity since the vector
 *  carries float-valued, normalised counts rather than independent bits.
 *
 *  Uses the Chem package's 7-family SMARTS (`pharmacophore-features.csv`),
 *  the same source the Pharmacophore Features info panel reads — keeps the
 *  feature definitions consistent with the rest of Chem.
 *
 *  Throws on script-not-registered / Python failure — caller catches and
 *  degrades gracefully (CATS column stays NaN; score falls back to ECFP4
 *  Tc only). */
async function computeCatsCosine(
  molecules: DG.Column, refSmiles: string, survivorIdxs: number[],
  catsCosine: Float32Array,
): Promise<void> {
  const inputSmiles = [refSmiles, ...survivorIdxs.map((i) => molecules.get(i) ?? '')];
  const inputCol = DG.Column.fromStrings('smiles', inputSmiles);
  inputCol.semType = DG.SEMTYPE.MOLECULE;
  const tempDf = DG.DataFrame.fromColumns([inputCol]);

  const featuresDf = await grok.data.loadTable(
    chemCommonRdKit.getRdKitWebRoot() + 'files/pharmacophore-features.csv');

  const catsFunc = DG.Func.find({name: 'CATSFingerprints', package: 'Chem'})[0];
  if (!catsFunc)
    throw new Error('Chem:CATSFingerprints script not registered. Re-run `grok api && grok publish` in packages/Chem.');

  await catsFunc.prepare({
    data: tempDf,
    smiles: inputCol,
    features: featuresDf,
  }).call();

  const fpCol = tempDf.col('cats_fp');
  if (!fpCol)
    throw new Error('CATS2D fingerprint column not produced — Python env may lack rdkit / numpy');

  const refVec = parseCatsVector(fpCol.get(0));
  if (refVec.length === 0)
    throw new Error(`CATS2D fingerprint for the reference is empty — RDKit could not parse "${refSmiles}" or no features matched`);
  const refNorm = vectorNorm(refVec);

  for (let s = 0; s < survivorIdxs.length; s++) {
    const rowIdx = survivorIdxs[s];
    const candVec = parseCatsVector(fpCol.get(s + 1));
    catsCosine[rowIdx] = cosineSimilarity(refVec, refNorm, candVec);
  }
}

/** Parses the space-separated float-vector string produced by
 *  `cats_fingerprints.py` into a Float32Array. Empty / malformed input
 *  yields a zero-length array — caller treats that as "no fingerprint". */
function parseCatsVector(s: string | null | undefined): Float32Array {
  if (!s) return new Float32Array(0);
  const parts = s.split(' ');
  const out = new Float32Array(parts.length);
  for (let i = 0; i < parts.length; i++) {
    const n = parseFloat(parts[i]);
    out[i] = Number.isFinite(n) ? n : 0;
  }
  return out;
}

function vectorNorm(v: Float32Array): number {
  let s = 0;
  for (let i = 0; i < v.length; i++) s += v[i] * v[i];
  return Math.sqrt(s);
}

/** Cosine similarity over two equal-length float vectors. Returns 0 if
 *  either vector has zero norm (no features matched). Pre-computed `aNorm`
 *  saves a √ on the reference side across the loop. */
function cosineSimilarity(a: Float32Array, aNorm: number, b: Float32Array): number {
  if (a.length === 0 || b.length === 0 || a.length !== b.length) return 0;
  if (aNorm === 0) return 0;
  let dot = 0;
  let bSq = 0;
  for (let i = 0; i < a.length; i++) {
    dot += a[i] * b[i];
    bSq += b[i] * b[i];
  }
  const bNorm = Math.sqrt(bSq);
  if (bNorm === 0) return 0;
  return dot / (aNorm * bNorm);
}

/** Computes the Maeda 2024 atom-ratio SH classifier for the top-N rows by
 *  composite score:
 *
 *      ratio_atom(reference, candidate) = |MCS|_atoms / |reference|_atoms
 *
 *  Maeda (J. Chem. Inf. Model. 2024, 64, 5557) defines a candidate as a
 *  scaffold-hopped compound iff `ratio_atom ≤ 0.4` against the query
 *  molecule (which here is the user-selected reference row). The flag is
 *  evaluated downstream in the orchestration's flag step.
 *
 *  Note on naming history: earlier versions of this code used `TcMCS`, the
 *  bond-Tanimoto Maeda also defines (`|MCS|_b / (|A|_b + |B|_b - |MCS|_b)`).
 *  That metric is used in Maeda for chemical-space-network *visualization*
 *  (edges drawn at TcMCS ≥ 0.4 — opposite direction, different formula),
 *  not for SH classification. The conflation has been corrected; `TcMCS`
 *  ≤ 0.4 was never a Maeda criterion.
 *
 *  Returns a map `rowIdx → Set<refAtomIdx>` so the marked-atoms refinement
 *  in the flag step can still test "the user-marked atoms must NOT all be
 *  inside the MCS". Sequential per-pair MCS computation due to the single-
 *  worker WASM bottleneck. */
async function computeTopMcsRatio(
  molecules: DG.Column, refSmiles: string, topRows: number[],
  mcsRatio: Float32Array, progress: DG.TaskBarProgressIndicator,
): Promise<Map<number, Set<number>>> {
  const out = new Map<number, Set<number>>();
  const rdKitService = await chemCommonRdKit.getRdKitService();
  const rdKitModule = chemCommonRdKit.getRdKitModule();

  // Reference heavy-atom count is the SH-ratio denominator — same for every
  // pair, so compute once.
  let refAtomCount = 0;
  let refMolForCount: any = null;
  try {
    refMolForCount = rdKitModule.get_mol(refSmiles);
    refAtomCount = refMolForCount?.get_num_atoms() ?? 0;
  } catch {/* leave 0 — every pair will fail safely */} finally {
    refMolForCount?.delete();
  }
  if (refAtomCount === 0) return out;

  for (let k = 0; k < topRows.length; k++) {
    const rowIdx = topRows[k];
    const candSmilesRaw = molecules.get(rowIdx);
    if (!candSmilesRaw) continue;
    const candSmiles = removeWaterAndSaltsSingle(candSmilesRaw);

    // FMCS bounds itself via the `Timeout` option (forwarded into the C++
    // side through `getMCS(..., timeoutSec)`). On timeout it returns the
    // best partial match, or empty if it found none — both treated as a
    // skip here. No JS-side `Promise.race` because that pattern leaves
    // worker[0] stuck on the previous pair while the next iteration queues.
    let mcsSmarts = '';
    try {
      mcsSmarts = await rdKitService.getMCS(
        [refSmiles, candSmiles], MCS_EXACT_ATOMS, MCS_EXACT_BONDS, MCS_TIMEOUT_SEC);
    } catch {
      continue;
    }
    if (!mcsSmarts) continue;

    let mcsMol: any = null;
    let refMol: any = null;
    try {
      mcsMol = rdKitModule.get_qmol(mcsSmarts);
      const mcsAtoms = mcsMol?.get_num_atoms() ?? 0;
      if (mcsAtoms === 0) continue;

      // Maeda 2024 SH classifier: ratio_atom = MCS_atoms / query_atoms.
      mcsRatio[rowIdx] = mcsAtoms / refAtomCount;

      // Capture which reference atoms are inside the MCS — used by the
      // marked-atoms refinement in the flag step.
      try {
        refMol = rdKitModule.get_mol(refSmiles);
        const matchJson = refMol.get_substruct_match(mcsMol);
        if (matchJson && matchJson !== '{}') {
          const parsed = JSON.parse(matchJson);
          const atoms: number[] = parsed?.atoms ?? [];
          out.set(rowIdx, new Set(atoms));
        }
      } catch {/* ignore */} finally {
        refMol?.delete();
      }
    } finally {
      mcsMol?.delete();
    }

    if ((k % 20) === 0) {
      const pct = 60 + Math.floor(35 * (k / Math.max(1, topRows.length)));
      progress.update(pct, `MCS ${k + 1}/${topRows.length}...`);
    }
  }
  return out;
}

function writeOutputColumns(
  table: DG.DataFrame,
  tanimoto: Float32Array, catsCosine: Float32Array,
  mcsRatio: Float32Array, score: Float32Array, isHop: Uint8Array,
): void {
  // Drop legacy columns from earlier versions so the table doesn't accumulate
  // orphan columns when re-running across upgrades.
  for (const legacy of LEGACY_PHARM_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);
  for (const legacy of LEGACY_MCS_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);

  replaceFloatColumn(table, COL_TANIMOTO, tanimoto);
  replaceFloatColumn(table, COL_CATS, catsCosine);
  replaceFloatColumn(table, COL_MCS_RATIO, mcsRatio);
  replaceFloatColumn(table, COL_SCORE, score);
  replaceBoolColumn(table, COL_FLAG, isHop);
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
  const tv = grok.shell.tv;
  if (!tv || tv.dataFrame !== table) return;

  const hops: number[] = [];
  const nonHopSurvivors: number[] = [];
  const dropped: number[] = [];
  for (let i = 0; i < table.rowCount; i++) {
    if (i === referenceRowIdx) continue;
    if (isHop[i]) hops.push(i);
    else if (!Number.isNaN(score[i])) nonHopSurvivors.push(i);
    else dropped.push(i);
  }
  const byScoreDesc = (a: number, b: number) => score[b] - score[a];
  const byTcDesc = (a: number, b: number) => tanimoto[b] - tanimoto[a];
  hops.sort(byScoreDesc);
  nonHopSurvivors.sort(byScoreDesc);
  dropped.sort(byTcDesc);

  const order = [referenceRowIdx, ...hops, ...nonHopSurvivors, ...dropped];
  try {
    tv.grid.setRowOrder(order);
  } catch {/* grid API mismatch — silent */}
}
