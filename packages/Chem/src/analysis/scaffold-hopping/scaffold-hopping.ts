import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as chemSearches from '../../chem-searches';
import * as chemCommonRdKit from '../../utils/chem-common-rdkit';
import {removeWaterAndSaltsSingle} from '../../utils/reactions/reactions';

/** Composite-score weights. Sum to 1. Three independent low-cost descriptors —
 *  ECFP4 Tanimoto (substructure / "shape" signal), Pharm2D Tanimoto
 *  (Gobbi-Poppinger pairs+triplets bit-vector), and CATS2D cosine
 *  (Schneider 1999 topological pharmacophore-pair float vector with feature-
 *  count normalisation). CATS pulls in low-Tc scaffold-hop candidates that
 *  ECFP4 misses, and the descriptors disagree often enough on hops to give
 *  the composite real signal. */
const W_TANIMOTO = 0.4;
const W_PHARM = 0.3;
const W_CATS = 0.3;

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
const COL_PHARM = 'Scaffold Hop Pharm Sim';
const COL_CATS = 'Scaffold Hop CATS Sim';
/** Legacy column names from earlier iterations. Cleaned up on re-run so users
 *  with old result tables don't end up with multiple pharmacophore columns:
 *  - 'Scaffold Hop Pharm Overlap' (v1: 7-bit family Jaccard)
 *  - 'Scaffold Hop Pharm Cosine' (v1.5 brief: ErG cosine, never shipped) */
const LEGACY_PHARM_COLS = ['Scaffold Hop Pharm Overlap', 'Scaffold Hop Pharm Cosine'];
/** Maeda 2024's actual SH classifier: ratio_atom = atoms(MCS) / atoms(query),
 *  flagged as a hop iff `ratio_atom ≤ 0.4`. (The earlier `TcMCS` column used
 *  the bond-Tanimoto formula and the wrong threshold direction — see the
 *  comment in `computeTopMcsRatio` for the citation correction.) */
const COL_MCS_RATIO = 'Scaffold Hop MCS Ratio';
const LEGACY_MCS_COLS = ['Scaffold Hop TcMCS'];
const COL_SCORE = 'Scaffold Hop Score';
const COL_FLAG = 'Scaffold Hop';

/** Orchestrates the 2D scaffold-hopping pipeline.
 *
 *  1. Standardize SMILES via `removeWaterAndSaltsSingle`.
 *  2. ECFP4 Tanimoto pre-filter — retain rows with `Tc ∈ [tcMin, tcMax]`.
 *  3. Pharmacophore Pharm2D Tanimoto for survivors — RDKit Gobbi-Poppinger
 *     Pharm2D fingerprint built from the Chem package's 7-family SMARTS
 *     with `maxPointCount=3` (pairs and triplets) at distance bins
 *     (0,2)(2,4)(4,6)(6,8)(8,12), Tanimoto over the resulting bit set.
 *  4. CATS2D cosine for survivors — Schneider 1999 topological-pharmacophore-
 *     pair descriptor (7×7×10 = 490-dim float vector with `count(A)+count(B)`
 *     normalisation per pair-and-distance bin). Captures soft scaffold-hop
 *     similarity that ECFP4 misses on isosteric replacements.
 *  5. Composite score = `0.4 × Tc + 0.3 × pharmTanimoto + 0.3 × catsCosine`.
 *  6. MCS atom-ratio top-200 — Maeda 2024 SH classifier:
 *     `ratio_atom = atoms(MCS) / atoms(query)`.
 *  7. Scaffold-hop flag: `ratio_atom ≤ mcsRatioMax (default 0.4 per Maeda) ` +
 *     `AND pharmTanimoto ≥ minPharmOverlap`.
 *     If user marked atoms, additionally require the candidate's MCS does NOT
 *     cover all marked atoms.
 *  8. Output 6 columns + apply pass/fail color coding + apply hops-first row
 *     order (reference → flagged hops by score desc → other survivors by score
 *     desc → dropped rows by Tc desc). */
export async function runScaffoldHopping(
  table: DG.DataFrame,
  molecules: DG.Column,
  referenceRowIdx: number,
  tanimotoMin: number,
  tanimotoMax: number,
  mcsRatioMax: number,
  minPharmOverlap: number,
  replaceableAtomsJson: string,
  useTcInFlag: boolean = true,
  usePharmInFlag: boolean = true,
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
    const pharmJaccard = new Float32Array(N).fill(NaN);
    const catsCosine = new Float32Array(N).fill(NaN);
    const mcsRatio = new Float32Array(N).fill(NaN);
    const score = new Float32Array(N).fill(NaN);
    const isHop = new Uint8Array(N);

    // Reference-row pinning. Hoisted BEFORE the no-survivors early-return so
    // the reference always shows 1.0 across every metric — without this the
    // narrow / inverted-Tc-range case left the reference row's Pharm / CATS
    // / Score / MCS as NaN, which looked like a per-metric failure rather
    // than "the reference is by definition perfectly similar to itself".
    score[referenceRowIdx] = 1.0;
    tanimoto[referenceRowIdx] = 1.0;
    pharmJaccard[referenceRowIdx] = 1.0;
    catsCosine[referenceRowIdx] = 1.0;
    mcsRatio[referenceRowIdx] = 1.0;

    if (survivorIdxs.length === 0) {
      grok.shell.warning(
        `No candidates passed the Tanimoto pre-filter [${tanimotoMin}, ${tanimotoMax}]. ` +
        `Try widening the bounds.`);
      writeOutputColumns(table, tanimoto, pharmJaccard, catsCosine, mcsRatio, score, isHop);
      applyColorCoding(table, tanimotoMin, tanimotoMax, mcsRatioMax, minPharmOverlap);
      applyHopsFirstOrder(table, referenceRowIdx, isHop, score, tanimoto);
      return;
    }

    progress.update(20, 'Computing Pharm2D pharmacophore fingerprints...');
    try {
      await computePharm2DTanimoto(molecules, refSmiles, survivorIdxs, pharmJaccard);
    } catch (e: any) {
      grok.shell.warning(
        `Pharm2D pharmacophore computation failed (${e?.message ?? e}). ` +
        `Pharmacophore similarity will be NaN; flag will rely on MCS atom-ratio only. ` +
        `Check that the Chem Python environment has rdkit available.`);
    }

    progress.update(35, 'Computing CATS2D pharmacophore-pair descriptors...');
    try {
      await computeCatsCosine(molecules, refSmiles, survivorIdxs, catsCosine);
    } catch (e: any) {
      grok.shell.warning(
        `CATS2D pharmacophore computation failed (${e?.message ?? e}). ` +
        `CATS similarity will be NaN; score will fall back to a 50/50 Tc + Pharm blend. ` +
        `Check that the Chem Python environment has rdkit and numpy available.`);
    }

    progress.update(50, 'Computing composite score...');
    // CATS-aware blend: when CATS is NaN (Python failure) fall back to the
    // earlier 0.5 / 0.5 Tc + Pharm composite so the run still produces useful
    // ranking. Ditto for Pharm: if Pharm is NaN we lean on Tc + CATS only.
    for (const i of survivorIdxs) {
      const tc = tanimoto[i];
      const pj = pharmJaccard[i];
      const ca = catsCosine[i];
      const havePharm = !Number.isNaN(pj);
      const haveCats = !Number.isNaN(ca);
      if (havePharm && haveCats) score[i] = W_TANIMOTO * tc + W_PHARM * pj + W_CATS * ca;
      else if (havePharm) score[i] = 0.5 * tc + 0.5 * pj;
      else if (haveCats) score[i] = 0.5 * tc + 0.5 * ca;
      else score[i] = tc;
    }

    progress.update(60, `MCS atom-ratio for top ${Math.min(survivorIdxs.length, TOP_N_FOR_MCS)} candidates...`);
    const topByScore = survivorIdxs.slice().sort((a, b) => score[b] - score[a]);
    const topMcsRows = topByScore.slice(0, TOP_N_FOR_MCS);
    const mcsResults = await computeTopMcsRatio(molecules, refSmiles, topMcsRows, mcsRatio, progress);

    // Flag composition. Maeda 2024 atom-ratio is ALWAYS applied (the paper's
    // canonical SH classifier). The two optional conditions — Tc window and
    // Pharm Sim — are user-toggleable via `useTcInFlag` / `usePharmInFlag`,
    // letting users drop into "strict-Maeda" mode by turning both off. The
    // marked-atoms refinement is implicit (only active when atoms were
    // marked) and is independent of the toggles.
    for (const i of survivorIdxs) {
      const tc = tanimoto[i];
      const pj = pharmJaccard[i];
      const m = mcsRatio[i];
      if (Number.isNaN(m)) continue;
      if (m > mcsRatioMax) continue;                                  // Maeda — always applied
      if (useTcInFlag && (tc < tanimotoMin || tc > tanimotoMax)) continue;
      if (usePharmInFlag && pj < minPharmOverlap) continue;
      if (replaceableAtoms.length > 0) {
        const refMcsAtoms = mcsResults.get(i);
        if (refMcsAtoms && replaceableAtoms.every((a) => refMcsAtoms.has(a)))
          continue;
      }
      isHop[i] = 1;
    }

    // (Reference pinning happens earlier — see the hoisted block above the
    // no-survivors early-return.)

    progress.update(95, 'Writing result columns...');
    writeOutputColumns(table, tanimoto, pharmJaccard, catsCosine, mcsRatio, score, isHop);
    applyColorCoding(table, tanimotoMin, tanimotoMax, mcsRatioMax, minPharmOverlap);
    applyHopsFirstOrder(table, referenceRowIdx, isHop, score, tanimoto);

    const flaggedCount = (isHop as any).reduce((s: number, x: number) => s + x, 0);
    grok.shell.info(
      `Scaffold Hopping done.\n` +
      `Reference: row ${referenceRowIdx} (${refSmiles})\n` +
      `Pre-filter survivors: ${survivorIdxs.length} / ${N - 1}\n` +
      `Top-200 MCS atom-ratio: ${mcsResults.size} pairs computed (${topMcsRows.length - mcsResults.size} timed out)\n` +
      `Flagged as scaffold hops: ${flaggedCount}\n` +
      `Replaceable region: ${replaceableAtoms.length === 0 ? 'auto (no atom marks)' :
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

/** Computes Pharm2D (Gobbi-Poppinger 1998) topological-pharmacophore
 *  fingerprints for the reference + every survivor in one server round-trip,
 *  then writes the Tanimoto similarity vs. the reference into `pharmTanimoto`
 *  for each survivor row.
 *
 *  Uses RDKit's `Chem.Pharm2D.SigFactory` driven by a feature factory built
 *  *from our 7-family SMARTS* (`pharmacophore-features.csv`) — same SMARTS
 *  the Pharmacophore Features info panel uses, including Halogen Bond. With
 *  `maxPointCount=3` the fingerprint encodes both feature *pairs* and
 *  *triplets* at binned topological distances, so it captures roughly where
 *  features sit relative to each other on the molecular graph.
 *
 *  Tanimoto over the resulting SparseBitVect is the standard Pharm2D
 *  similarity. The Python script `Chem:Pharm2DFingerprints` returns the
 *  on-bit indices as a comma-separated string per molecule, which we parse
 *  into `Set<number>` here for fast intersection / union counts.
 *
 *  Throws on script-not-registered / Python failure — caller catches and
 *  degrades gracefully (pharm column stays NaN; flag relies on TcMCS only). */
async function computePharm2DTanimoto(
  molecules: DG.Column, refSmiles: string, survivorIdxs: number[],
  pharmTanimoto: Float32Array,
): Promise<void> {
  const inputSmiles = [refSmiles, ...survivorIdxs.map((i) => molecules.get(i) ?? '')];
  const inputCol = DG.Column.fromStrings('smiles', inputSmiles);
  inputCol.semType = DG.SEMTYPE.MOLECULE;
  const tempDf = DG.DataFrame.fromColumns([inputCol]);

  // Load our 7-family SMARTS — same source as the Pharmacophore Features
  // info panel, so this stays consistent with the rest of Chem.
  const featuresDf = await grok.data.loadTable(
    chemCommonRdKit.getRdKitWebRoot() + 'files/pharmacophore-features.csv');

  const pharm2dFunc = DG.Func.find({name: 'Pharm2DFingerprints', package: 'Chem'})[0];
  if (!pharm2dFunc)
    throw new Error('Chem:Pharm2DFingerprints script not registered. Re-run `grok api && grok publish` in packages/Chem.');

  await pharm2dFunc.prepare({
    data: tempDf,
    smiles: inputCol,
    features: featuresDf,
  }).call();

  const fpCol = tempDf.col('pharm2d_fp');
  if (!fpCol)
    throw new Error('Pharm2D fingerprint column not produced — Python env may lack rdkit.Chem.Pharm2D');

  const refOnBits = parsePharm2DOnBits(fpCol.get(0));
  if (refOnBits.size === 0)
    throw new Error(`Pharm2D fingerprint for the reference is empty — RDKit could not parse "${refSmiles}" or no features matched`);

  for (let s = 0; s < survivorIdxs.length; s++) {
    const rowIdx = survivorIdxs[s];
    const candOnBits = parsePharm2DOnBits(fpCol.get(s + 1));
    pharmTanimoto[rowIdx] = bitTanimoto(refOnBits, candOnBits);
  }
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
 *  Same 7-family SMARTS as the Pharm2D step, so the two descriptors share
 *  feature definitions — the CATS view emphasises pairwise distance
 *  patterns, the Pharm2D view emphasises pair- and triplet-bit overlap.
 *  Throws on script-not-registered / Python failure — caller catches and
 *  degrades gracefully (CATS column stays NaN; score falls back to 0.5/0.5
 *  Tc + Pharm). */
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

/** Parses the comma-separated on-bit-indices string produced by
 *  `pharm2d_fingerprints.py` into a Set<number>. Empty string / malformed
 *  input yields an empty set — caller treats that as "no fingerprint". */
function parsePharm2DOnBits(s: string | null | undefined): Set<number> {
  if (!s) return new Set();
  const out = new Set<number>();
  for (const part of s.split(',')) {
    const n = parseInt(part, 10);
    if (Number.isFinite(n)) out.add(n);
  }
  return out;
}

/** Tanimoto similarity over two on-bit sets: |A ∩ B| / |A ∪ B|. Returns 0
 *  if either set is empty — degrades safely for the downstream flag check. */
function bitTanimoto(a: Set<number>, b: Set<number>): number {
  if (a.size === 0 || b.size === 0) return 0;
  let intersection = 0;
  const [smaller, larger] = a.size < b.size ? [a, b] : [b, a];
  for (const x of smaller) if (larger.has(x)) intersection++;
  const union = a.size + b.size - intersection;
  return union === 0 ? 0 : intersection / union;
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
  tanimoto: Float32Array, pharmJaccard: Float32Array, catsCosine: Float32Array,
  mcsRatio: Float32Array, score: Float32Array, isHop: Uint8Array,
): void {
  // Drop legacy columns from earlier versions so the table doesn't accumulate
  // orphan columns when re-running across upgrades.
  for (const legacy of LEGACY_PHARM_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);
  for (const legacy of LEGACY_MCS_COLS)
    if (table.col(legacy)) table.columns.remove(legacy);

  replaceFloatColumn(table, COL_TANIMOTO, tanimoto);
  replaceFloatColumn(table, COL_PHARM, pharmJaccard);
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
  mcsRatioMax: number, minPharmOverlap: number,
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

  const pharmCol = table.col(COL_PHARM);
  pharmCol?.meta.colors.setConditional({
    [`< ${fmt(minPharmOverlap)}`]: FAIL,
    [`${fmt(minPharmOverlap)}-1`]: PASS,
  });
  setTextColor(pharmCol);

  // CATS isn't gated by a user threshold — it only contributes to the
  // composite score — so a smooth gradient communicates "more is better"
  // without implying a hard pass/fail line.
  const catsCol = table.col(COL_CATS);
  catsCol?.meta.colors.setLinear([FAIL, MID, PASS]);
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

  const flagCol = table.col(COL_FLAG);
  flagCol?.meta.colors.setCategorical({'true': PASS, 'false': FAIL});
  setTextColor(flagCol);
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
