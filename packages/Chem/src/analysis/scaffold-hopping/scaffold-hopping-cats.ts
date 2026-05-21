/** CATS2D (Schneider 1999) pharmacophore-pair similarity step for scaffold
 *  hopping, plus the shared `pharmacophore-features.csv` cache that both
 *  CATS and the ErG matching path consume.
 *
 *  Extracted from `scaffold-hopping.ts` as part of the multi-module split;
 *  no behaviour change. The CATS step is a server round-trip into the
 *  Python `CATSFingerprints` script — runtime cost dominates this file's
 *  function call, so the JS-side parsing and cosine math are kept tight. */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as chemCommonRdKit from '../../utils/chem-common-rdkit';
import {removeWaterAndSaltsSingle} from '../../utils/reactions/reactions';

/** Module-level cache for `files/pharmacophore-features.csv`. Same 9 KB CSV
 *  is consumed by two distinct paths (CATS cosine + ErG matching), and a
 *  non-Local run with marked atoms used to fetch it twice per invocation.
 *  Caching the Promise (rather than the DataFrame) lets concurrent callers
 *  share an in-flight request without races.
 *
 *  Rejection handling: a Promise that rejects stays rejected forever, so a
 *  naive `if (!cache) cache = fetch()` would PERMANENTLY disable CATS for
 *  the rest of the session after a single network blip. The `.catch()`
 *  below evicts the failed Promise from the cache so the next caller gets
 *  a fresh fetch attempt, then re-throws so the current caller still sees
 *  the failure (rather than getting a silently undefined DataFrame). */
let _pharmacophoreFeaturesCache: Promise<DG.DataFrame> | null = null;
export function getPharmacophoreFeatures(): Promise<DG.DataFrame> {
  if (!_pharmacophoreFeaturesCache) {
    // Bind to a local first so the catch handler can compare identity
    // against the same reference we'll store in the module cache. (If we
    // chained `.catch()` directly off `loadTable(...)` and stored the
    // result, the local `loadTable(...)` Promise and the stored
    // `.catch()`-wrapped Promise would be different objects and the
    // self-eviction check would never match.)
    const wrapped: Promise<DG.DataFrame> = grok.data.loadTable(
      chemCommonRdKit.getRdKitWebRoot() + 'files/pharmacophore-features.csv',
    ).catch((e) => {
      // Evict ONLY if the cache still points at THIS in-flight Promise.
      // If a later caller already swapped in a new one (race window
      // between this rejection settling and our self-eviction running),
      // leave it alone.
      if (_pharmacophoreFeaturesCache === wrapped) _pharmacophoreFeaturesCache = null;
      throw e;
    });
    _pharmacophoreFeaturesCache = wrapped;
  }
  return _pharmacophoreFeaturesCache;
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
export async function computeCatsCosine(
  molecules: DG.Column, refSmiles: string, survivorIdxs: number[],
  catsCosine: Float32Array,
): Promise<void> {
  // Salt-strip every input. The reference is already pre-stripped by
  // the caller (line where refSmiles is built), but candidates come
  // straight from the column. Without this, counterions (Cl-, Na+, etc.)
  // would contribute extra pharmacophore features and shift cosine
  // similarity by enough to misrank near-tie hops. Same fix is applied
  // on the ECFP4 path below — keeps the three descriptors (Tc / CATS /
  // MCS) operating on identical inputs.
  const inputSmiles = [refSmiles, ...survivorIdxs.map(
    (i) => removeWaterAndSaltsSingle(molecules.get(i) ?? ''))];
  const inputCol = DG.Column.fromStrings('smiles', inputSmiles);
  inputCol.semType = DG.SEMTYPE.MOLECULE;
  const tempDf = DG.DataFrame.fromColumns([inputCol]);

  const featuresDf = await getPharmacophoreFeatures();

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
  if (refVec.length === 0) {
    throw new Error(`CATS2D fingerprint for the reference is empty — ` +
      `RDKit could not parse "${refSmiles}" or no features matched`);
  }
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
export function parseCatsVector(s: string | null | undefined): Float32Array {
  if (!s) return new Float32Array(0);
  const parts = s.split(' ');
  const out = new Float32Array(parts.length);
  for (let i = 0; i < parts.length; i++) {
    const n = parseFloat(parts[i]);
    out[i] = Number.isFinite(n) ? n : 0;
  }
  return out;
}

export function vectorNorm(v: Float32Array): number {
  let s = 0;
  for (let i = 0; i < v.length; i++) s += v[i] * v[i];
  return Math.sqrt(s);
}

/** Cosine similarity over two equal-length float vectors. Returns NaN
 *  when similarity is mathematically undefined (length mismatch, length
 *  zero, or either vector has zero norm — typically a tiny molecule with
 *  no detectable pharmacophore features).
 *
 *  Why NaN and not 0: callers compose this into the composite score
 *  `score = NaN(ca) ? tc : 0.4·tc + 0.6·ca`. A zero-norm CATS returning
 *  0 would compute `0.4·tc + 0.6·0` → artificially penalize feature-less
 *  molecules below their structural similarity. Returning NaN routes the
 *  score through the Tc-only fallback, which is the correct behaviour
 *  ("no CATS signal, fall back to structural").
 *
 *  Pre-computed `aNorm` saves a √ on the reference side across the loop. */
export function cosineSimilarity(a: Float32Array, aNorm: number, b: Float32Array): number {
  if (a.length === 0 || b.length === 0 || a.length !== b.length) return NaN;
  if (aNorm === 0) return NaN;
  let dot = 0;
  let bSq = 0;
  for (let i = 0; i < a.length; i++) {
    dot += a[i] * b[i];
    bSq += b[i] * b[i];
  }
  const bNorm = Math.sqrt(bSq);
  if (bNorm === 0) return NaN;
  return dot / (aNorm * bNorm);
}
