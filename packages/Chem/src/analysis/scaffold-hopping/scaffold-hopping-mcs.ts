/** ECFP4 Tanimoto pre-filter + FMCS sweep with Murcko-equality skip — the
 *  expensive structural-similarity steps of the scaffold-hopping pipeline.
 *
 *  Extracted from `scaffold-hopping.ts` as part of the multi-module split;
 *  no behaviour change. Constants `TOP_N_FOR_MCS` and friends moved with
 *  their consumer. */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import * as chemSearches from '../../chem-searches';
import * as chemCommonRdKit from '../../utils/chem-common-rdkit';
import {removeWaterAndSaltsSingle} from '../../utils/reactions/reactions';

import {
  buildFamilyQmols, disposeFamilyQmols, computeErgSharedAtoms, FamilyQmols,
  labelAtomsByPharmacophore,
} from './erg-scaffold-match';
import {
  extractFragmentSmiles, parseV2000Adjacency, extractFragmentMolblockWithRGroups,
  findRingSystems, selectReplacementRingAtoms,
} from './scaffold-hopping-molblock';
import {getPharmacophoreFeatures} from './scaffold-hopping-cats';

/** Top-N rows by composite score that get the (expensive) per-pair MCS step. */
export const TOP_N_FOR_MCS = 200;

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

/** Green used to colour the shared region (atoms IN the MCS) when rendered
 *  on the full candidate molecule in the `Scaffold Hop Replacement` column. */
const SHARED_HIGHLIGHT_COLOR = '#5BC85B';

export async function computeTanimoto(
  molecules: DG.Column, refSmiles: string, N: number,
): Promise<Float32Array> {
  const out = new Float32Array(N);
  // Salt-strip every candidate before fingerprinting — keeps ECFP4 / CATS /
  // MCS operating on identical inputs (counterions otherwise inflate both
  // bit counts on the ECFP side and pharmacophore counts on the CATS side,
  // silently shifting similarity scores). The reference is already
  // stripped by the caller (refSmiles).
  const stripped: string[] = new Array(N);
  for (let i = 0; i < N; i++)
    stripped[i] = removeWaterAndSaltsSingle(molecules.get(i) ?? '');
  const strippedCol = DG.Column.fromStrings(molecules.name, stripped);
  strippedCol.semType = DG.SEMTYPE.MOLECULE;
  const tcCol = await chemSearches.chemGetSimilarities(strippedCol, refSmiles);
  if (!tcCol) {
    grok.shell.error('Failed to compute ECFP4 Tanimoto similarities');
    out.fill(NaN);
    return out;
  }
  for (let i = 0; i < N; i++) out[i] = tcCol.get(i) ?? NaN;
  return out;
}

// CATS step (computeCatsCosine, parseCatsVector, vectorNorm,
// cosineSimilarity) moved to ./scaffold-hopping-cats.ts

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
 *  Returns three per-row maps:
 *  - `refMcsAtoms`: `rowIdx → Set<refAtomIdx>` — which reference atoms are
 *    inside the strict MCS, used by the marked-atoms refinement in the
 *    flag step (Maeda's classifier needs this).
 *  - `sharedSubstructs`: `rowIdx → ISubstruct` — candidate-side atoms in
 *    the **ErG-matched shared region** (NOT the strict MCS), packaged
 *    for the `Scaffold Hop Replacement` column's per-row green highlight.
 *    ErG (Stiefl 2006) collapses each ring system to a pharmacophore-
 *    labelled node and matches by label compatibility, so pyrimidine ↔
 *    pyrimidine, pyridyl ↔ thiazole, tolyl ↔ chloroaryl all match — the
 *    chemically meaningful "this is what plays the same role" view.
 *  - `sharedFragments`: `rowIdx → string` — same ErG-matched atoms
 *    extracted as a standalone fragment SMILES for the
 *    `Scaffold Hop Shared Region` column.
 *
 *  Note: the SH flag (Maeda atom-ratio) still uses the *strict* MCS, so
 *  the paper-faithful classifier behaviour is preserved. ErG only drives
 *  the per-row Replacement / Shared Region visualisations.
 *
 *  Sequential per-pair MCS computation due to the single-worker WASM
 *  bottleneck. */
export async function computeTopMcsRatio(
  molecules: DG.Column, refSmiles: string, topRows: number[],
  mcsRatio: Float32Array, markedRefAtoms: Set<number>,
  useRGroup: boolean,
  progress: DG.TaskBarProgressIndicator,
): Promise<{
  refMcsAtoms: Map<number, Set<number>>;
  sharedSubstructs: Map<number, ISubstruct>;
  sharedFragments: Map<number, string>;
  /** Per-row ref-atom-index → cand-atom-index mapping from the strict MCS.
   *  Used by the topology-aware "marked region preserved" check in the flag
   *  step: composition (atom-in-MCS) is necessary but not sufficient —
   *  connectivity must also be preserved. Without this map the check
   *  miscalls "preserved" cases where the same marked atoms are wired to
   *  the unmarked region at different attachment positions (e.g. compound
   *  30 of the GSK650394 series: same pyrrolopyridine atoms, swapped
   *  2,4 ↔ 3,5 substitution pattern). Populated only when atoms are
   *  marked (`markedRefAtoms.size > 0`). */
  refToCandByRow: Map<number, Map<number, number>>;
  /** Per-row candidate-atom adjacency map (0-indexed). Sister to
   *  `refToCandByRow` — same connectivity check needs both ref and cand
   *  graphs to compare neighbour identities. Populated only when atoms
   *  are marked. */
  candAdjByRow: Map<number, Map<number, number[]>>;
  /** Reference-atom adjacency map (0-indexed). Built once outside the
   *  per-row loop since the reference is fixed across all candidates.
   *  Populated only when atoms are marked. Empty map when no marks. */
  refAdj: Map<number, number[]>;
}> {
  const out = new Map<number, Set<number>>();
  const sharedSubstructs = new Map<number, ISubstruct>();
  const sharedFragments = new Map<number, string>();
  const refToCandByRow = new Map<number, Map<number, number>>();
  const candAdjByRow = new Map<number, Map<number, number[]>>();
  const rdKitService = await chemCommonRdKit.getRdKitService();
  const rdKitModule = chemCommonRdKit.getRdKitModule();
  const sharedColorRgb = DG.Color.hexToPercentRgb(SHARED_HIGHLIGHT_COLOR);

  // Replacement-column computation mode. Two paths:
  // - `useRGroup = true` (Local preset): per-pair STRICT MCS extraction.
  //   For each candidate, the MCS-matched atoms on the candidate side
  //   are the "preserved" region; everything else is the "replaced"
  //   region. No ErG ring-system matching, no chain-atom matching, no
  //   path-completion — just MCS atoms in vs. MCS atoms out. Tolerates
  //   slight variations in the kept region (MCS finds the largest
  //   common subgraph, not exact-core identity, so the Aurora-A-style
  //   case where candidates have slightly different propionamide /
  //   pyrazole substituents still produces clean output).
  // - `useRGroup = false` (Easy / Middle / Hard): per-row ErG
  //   pharmacophore matching. Pharmacophore-equivalent ring matches
  //   (pyrimidine ↔ thiazole etc.) handled, at the cost of fuzzier
  //   boundaries.
  // The strict-MCS atom-ratio computation that drives the Maeda flag
  // is unaffected — it runs in both modes (same MCS call is reused).
  let familyQmols: FamilyQmols | null = null;
  if (markedRefAtoms.size > 0 && !useRGroup) {
    try {
      const featuresDf = await getPharmacophoreFeatures();
      familyQmols = buildFamilyQmols(featuresDf, rdKitModule);
    } catch {/* ErG matching falls back to no-op if features can't load */}
  }

  // Reference molblock — used by the Local-mode path to compute ref
  // adjacency for "edge anchor" detection AND by the topology-aware
  // "marked region preserved" check across all modes when atoms are
  // marked. Fetched once outside the per-row loop since the reference
  // doesn't change between rows.
  let refMolblockForLocal = '';
  if (markedRefAtoms.size > 0) {
    let tmpRefMol: any = null;
    try {
      tmpRefMol = rdKitModule.get_mol(refSmiles);
      refMolblockForLocal = tmpRefMol?.get_molblock?.() ?? '';
    } catch {/* leave empty — Local-mode falls back to largest-component */} finally {tmpRefMol?.delete();}
  }
  // Reference adjacency — same precondition (marks present); cheap to
  // parse once and reused for every candidate's connectivity comparison
  // AND for the Local-mode edge-anchor detection inside the loop.
  const refAdjMap: Map<number, number[]> = markedRefAtoms.size > 0 ?
    parseV2000Adjacency(refMolblockForLocal) : new Map();

  // Pharmacophore-based ring matching — pre-compute the marked ref ring
  // CATS-like signatures (count of atoms per pharmacophore family) once
  // here, so the per-row substruct fallback can pair cand rings to ref
  // rings by pharmacophore similarity (not just by structural match).
  // This catches the case where a cand has a ring playing the
  // "marked-region pharmacophore role" but structurally resembles a
  // non-marked ref ring — e.g. sorafenib's pyridyl plays the
  // heteroaromatic-hinge-binder role of imatinib's MARKED pyrimidine
  // (2 N's), but structurally it's a 1-N pyridine that pure-substruct
  // matches imatinib's NON-MARKED pyridyl decoration. CATS pairs by
  // pharmacophore-feature-similarity, so pyridyl → pyrimidine (both
  // heteroaromatic with HBA N) wins over pyridyl → benzene.
  const markedRefRingSigs: Map<string, number>[] = [];
  const markedRefRings: Set<number>[] = [];
  if (familyQmols && markedRefAtoms.size > 0 && refMolblockForLocal) {
    let refMolForLabels: any = null;
    try {
      refMolForLabels = rdKitModule.get_mol(refSmiles);
      if (refMolForLabels) {
        const refLabels = labelAtomsByPharmacophore(refMolForLabels, familyQmols);
        const allRefRings = findRingSystems(refAdjMap);
        for (const ring of allRefRings) {
          let overlapsMarked = false;
          for (const a of ring)
            if (markedRefAtoms.has(a)) {overlapsMarked = true; break;}

          if (!overlapsMarked) continue;
          const sig = new Map<string, number>();
          for (const a of ring) {
            const fams = refLabels.get(a) ?? new Set<string>();
            for (const f of fams) sig.set(f, (sig.get(f) ?? 0) + 1);
          }
          markedRefRingSigs.push(sig);
          markedRefRings.push(ring);
        }
      }
    } catch {/* CATS pre-compute failed — falls back to structural-only */} finally {refMolForLabels?.delete?.();}
  }

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
  if (refAtomCount === 0) {
    return {
      refMcsAtoms: out, sharedSubstructs, sharedFragments,
      refToCandByRow, candAdjByRow, refAdj: refAdjMap,
    };
  }

  // Murcko-equality fast path. Batch-compute Bemis-Murcko scaffold SMILES
  // for the reference + every top-row candidate via `Chem:MurckoScaffolds`
  // (one server round-trip on ~200 molecules, sub-second on typical
  // datasets). For any candidate whose Murcko matches the reference's
  // exactly, the two molecules share their entire ring + linker
  // topology — the only difference is R-group decoration. Skip FMCS and
  // mark mcsRatio = 1.0 (full coverage by construction), which
  // correctly classifies them as not-a-hop without paying the 8 s
  // worst-case FMCS cost. This:
  //   - cuts wall-clock for SAR-style series where many candidates share
  //     the reference's scaffold (~50-80 % of FMCS calls skipped),
  //   - eliminates the FMCS-timeout failure mode on those same pairs
  //     (timeout → partial MCS → falsely low atom-ratio → false hop
  //     flag), which is the more important correctness benefit on
  //     complex / symmetric molecules.
  //
  // Gated on Local preset (`useRGroup`). Empirically on Easy/Middle/Hard
  // the upper Tc bound (≤ 0.50/0.70) already pre-filters out
  // same-scaffold candidates — they all have Tc > 0.6 by construction —
  // so the Murcko batch returns 0 matches and the ~500 ms call buys
  // nothing. Local's Tc window (0.50-0.95) is exactly where
  // same-scaffold candidates live, so the optimisation pays off
  // (50-80 % skip rate measured on quinazoline / kinase-inhibitor SAR
  // datasets).
  //
  // Skipped candidates don't populate sharedSubstructs / sharedFragments —
  // those columns are only meaningful when the user marked atoms, and a
  // Murcko-equal pair has no marked-region distinction to highlight.
  const murckoByRow = new Map<number, string>();
  let refMurcko = '';
  if (!useRGroup) {
    console.log('[scaffold-hopping] Murcko skip disabled on non-Local presets ' +
      '(upper Tc bound pre-filters same-scaffold candidates).');
  } else if (markedRefAtoms.size > 0) {
    // When atoms are marked, the marked-region-changed detector needs the
    // strict-MCS atom mapping (refMcsAtoms / refToCand / candAdj) on every
    // candidate. The Murcko-skip path produces NONE of these (it sets
    // `mcsRatio = 1.0` and `continue`s before the MCS extraction block),
    // so marked-region info would be missing on same-scaffold candidates
    // — exactly the population where local-hop detection matters (a
    // candidate sharing the reference's Bemis-Murcko scaffold but with
    // substituents rewired between marked-region atoms IS a local hop;
    // pre-mark code couldn't distinguish that case and treated all same-
    // scaffold candidates as "obviously not hops" via the 1.0 ceiling).
    // Trading ~50-80% of FMCS skip rate for correctness on the marked-
    // atoms population — Local-with-marks workflows are exactly where
    // users want the strictest analysis anyway.
    console.log('[scaffold-hopping] Murcko skip disabled when marked atoms ' +
      'are given (marked-region-changed check needs full FMCS).');
  } else {
    try {
      const murckoInputSmiles = [refSmiles];
      const murckoInputRowIdxs: number[] = [];
      for (const rowIdx of topRows) {
        const candRaw = molecules.get(rowIdx);
        if (!candRaw) continue;
        murckoInputSmiles.push(removeWaterAndSaltsSingle(candRaw));
        murckoInputRowIdxs.push(rowIdx);
      }
      const inputCol = DG.Column.fromStrings('smiles', murckoInputSmiles);
      inputCol.semType = DG.SEMTYPE.MOLECULE;
      const inputDf = DG.DataFrame.fromColumns([inputCol]);
      const murckoFunc = DG.Func.find({name: 'MurckoScaffolds', package: 'Chem'})[0];
      if (!murckoFunc)
        console.warn('[scaffold-hopping] Chem:MurckoScaffolds not found, FMCS skip disabled.');
      if (murckoFunc) {
      // `MurckoScaffolds` declares `action:join(data)` on its output, so the
      // scaffolds get joined into the INPUT dataframe (`inputDf`) as a new
      // `scaffolds` column — they are NOT in `call.getOutputParamValue()`,
      // which returns an empty join-action sentinel dataframe.
        const t0 = performance.now();
        await murckoFunc.prepare({data: inputDf, smiles: inputCol}).call();
        const t1 = performance.now();
        const scaffoldsCol = inputDf.col('scaffolds');
        if (scaffoldsCol) {
          refMurcko = scaffoldsCol.get(0) ?? '';
          for (let i = 0; i < murckoInputRowIdxs.length; i++) {
            const m = scaffoldsCol.get(i + 1);
            if (m) murckoByRow.set(murckoInputRowIdxs[i], m);
          }
          const willSkip = refMurcko ? Array.from(murckoByRow.values())
            .filter((m) => m === refMurcko).length : 0;
          console.log(`[scaffold-hopping] Murcko batch: ${Math.round(t1 - t0)} ms, ` +
          `ref scaffold "${refMurcko.slice(0, 60)}${refMurcko.length > 60 ? '…' : ''}", ` +
          `${willSkip}/${murckoByRow.size} candidates will skip FMCS.`);
        }
      }
    } catch (e) {
      console.warn('[scaffold-hopping] Murcko batch failed, falling back to FMCS-only path:', e);
    // Empty refMurcko / empty murckoByRow means the equality check below
    // never short-circuits — every pair goes through FMCS as before.
    }
  }

  for (let k = 0; k < topRows.length; k++) {
    // Cooperative cancel — each FMCS call can take up to 8 s, so an
    // uncancelable sweep over 200 candidates is the worst-case wall clock
    // (~27 min). Checking before each pair gives the user an out without
    // adding measurable overhead.
    if (progress.canceled) throw new Error('Scaffold hopping cancelled by user.');
    const rowIdx = topRows[k];
    const candSmilesRaw = molecules.get(rowIdx);
    if (!candSmilesRaw) continue;
    const candSmiles = removeWaterAndSaltsSingle(candSmilesRaw);

    // Murcko-equality short-circuit. Empty Murcko (acyclic molecules)
    // never short-circuits: '' === '' would over-match every fully-
    // acyclic pair as "same scaffold," which is wrong — acyclic
    // molecules need real FMCS to decide. The `refMurcko` non-empty
    // guard catches that.
    const candMurcko = murckoByRow.get(rowIdx);
    if (refMurcko && candMurcko && candMurcko === refMurcko) {
      mcsRatio[rowIdx] = 1.0;
      // No need to set out (refMcsAtoms) — it's only consumed by the
      // marked-atoms flag refinement, and a Murcko-equal pair fails the
      // hop criterion at the atom-ratio gate (1.0 > any preset's cap),
      // so the marked-atoms branch never runs for this row.
      //
      // BUT in Local mode (preset=Local) the Replacement / Replaced
      // Region columns get populated for every survivor, regardless of
      // whether the candidate flags as a hop — Local users explicitly
      // care about "what's preserved" on close analogs. The Murcko
      // short-circuit used to leave those cells blank because it
      // bypassed the FMCS extraction step that fills sharedFragments.
      // Emit the conserved Murcko scaffold as a placeholder fragment so
      // the column communicates "this candidate is a same-scaffold
      // analog of the reference (substituent-level differences not
      // extracted here)" rather than appearing empty. Non-Local users
      // still see the placeholder; it's never misleading.
      sharedFragments.set(rowIdx, refMurcko);
      continue;
    }

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

      // Capture which reference atoms are inside the strict MCS. The
      // marked-atoms refinement in the flag step needs this (it tests
      // "did the candidate's MCS cover ALL marked atoms" against Maeda's
      // strict atom-by-atom MCS, NOT the looser ErG match).
      //
      // When atoms are marked, ALSO build the ref→cand atom mapping and
      // the candidate adjacency map for the topology-aware connectivity
      // check that downstream code uses to distinguish true preservation
      // (same atoms, same wiring) from substitution-pattern swaps (same
      // atoms, different wiring — e.g. compound 30 in the GSK650394
      // series: same pyrrolopyridine, swapped 2,4 ↔ 3,5 substitution).
      let candMolForMapping: any = null;
      try {
        refMol = rdKitModule.get_mol(refSmiles);
        const refMatchJson = refMol.get_substruct_match(mcsMol);
        let refMatchAtoms: number[] = [];
        if (refMatchJson && refMatchJson !== '{}') {
          const parsed = JSON.parse(refMatchJson);
          refMatchAtoms = parsed?.atoms ?? [];
          out.set(rowIdx, new Set(refMatchAtoms));
        }
        if (markedRefAtoms.size > 0 && refMatchAtoms.length > 0) {
          // Use the already-salt-stripped `candSmiles` (line 354), not the
          // raw `molecules.get(rowIdx)`. mcsMol was built from FMCS over
          // the salt-stripped pair, so substruct-matching against an
          // unstripped molecule would return atom indices that include
          // salt-fragment positions — wrong highlights, wrong refToCand
          // mapping, wrong adjacency for the connectivity check.
          if (candSmiles) {
            candMolForMapping = rdKitModule.get_mol(candSmiles);
            if (candMolForMapping) {
              const candMatchJson = candMolForMapping.get_substruct_match(mcsMol);
              if (candMatchJson && candMatchJson !== '{}') {
                const candMatchAtoms: number[] = JSON.parse(candMatchJson)?.atoms ?? [];
                const pairLen = Math.min(refMatchAtoms.length, candMatchAtoms.length);
                if (pairLen > 0) {
                  const refToCand = new Map<number, number>();
                  for (let i = 0; i < pairLen; i++)
                    refToCand.set(refMatchAtoms[i], candMatchAtoms[i]);
                  refToCandByRow.set(rowIdx, refToCand);
                  // Cand adjacency, parsed from the candidate's molblock.
                  const candMolblock = candMolForMapping.get_molblock?.() ?? '';
                  if (candMolblock)
                    candAdjByRow.set(rowIdx, parseV2000Adjacency(candMolblock));
                }
              }
            }
          }
        }
      } catch {/* ignore */} finally {
        refMol?.delete();
        candMolForMapping?.delete();
      }

      // Replacement-column data for this row.
      if (useRGroup) {
        // Local-mode path: MCS-mapping-based candidate extraction.
        //
        // The MCS computed above (mcsMol) maps reference atoms to
        // candidate atoms pair-wise: mcsMol atom k ↔ ref atom
        // refMatch[k] ↔ cand atom candMatch[k]. The "preserved" region
        // on the candidate side is the set of candidate atoms that pair
        // with NON-MARKED reference atoms — i.e. atoms the user said
        // they want to KEEP. The "replaced" region is everything else:
        // candidate atoms that pair with a marked ref atom (the MCS
        // happened to find a cross-match that doesn't matter for the
        // user's intent), or candidate atoms not in the MCS at all
        // (genuinely new chemistry).
        //
        // This is sharper than "candidate atoms not in MCS" because it
        // respects the user's mark/unmark boundary even when the MCS
        // over-matches across ring/non-ring types.
        let candMolForLocal: any = null;
        let refMolForLocal: any = null;
        try {
          // Reuse the salt-stripped `candSmiles` so substruct-match atom
          // indices align with the MCS (which was computed on the
          // stripped pair). Using raw `molecules.get(rowIdx)` here used
          // to produce shifted highlight atoms on salt-bearing molecules.
          if (candSmiles && sharedColorRgb) {
            candMolForLocal = rdKitModule.get_mol(candSmiles);
            refMolForLocal = rdKitModule.get_mol(refSmiles);
            if (candMolForLocal && refMolForLocal) {
              // Per-side atom mappings. Each `atoms[i]` is the side's
              // atom matching mcsMol's atom i, so refAtoms[i] ↔
              // candAtoms[i] is a pair.
              let refMcsAtomsArr: number[] = [];
              let candMcsAtomsArr: number[] = [];
              try {
                const refJson = refMolForLocal.get_substruct_match(mcsMol);
                if (refJson && refJson !== '{}')
                  refMcsAtomsArr = JSON.parse(refJson)?.atoms ?? [];
              } catch {/* leave empty */}
              try {
                const candJson = candMolForLocal.get_substruct_match(mcsMol);
                if (candJson && candJson !== '{}')
                  candMcsAtomsArr = JSON.parse(candJson)?.atoms ?? [];
              } catch {/* leave empty */}

              // Candidate atoms paired to *non-marked* ref atoms =
              // "preserved" (the user wants these kept).
              const candPreserved = new Set<number>();
              const pairLen = Math.min(refMcsAtomsArr.length, candMcsAtomsArr.length);
              if (markedRefAtoms.size > 0) {
                for (let i = 0; i < pairLen; i++) {
                  if (!markedRefAtoms.has(refMcsAtomsArr[i]))
                    candPreserved.add(candMcsAtomsArr[i]);
                }
              } else {
                // No marks → fall back to the simpler definition
                // "every MCS-paired candidate atom is preserved." The
                // replaced region is then just whatever's outside the
                // MCS on the candidate side.
                for (let i = 0; i < pairLen; i++)
                  candPreserved.add(candMcsAtomsArr[i]);
              }

              const candNumAtoms = candMolForLocal.get_num_atoms?.() ?? 0;
              const candMolblock = candMolForLocal.get_molblock?.() ?? '';

              // All candidate atoms NOT paired to a non-marked ref atom.
              // Includes (a) atoms paired to MARKED ref atoms — direct
              // anchors of the replacement region — and (b) candidate
              // atoms with no MCS pair, which can be EITHER the genuine
              // replacement (e.g., a piperidine swapped for the marked
              // piperazine where exact-atom MCS rejected the cross-ring
              // match) OR noise from MCS imperfections elsewhere in the
              // molecule (e.g., an acetamide-tail mismatch on the OTHER
              // side of the molecule from the marked region). To
              // distinguish, we score connected components by whether
              // they touch a marked-region anchor.
              const allNonPreserved = new Set<number>();
              for (let a = 0; a < candNumAtoms; a++)
                if (!candPreserved.has(a)) allNonPreserved.add(a);

              // Full candidate adjacency (every atom, not just
              // non-preserved) — we need it to check whether components
              // are adjacent to anchor atoms in the preserved region.
              const candAdj = new Map<number, number[]>();
              for (let a = 0; a < candNumAtoms; a++) candAdj.set(a, []);
              let parsedNA = 0;
              let parsedNB = 0;
              if (candMolblock) {
                const lines = candMolblock.split('\n');
                const counts = lines[3];
                if (counts) {
                  parsedNA = parseInt(counts.substring(0, 3));
                  parsedNB = parseInt(counts.substring(3, 6));
                  if (Number.isFinite(parsedNA) && Number.isFinite(parsedNB)) {
                    for (let j = 0; j < parsedNB; j++) {
                      const bondLine = lines[4 + parsedNA + j];
                      if (!bondLine || bondLine.length < 6) continue;
                      const a = parseInt(bondLine.substring(0, 3)) - 1;
                      const b = parseInt(bondLine.substring(3, 6)) - 1;
                      candAdj.get(a)?.push(b);
                      candAdj.get(b)?.push(a);
                    }
                  }
                }
              }

              // Direct anchors: candidate atoms paired (via MCS) to
              // MARKED ref atoms. These ARE the replacement region when
              // MCS happens to find cross-ring atom matches.
              const directAnchors = new Set<number>();
              for (let i = 0; i < pairLen; i++) {
                if (markedRefAtoms.has(refMcsAtomsArr[i]))
                  directAnchors.add(candMcsAtomsArr[i]);
              }

              // Edge anchors: candidate atoms paired to UNMARKED ref
              // atoms that are immediate neighbours (in ref graph) of
              // marked ref atoms. These sit at the boundary of the
              // marked region on the candidate side. Used when the MCS
              // doesn't pair across the ring-type mismatch (e.g.,
              // piperazine ↔ piperidine, where exact-atom matching
              // refuses the N↔C swap), so directAnchors would be empty
              // but the replacement region is still topologically next
              // to the preserved core.
              // Reuses the hoisted `refAdjMap` (built once outside the
              // per-row loop) so we don't re-parse the reference molblock
              // ~200 times.
              const refBorderAtoms = new Set<number>();
              for (const m of markedRefAtoms) {
                for (const nbr of refAdjMap.get(m) ?? [])
                  if (!markedRefAtoms.has(nbr)) refBorderAtoms.add(nbr);
              }
              const edgeAnchors = new Set<number>();
              for (let i = 0; i < pairLen; i++) {
                if (refBorderAtoms.has(refMcsAtomsArr[i]))
                  edgeAnchors.add(candMcsAtomsArr[i]);
              }

              // Ring-substitution selection (primary path). The user's
              // intent when marking a ring in the reference is "find me
              // candidates whose RING at that position is different" —
              // not "find me candidates whose entire half of the molecule
              // is different". The previous BFS-on-non-preserved-atoms
              // approach unioned everything not in the MCS-preserved core
              // and produced too-large Replaced Region cells on Middle/
              // global hops (e.g. DLK pyrimidine → pyrazole-with-
              // piperidinyl: pyrazole IS the ring swap; piperidinyl is a
              // separate substituent change).
              //
              // Selection criterion: a ring system in the candidate is a
              // "replacement" iff (a) it isn't entirely inside the
              // MCS-preserved set (= it isn't part of the conserved
              // core), AND (b) at least one of its atoms has a neighbour
              // OUTSIDE the ring that is an EDGE anchor (= a cand atom
              // paired to a ref atom that bordered the marked region in
              // ref).
              //
              // Why "bonded to edge anchor outside the ring" and not
              // "contains a direct anchor": the FMCS atom-by-atom
              // pairing is brittle. It can place cand piperidine atoms
              // into `directAnchors` (paired to ref pyrimidine atoms by
              // element only) even when the chemically correct pairing
              // would be cand-piperidine ↔ ref-piperidine. The
              // direct-anchor check then wrongly selects the cand
              // piperidine ring as a "replacement". Edge-anchor
              // adjacency is more robust because (i) edge anchors are
              // paired to ref atoms IMMEDIATELY OUTSIDE the marked
              // region (= the substituent atoms bordering the marked
              // ring in ref), so spoofing them requires the candidate
              // to actually be wired to the chemically equivalent
              // boundary position, and (ii) the "outside the ring"
              // constraint ensures we only count true cross-ring bonds,
              // not intra-ring edges.
              //
              // Why "bonded to ANY candPreserved" is also wrong: the
              // oxetane on cand 19's piperidine is bonded to the
              // piperidine N (in candPreserved), but the piperidine N
              // in ref is NOT adjacent to the marked pyrimidine (it's
              // the far end of piperidine, holding the acetyl group).
              // So the piperidine N is NOT an edgeAnchor — its match
              // sits two bonds inside the conserved core, not at the
              // boundary. The edge-anchor test correctly excludes
              // oxetane.
              //
              // Falls back to the legacy BFS-on-non-preserved approach
              // when no ring system qualifies (acyclic marked region,
              // or no rings near the boundary).
              const ringSystems = findRingSystems(candAdj);
              const selectedRingSystems: Set<number>[] = [];
              for (const ring of ringSystems) {
                let allPreserved = true;
                for (const a of ring)
                  if (!candPreserved.has(a)) {allPreserved = false; break;}

                if (allPreserved) continue;

                let bondedToEdge = false;
                for (const a of ring) {
                  for (const nbr of candAdj.get(a) ?? []) {
                    if (ring.has(nbr)) continue;
                    if (edgeAnchors.has(nbr)) {bondedToEdge = true; break;}
                  }
                  if (bondedToEdge) break;
                }
                if (bondedToEdge) selectedRingSystems.push(ring);
              }

              const chosen = new Set<number>();
              if (selectedRingSystems.length > 0) {
                // Ring-substitution path. Union of selected ring systems,
                // restricted to non-preserved atoms (exclude any ring
                // atom that's actually in the MCS-preserved set, which
                // can happen at the boundary of a fused ring system).
                for (const ring of selectedRingSystems) {
                  for (const a of ring)
                    if (!candPreserved.has(a)) chosen.add(a);
                }
              }

              // Fallback path — original BFS-on-non-preserved-atoms
              // logic, used only when ring-substitution picked NOTHING
              // (acyclic marked region case, no rings at the anchor
              // position, etc.). Same algorithm as before: BFS connected
              // components, include those touching an anchor.
              if (chosen.size === 0) {
                const components: Set<number>[] = [];
                const visited = new Set<number>();
                for (const start of allNonPreserved) {
                  if (visited.has(start)) continue;
                  const component = new Set<number>();
                  const stack: number[] = [start];
                  while (stack.length > 0) {
                    const x = stack.pop()!;
                    if (visited.has(x)) continue;
                    visited.add(x);
                    component.add(x);
                    for (const y of candAdj.get(x) ?? [])
                      if (!visited.has(y) && allNonPreserved.has(y)) stack.push(y);
                  }
                  components.push(component);
                }
                const includedComponents: Set<number>[] = [];
                for (const c of components) {
                  let touches = false;
                  for (const a of c) {
                    if (directAnchors.has(a)) {touches = true; break;}
                    let adjacentEdge = false;
                    for (const nbr of candAdj.get(a) ?? [])
                      if (edgeAnchors.has(nbr)) {adjacentEdge = true; break;}

                    if (adjacentEdge) {touches = true; break;}
                  }
                  if (touches) includedComponents.push(c);
                }
                if (includedComponents.length > 0)
                  for (const c of includedComponents) for (const a of c) chosen.add(a);
                else if (components.length > 0) {
                  // Last-resort: largest component, so the cell isn't empty.
                  let largest = components[0];
                  for (const c of components) if (c.size > largest.size) largest = c;
                  for (const a of largest) chosen.add(a);
                }
              }
              const replacedAtoms = Array.from(chosen);
              if (replacedAtoms.length > 0) {
                const replacedSet = chosen;
                // Bonds inside the chosen component(s) for highlighting.
                const replacedBonds: number[] = [];
                if (candMolblock && Number.isFinite(parsedNA) && Number.isFinite(parsedNB)) {
                  const lines = candMolblock.split('\n');
                  for (let j = 0; j < parsedNB; j++) {
                    const bondLine = lines[4 + parsedNA + j];
                    if (!bondLine || bondLine.length < 6) continue;
                    const a = parseInt(bondLine.substring(0, 3)) - 1;
                    const b = parseInt(bondLine.substring(3, 6)) - 1;
                    if (replacedSet.has(a) && replacedSet.has(b))
                      replacedBonds.push(j);
                  }
                }

                const highlightAtomColors: {[k: number]: number[]} = {};
                for (const a of replacedAtoms)
                  highlightAtomColors[a] = [...sharedColorRgb];
                const highlightBondColors: {[k: number]: number[]} = {};
                for (const b of replacedBonds)
                  highlightBondColors[b] = [...sharedColorRgb];
                sharedSubstructs.set(rowIdx, {
                  atoms: replacedAtoms,
                  bonds: replacedBonds,
                  highlightAtomColors,
                  highlightBondColors,
                });

                if (candMolblock) {
                  // Local mode uses the with-RGroups molblock extractor:
                  // attachment points render as `R1`, `R2` (V2000
                  // M RGP records the grid cell renderer picks up
                  // automatically). A single-atom replacement renders
                  // as e.g. `C-R1` instead of "CH4".
                  const frag = extractFragmentMolblockWithRGroups(
                    candMolblock, replacedSet, rdKitModule);
                  if (frag) sharedFragments.set(rowIdx, frag);
                }
              }
            }
          }
        } catch {/* Local-mode extraction failure — leave the row empty */} finally {
          candMolForLocal?.delete();
          refMolForLocal?.delete();
        }
      } else {
        // Non-Local (Easy / Middle / Hard) path. When atoms are marked,
        // FIRST try the same MCS-based ring-substitution selection the
        // Local path uses — it directly answers "which cand ring sits
        // at the position of the marked ref ring" using the strict-MCS
        // atom pairing (`refToCandByRow`) instead of pharmacophore-
        // family matching. This gives a chemically meaningful answer
        // for the DLK pyrimidine → pyridine / pyrazole / thiazole
        // / indazole series. The ErG path only fires as a fallback,
        // when no ring system qualifies — typically because the marked
        // region was acyclic, or the MCS didn't extend through the
        // boundary on the cand side.
        let candMolForErg: any = null;
        let refMolForErg: any = null;
        try {
          // Reuse the salt-stripped `candSmiles` for the same reason as
          // the Local path above: ErG matches against an atom-indexed
          // structure whose indices must align with the MCS-derived
          // anchors we computed earlier on the stripped pair.
          if (!candSmiles || !sharedColorRgb) {/* skip */} else {
            candMolForErg = rdKitModule.get_mol(candSmiles);
            refMolForErg = rdKitModule.get_mol(refSmiles);
            if (candMolForErg && refMolForErg) {
              let replacedAtoms = new Set<number>();
              let replacedBonds: number[] = [];
              const candMolblock = candMolForErg.get_molblock?.() ?? '';

              // Primary: MCS-based ring-substitution selection.
              const refToCand = refToCandByRow.get(rowIdx);
              const candAdj = candAdjByRow.get(rowIdx);
              if (markedRefAtoms.size > 0 && refToCand && candAdj &&
                  refAdjMap.size > 0) {
                // Build candPreserved (cand atoms paired to non-marked ref).
                const candPreserved = new Set<number>();
                for (const [refA, candA] of refToCand.entries())
                  if (!markedRefAtoms.has(refA)) candPreserved.add(candA);

                // Build edge anchors (cand atoms paired to ref atoms
                // immediately adjacent to the marked region).
                const refBorderAtoms = new Set<number>();
                for (const m of markedRefAtoms) {
                  for (const nbr of refAdjMap.get(m) ?? [])
                    if (!markedRefAtoms.has(nbr)) refBorderAtoms.add(nbr);
                }
                const edgeAnchors = new Set<number>();
                for (const [refA, candA] of refToCand.entries())
                  if (refBorderAtoms.has(refA)) edgeAnchors.add(candA);

                replacedAtoms = selectReplacementRingAtoms(
                  candAdj, candPreserved, edgeAnchors);

                // Second fallback within the marked-mode path: when the
                // global MCS optimised for OTHER parts of the molecule
                // and missed the marked region entirely (= no edge
                // anchors), substructure-match the marked-region atoms
                // directly against the candidate. This recovers cases
                // like Imatinib vs Dasatinib where FMCS pairs the
                // amide/piperazine arms (~11 atoms) and SKIPS the
                // pyrimidine entirely — leaving the user's marked
                // region with no MCS image.
                //
                // Approach: extract the marked-region sub-SMILES from
                // the reference, get_qmol it, substruct-match in cand,
                // and use the matched cand atoms as directAnchors. Then
                // include the ring system containing directAnchors PLUS
                // any ring system within 2 bonds of directAnchors that
                // isn't in candPreserved (e.g. dasatinib's thiazole,
                // separated from pyrimidine by one NH linker).
                // Always run the substruct-match fallback when atoms are
                // marked. The MCS-based ring-substitution above can
                // produce small / partial results when FMCS spuriously
                // pairs cand atoms to ref atoms outside the marked
                // region (e.g. dasatinib's thiazole carbons pair-up
                // with imatinib's pyridine carbons by element, dropping
                // them into candPreserved and shrinking the output).
                // Substruct-match doesn't depend on MCS — it directly
                // searches for the marked-region pattern in the
                // candidate. We UNION the two results so the output
                // contains everything either path finds.
                if (markedRefAtoms.size > 0) {
                  try {
                    const refMolblock = refMolForErg.get_molblock?.() ?? '';
                    // Build the set of "marked ref atoms unioned with
                    // any candidate ring" candidates. The marked region
                    // may span MULTIPLE ring systems (e.g. user marks
                    // pyrimidine + tolyl phenyl). Decompose the marked
                    // region per-ring-system so a candidate that has
                    // only ONE of the marked rings (e.g. dasatinib has
                    // pyrimidine but no tolyl) still scores some
                    // directAnchors. Plus try the FULL marked-region
                    // SMILES first as the strongest signal.
                    const directAnchors = new Set<number>();
                    // Pre-compute conservedByMatch: cand atoms covered by
                    // any NON-marked ref ring substruct. Used (a) to
                    // reject anchor candidates that land entirely on
                    // conserved atoms (e.g. an extracted phenyl ring of
                    // the marked aminotolyl matches bosutinib's
                    // chlorobenzene, but chlorobenzene is also matched
                    // by imatinib's benzamide phenyl — it's CONSERVED,
                    // not novel) and (b) to seed the no-anchor fallback
                    // when no novel anchor survives.
                    const conservedByMatchPre = new Set<number>();
                    {
                      const refRingsPre = findRingSystems(refAdjMap);
                      for (const refRing of refRingsPre) {
                        let overlapsMarked = false;
                        for (const a of refRing)
                          if (markedRefAtoms.has(a)) {overlapsMarked = true; break;}

                        if (overlapsMarked) continue;
                        const refRingFrag = extractFragmentSmiles(
                          refMolblock, refRing, rdKitModule);
                        if (!refRingFrag) continue;
                        let refRingQmol: any = null;
                        try {
                          refRingQmol = rdKitModule.get_qmol(refRingFrag);
                          if (!refRingQmol) continue;
                          const j = candMolForErg.get_substruct_match(refRingQmol);
                          if (j && j !== '{}') {
                            const ms: number[] = JSON.parse(j)?.atoms ?? [];
                            for (const a of ms) conservedByMatchPre.add(a);
                          }
                        } catch {/* ignore */} finally {refRingQmol?.delete?.();}
                      }
                    }
                    // Reject a candidate match if EVERY matched atom is
                    // conserved. Such matches land on the conserved core
                    // and would produce a spurious "ring substitution"
                    // result (e.g. bosutinib's phenyl ring matched by
                    // imatinib's tolyl-phenyl, but it's the SAME ring
                    // role as the conserved benzamide phenyl).
                    const isNovelMatch = (atoms: number[]): boolean => {
                      for (const a of atoms)
                        if (!conservedByMatchPre.has(a)) return true;
                      return false;
                    };
                    // Primary: full marked-region SMILES.
                    const fullMarkedSub = extractFragmentSmiles(
                      refMolblock, markedRefAtoms, rdKitModule);
                    // `markedRegionPreserved` = the candidate contains the
                    // FULL marked-region SMILES as a substructure. When
                    // this fires (close-analog cases like erlotinib /
                    // lapatinib / afatinib that all share gefitinib's
                    // anilinoquinazoline), we DON'T extend to adjacent
                    // novel rings or run CATS — the user marked region IS
                    // the analog, full stop. Otherwise the output bloats
                    // with substituent rings (lapatinib's furan + benzyl-
                    // fluorobenzyl arm) that weren't in the marked region.
                    let markedRegionPreserved = false;
                    if (fullMarkedSub) {
                      let fullQmol: any = null;
                      try {
                        fullQmol = rdKitModule.get_qmol(fullMarkedSub);
                        if (fullQmol) {
                          const j = candMolForErg.get_substruct_match(fullQmol);
                          if (j && j !== '{}') {
                            const atoms: number[] = JSON.parse(j)?.atoms ?? [];
                            if (isNovelMatch(atoms)) {
                              for (const a of atoms) directAnchors.add(a);
                              markedRegionPreserved = true;
                            }
                          }
                        }
                      } catch {/* ignore */} finally {fullQmol?.delete?.();}
                    }
                    // Per-ring fallback: each ref ring that contains
                    // marked atoms gets its own SMARTS, matched
                    // independently. Multi-ring marked regions (e.g.
                    // imatinib's pyrimidine + tolyl) need this so a
                    // candidate that has only ONE of the marked rings
                    // (dasatinib has pyrimidine but no tolyl) still
                    // scores anchors. The HETEROATOM-COUNT-DESC sort
                    // ensures the MOST SPECIFIC ring (pyrimidine, 2 N)
                    // becomes the "primary" anchor. Later rings are
                    // accepted only if their match in the candidate is
                    // within K=2 graph hops of the primary anchor — so
                    // a tolyl substruct that lands on a chemically
                    // DISTANT phenyl ring in the candidate (e.g.
                    // dasatinib's chloro-methyl phenyl, 4+ bonds from
                    // pyrimidine) is rejected, while a tolyl substruct
                    // INSIDE the reference's own marked region (2 bonds
                    // from pyrimidine via NH) is accepted.
                    if (directAnchors.size === 0) {
                      const refRings = findRingSystems(refAdjMap);
                      const markedRefRings = refRings.filter((r) => {
                        for (const a of r) if (markedRefAtoms.has(a)) return true;
                        return false;
                      });
                      // Heteroatom-count descending so the most
                      // discriminating ring is tried first.
                      const refAtomEl: Map<number, string> = (() => {
                        const m = new Map<number, string>();
                        const lines = refMolblock.split('\n');
                        const counts = lines[3];
                        const nA = parseInt(counts.substring(0, 3));
                        if (Number.isFinite(nA)) {
                          for (let i = 0; i < nA; i++) {
                            const ln = lines[4 + i];
                            if (ln && ln.length >= 34) m.set(i, ln.substring(31, 34).trim());
                          }
                        }
                        return m;
                      })();
                      const heteroCount = (ring: Set<number>) => {
                        let n = 0;
                        for (const a of ring) {
                          const el = refAtomEl.get(a) ?? 'C';
                          if (el !== 'C' && el !== 'H') n++;
                        }
                        return n;
                      };
                      markedRefRings.sort((a, b) => heteroCount(b) - heteroCount(a));
                      // First ring becomes primary; later rings are
                      // distance-gated.
                      const primaryAnchors = new Set<number>();
                      for (let r = 0; r < markedRefRings.length; r++) {
                        const refRing = markedRefRings[r];
                        const ringSubSmi = extractFragmentSmiles(
                          refMolblock, refRing, rdKitModule);
                        if (!ringSubSmi) continue;
                        let ringQmol: any = null;
                        let matched: number[] = [];
                        try {
                          ringQmol = rdKitModule.get_qmol(ringSubSmi);
                          if (!ringQmol) continue;
                          const j = candMolForErg.get_substruct_match(ringQmol);
                          if (j && j !== '{}') matched = JSON.parse(j)?.atoms ?? [];
                        } catch {/* ignore */} finally {ringQmol?.delete?.();}
                        if (matched.length === 0) continue;
                        // Reject matches that land entirely on conserved
                        // atoms (would produce a spurious anchor on the
                        // conserved core rather than a novel ring).
                        if (!isNovelMatch(matched)) continue;
                        if (primaryAnchors.size === 0) {
                          // First successful match becomes primary.
                          for (const a of matched) {
                            primaryAnchors.add(a);
                            directAnchors.add(a);
                          }
                        } else {
                          // Distance check: at least one matched atom
                          // must be within K=2 hops of primaryAnchors in
                          // the cand graph. Otherwise this match lives at
                          // a chemically unrelated position (dasatinib's
                          // chloro-methyl-phenyl when tolyl is matched).
                          const reachable = new Set<number>(primaryAnchors);
                          let frontier = new Set<number>(primaryAnchors);
                          for (let d = 0; d < 2; d++) {
                            const next = new Set<number>();
                            for (const a of frontier) {
                              for (const nbr of candAdj.get(a) ?? []) {
                                if (reachable.has(nbr)) continue;
                                reachable.add(nbr);
                                next.add(nbr);
                              }
                            }
                            frontier = next;
                          }
                          let close = false;
                          for (const a of matched)
                            if (reachable.has(a)) {close = true; break;}

                          if (close) for (const a of matched) directAnchors.add(a);
                        }
                      }
                    }
                    // FMCS pairing: seed additional anchors from FMCS-
                    // paired cand atoms (every cand atom paired with a
                    // marked ref atom via FMCS becomes an anchor). This
                    // catches close analogs where the marked-region
                    // SMILES is too specific to substruct-match (e.g.
                    // nilotinib's tolyl has only ONE ring N vs imatinib's
                    // two — exact SMILES match fails, but FMCS still
                    // pairs the atoms cleanly).
                    //
                    // GRAPH-DISTANCE GATE: when we already have "core"
                    // anchors from substruct match (per-ring or full
                    // marked SMILES), only accept FMCS-paired atoms
                    // within K=3 bonds of those core anchors in the cand
                    // graph. This rejects FMCS pairings that land on
                    // chemically-distant parts of the candidate — e.g.
                    // FMCS pairs imatinib's tolyl (in the marked region)
                    // with dasatinib's chlorotolyl (in the AMIDE ARM,
                    // ~7 bonds from dasatinib's pyrimidine), but those
                    // belong to a different pharmacophore position and
                    // would bloat the output. When no core anchors exist
                    // yet, accept all FMCS pairings.
                    if (refToCand) {
                      const coreAnchors = new Set(directAnchors);
                      let nearCore: Set<number> | null = null;
                      if (coreAnchors.size > 0) {
                        nearCore = new Set(coreAnchors);
                        let frontier = new Set(coreAnchors);
                        for (let d = 0; d < 3; d++) {
                          const next = new Set<number>();
                          for (const a of frontier) {
                            for (const nbr of candAdj.get(a) ?? []) {
                              if (nearCore.has(nbr)) continue;
                              nearCore.add(nbr);
                              next.add(nbr);
                            }
                          }
                          frontier = next;
                          if (frontier.size === 0) break;
                        }
                      }
                      for (const [refA, candA] of refToCand.entries()) {
                        if (!markedRefAtoms.has(refA)) continue;
                        if (nearCore && !nearCore.has(candA)) continue;
                        directAnchors.add(candA);
                      }
                    }
                    // CATS-based ring pairing: for each MARKED ref ring,
                    // find the cand ring whose pharmacophore-feature
                    // signature is most similar (cosine over family
                    // atom-counts), and add its atoms to directAnchors.
                    // Catches role-equivalence that structural matching
                    // misses — e.g. sorafenib's pyridyl (1 N, aromatic)
                    // is the best CATS match for imatinib's marked
                    // pyrimidine (2 N, aromatic) even though structurally
                    // the pyridyl substruct-matches imatinib's NON-marked
                    // pyridyl decoration.
                    //
                    // GATED on MULTI-RING marked regions (≥2 ring systems).
                    // For single-ring marked regions (e.g. DLK's
                    // pyrimidine-only), CATS would add a SECOND ring
                    // (best-pharmacophore-match) to directAnchors, which
                    // bloats the output beyond the clean ring-swap the
                    // user expects. Multi-ring marked regions already
                    // imply the user wants multiple ring systems shown.
                    const catsPickedAtoms = new Set<number>();
                    if (markedRefRingSigs.length >= 2 && familyQmols && !markedRegionPreserved) {
                      let candAtomLabels: Map<number, Set<string>> | null = null;
                      try {
                        candAtomLabels = labelAtomsByPharmacophore(
                          candMolForErg, familyQmols);
                      } catch {/* CATS labelling failed — skip pass */}
                      if (candAtomLabels) {
                        const candRings = findRingSystems(candAdj);
                        const candRingSigs = candRings.map((r) => {
                          const sig = new Map<string, number>();
                          for (const a of r) {
                            const fams = candAtomLabels!.get(a) ?? new Set<string>();
                            for (const f of fams) sig.set(f, (sig.get(f) ?? 0) + 1);
                          }
                          return sig;
                        });
                        const cosineSig = (
                          a: Map<string, number>, b: Map<string, number>,
                        ): number => {
                          let dot = 0; let na = 0; let nb = 0;
                          for (const [k, v] of a.entries()) {
                            na += v * v;
                            if (b.has(k)) dot += v * b.get(k)!;
                          }
                          for (const v of b.values()) nb += v * v;
                          if (na === 0 || nb === 0) return 0;
                          return dot / Math.sqrt(na * nb);
                        };
                        // For each marked ref ring, find best cand ring.
                        // Use ring-size proximity as a tiebreaker (similar
                        // size = better match for fused-vs-fused vs.
                        // simple-vs-fused mismatch).
                        //
                        // SKIP marked ref rings that are ALREADY covered
                        // by FMCS — track which ref atoms are paired with
                        // anchors and skip rings whose majority of atoms
                        // are already mapped. This prevents CATS from
                        // adding a SECOND cand ring per marked ref ring
                        // for close analogs (e.g. nilotinib already has
                        // its tolyl FMCS-paired with imatinib's marked
                        // tolyl; CATS would otherwise add CF3-phenyl as
                        // "best pharmacophore match" too).
                        const candToRef = new Map<number, number>();
                        if (refToCand) {
                          for (const [refA, candA] of refToCand.entries())
                            candToRef.set(candA, refA);
                        }
                        const usedCandRings = new Set<number>();
                        for (let r = 0; r < markedRefRingSigs.length; r++) {
                          // Coverage check: count this marked ref ring's
                          // atoms that have a cand pair in directAnchors.
                          let covered = 0;
                          for (const refA of markedRefRings[r]) {
                            for (const candA of directAnchors)
                              if (candToRef.get(candA) === refA) {covered++; break;}
                          }
                          if (covered * 2 >= markedRefRings[r].size) continue;
                          const refSig = markedRefRingSigs[r];
                          const refSize = markedRefRings[r].size;
                          let bestIdx = -1;
                          let bestScore = -1;
                          for (let c = 0; c < candRings.length; c++) {
                            if (usedCandRings.has(c)) continue;
                            // Skip cand rings already containing FMCS-
                            // paired anchors (avoid double-counting).
                            let hasAnchor = false;
                            for (const a of candRings[c])
                              if (directAnchors.has(a)) {hasAnchor = true; break;}

                            if (hasAnchor) continue;
                            const sim = cosineSig(refSig, candRingSigs[c]);
                            const sizePenalty = Math.abs(
                              candRings[c].size - refSize) / Math.max(refSize, 1);
                            const score = sim - 0.2 * sizePenalty;
                            if (score > bestScore) {
                              bestScore = score;
                              bestIdx = c;
                            }
                          }
                          // Require minimum similarity to avoid spurious
                          // pairings on chemotype mismatches.
                          if (bestIdx >= 0 && bestScore >= 0.5) {
                            // DISTANCE GATE: the CATS-picked ring must be
                            // CLOSE (<=2 bonds) to an existing anchor.
                            // Otherwise CATS finds the right pharmacophore
                            // signature but at the WRONG position — e.g.
                            // dasatinib's chlorotolyl resembles imatinib's
                            // marked aminotolyl by CATS sig, but it's 4-5
                            // bonds from the pyrimidine anchor (via the
                            // thiazole-amide chain). Same issue with
                            // ponatinib's left pyridyl 3 bonds via alkyne.
                            // K=2 keeps sorafenib's pyridyl (1 bond from
                            // middle phenyl via ether O) while rejecting
                            // these distant matches.
                            if (directAnchors.size > 0) {
                              const reach = new Set<number>(directAnchors);
                              let front = new Set<number>(directAnchors);
                              for (let d = 0; d < 2; d++) {
                                const next = new Set<number>();
                                for (const a of front) {
                                  for (const nbr of candAdj.get(a) ?? []) {
                                    if (reach.has(nbr)) continue;
                                    reach.add(nbr);
                                    next.add(nbr);
                                  }
                                }
                                front = next;
                                if (front.size === 0) break;
                              }
                              let close = false;
                              for (const a of candRings[bestIdx])
                                if (reach.has(a)) {close = true; break;}

                              if (!close) continue;
                            }
                            usedCandRings.add(bestIdx);
                            for (const a of candRings[bestIdx]) {
                              directAnchors.add(a);
                              catsPickedAtoms.add(a);
                            }
                          }
                        }
                        // POST-CATS CLUSTER PRUNING: when CATS picked at
                        // least one cand ring, prune FMCS-paired anchors
                        // that are FAR from the CATS pick. Rationale:
                        // CATS identifies the "correct" pharmacophore
                        // position (e.g. sorafenib's left-half pyridyl
                        // analogous to imatinib's marked pyrimidine).
                        // FMCS pairings that landed on a DIFFERENT region
                        // (e.g. sorafenib's diaryl-urea right half) lie
                        // outside the CATS-endorsed cluster and would
                        // bloat the output. Keep atoms within K=3 bonds
                        // of CATS picks OR in ring systems that touch
                        // that K=3 reach.
                        if (catsPickedAtoms.size > 0) {
                          const catsReach = new Set<number>(catsPickedAtoms);
                          let cfront = new Set<number>(catsPickedAtoms);
                          for (let d = 0; d < 3; d++) {
                            const nxt = new Set<number>();
                            for (const a of cfront) {
                              for (const nbr of candAdj.get(a) ?? []) {
                                if (catsReach.has(nbr)) continue;
                                catsReach.add(nbr);
                                nxt.add(nbr);
                              }
                            }
                            cfront = nxt;
                            if (cfront.size === 0) break;
                          }
                          // Whole-ring extension — if any ring atom is
                          // in catsReach, include the WHOLE ring system.
                          const validAtoms = new Set<number>(catsReach);
                          for (const ring of candRings) {
                            let touchesReach = false;
                            for (const a of ring)
                              if (catsReach.has(a)) {touchesReach = true; break;}

                            if (touchesReach)
                              for (const a of ring) validAtoms.add(a);
                          }
                          // Filter directAnchors — keep only valid atoms.
                          const kept = new Set<number>();
                          for (const a of directAnchors)
                            if (validAtoms.has(a)) kept.add(a);

                          directAnchors.clear();
                          for (const a of kept) directAnchors.add(a);
                        }
                      }
                    }
                    // No-anchor fallback: when NEITHER the full marked-region
                    // SMILES NOR any individual marked ring substruct-matched
                    // in the candidate NOR FMCS-pairing produced anchors
                    // (very distant scaffold hops — e.g. bosutinib's
                    // quinoline doesn't match imatinib's pyrimidine,
                    // ponatinib's imidazo-pyridazine doesn't match either),
                    // seed directAnchors from the LARGEST non-conserved ring
                    // system in the candidate. Rationale: the user's marked
                    // region IS what was replaced; the candidate's novel
                    // ring system (= cand ring atoms not matched by any
                    // non-marked ref ring) plays the same structural role
                    // and is what we want to show.
                    if (directAnchors.size === 0) {
                      const candRingsForFallback = findRingSystems(candAdj);
                      const conservedForFallback = new Set<number>();
                      const refRingsForFallback = findRingSystems(refAdjMap);
                      for (const refRing of refRingsForFallback) {
                        let overlapsMarked = false;
                        for (const a of refRing)
                          if (markedRefAtoms.has(a)) {overlapsMarked = true; break;}

                        if (overlapsMarked) continue;
                        const refRingFrag = extractFragmentSmiles(
                          refMolblock, refRing, rdKitModule);
                        if (!refRingFrag) continue;
                        let refRingQmol: any = null;
                        try {
                          refRingQmol = rdKitModule.get_qmol(refRingFrag);
                          if (!refRingQmol) continue;
                          const j = candMolForErg.get_substruct_match(refRingQmol);
                          if (j && j !== '{}') {
                            const ms: number[] = JSON.parse(j)?.atoms ?? [];
                            for (const a of ms) conservedForFallback.add(a);
                          }
                        } catch {/* ignore */} finally {refRingQmol?.delete?.();}
                      }
                      // Score each cand ring system by non-conserved atoms.
                      let bestRing: Set<number> | null = null;
                      let bestScore = 0;
                      for (const ring of candRingsForFallback) {
                        let novel = 0;
                        for (const a of ring) if (!conservedForFallback.has(a)) novel++;
                        // Require at least half the ring atoms to be novel
                        // and a minimum of 4 atoms — avoids picking a
                        // partly-conserved ring as the "novel" anchor.
                        if (novel >= 4 && novel > bestScore && novel * 2 >= ring.size) {
                          bestScore = novel;
                          bestRing = ring;
                        }
                      }
                      if (bestRing) for (const a of bestRing) directAnchors.add(a);
                    }
                    if (directAnchors.size > 0) {
                      {
                        // Ring-system-based selection:
                        //  - The "marked ring" (directRing): the
                        //    cand ring system containing the
                        //    substruct-matched directAnchors.
                        //  - "Adjacent new rings": ring systems
                        //    that are CONNECTED to the marked
                        //    ring via the cand graph within K
                        //    bonds (we use K=3 to allow for
                        //    pyrimidine-NH-thiazole-style 2-hop
                        //    cases) AND aren't substruct-matched
                        //    by any non-marked ref ring system
                        //    (= aren't conserved).
                        // The linker ATOMS between rings stay
                        // out of the result so the cell shows
                        // just the ring SMILES, no dangling NH /
                        // methylene atoms.
                        const rings = findRingSystems(candAdj);
                        // Build per-atom ring-system index.
                        const atomToRing = new Map<number, number>();
                        for (let r = 0; r < rings.length; r++)
                          for (const a of rings[r]) atomToRing.set(a, r);

                        // Identify "conserved" cand rings via
                        // substruct match against each non-marked
                        // ref ring system. A ref ring system is
                        // the connected component of ref's
                        // non-marked atoms restricted to ring
                        // atoms; we substruct-match each one's
                        // SMILES against cand. Matched cand atoms
                        // are added to a `conservedByMatch` set.
                        const conservedByMatch = new Set<number>();
                        const refRings = findRingSystems(refAdjMap);
                        for (const refRing of refRings) {
                          // Skip ref rings that overlap with the
                          // marked region (those are the marked
                          // region itself, not the conserved core).
                          let overlapsMarked = false;
                          for (const a of refRing)
                            if (markedRefAtoms.has(a)) {overlapsMarked = true; break;}

                          if (overlapsMarked) continue;
                          const refRingFrag = extractFragmentSmiles(
                            refMolblock, refRing, rdKitModule);
                          if (!refRingFrag) continue;
                          let refRingQmol: any = null;
                          try {
                            refRingQmol = rdKitModule.get_qmol(refRingFrag);
                            if (!refRingQmol) continue;
                            const j = candMolForErg.get_substruct_match(refRingQmol);
                            if (j && j !== '{}') {
                              const ms: number[] = JSON.parse(j)?.atoms ?? [];
                              for (const a of ms) conservedByMatch.add(a);
                            }
                          } catch {/* ignore */} finally {refRingQmol?.delete?.();}
                        }
                        // Atom-level BFS — walk through the cand
                        // graph from directAnchors, skipping
                        // conserved atoms. This naturally walks
                        // through LINKER atoms (NH, CH2 etc.)
                        // that aren't in any ring but aren't
                        // conserved either, reaching ring systems
                        // 2-3 bonds away (e.g. dasatinib's
                        // thiazole, separated from pyrimidine by
                        // one NH linker).
                        // Count ring systems the MARKED REGION
                        // spans — drives linker inclusion (only
                        // 2+ ring marked regions add linker atoms,
                        // so DLK-style single-ring swaps stay
                        // clean).
                        const markedRingSystemCount = (() => {
                          const refRingsForCount = findRingSystems(refAdjMap);
                          let n = 0;
                          for (const refRing of refRingsForCount) {
                            for (const a of refRing)
                              if (markedRefAtoms.has(a)) {n++; break;}
                          }
                          return n;
                        })();
                        const reachable = new Set<number>(directAnchors);
                        let frontier = new Set<number>(directAnchors);
                        for (let depth = 0; depth < 3; depth++) {
                          const next = new Set<number>();
                          for (const a of frontier) {
                            for (const nbr of candAdj.get(a) ?? []) {
                              if (reachable.has(nbr)) continue;
                              if (conservedByMatch.has(nbr)) continue;
                              reachable.add(nbr);
                              next.add(nbr);
                            }
                          }
                          frontier = next;
                          if (frontier.size === 0) break;
                        }
                        // Output: union of ring atoms from any ring
                        // system that has at least one reachable
                        // atom PLUS linker atoms ON THE PATH
                        // between those rings. Pure linker atoms
                        // (NH, CH2 etc. connecting the rings) are
                        // INCLUDED so the cell shows a single
                        // connected fragment instead of `.`-
                        // separated pieces — e.g. dasatinib's
                        // pyrimidine-NH-thiazole renders as one
                        // SMILES rather than `pyrimidine . thiazole`.
                        //
                        // Rule: include an atom IF either
                        //   (a) it's in any selected ring system, OR
                        //   (b) it's a non-ring atom that's
                        //       reachable AND not conservedByMatch.
                        // The non-ring constraint excludes atoms
                        // that happen to be on the conserved core's
                        // boundary; they pass the BFS skip-check
                        // but shouldn't drag the conserved core
                        // into the result.
                        // Build cand atom-element map AND in-ring
                        // double-bond pairs for the adjacent-ring
                        // significance check. RDKit's molblock
                        // output is KEKULÉ (no type-4 aromatic
                        // bonds) — so aromaticity is detected
                        // structurally: a ring is aromatic if it
                        // has ~ring.size/2 double bonds (3 for
                        // benzene, 2 for pyrrole/thiophene).
                        const candAtomEl: Map<number, string> = new Map();
                        const candDoubleBonds = new Set<string>();
                        const candBondKey = (a: number, b: number) =>
                          a < b ? `${a}-${b}` : `${b}-${a}`;
                        if (candMolblock) {
                          const lines = candMolblock.split('\n');
                          const counts = lines[3];
                          const nA = parseInt(counts?.substring(0, 3) ?? '0');
                          const nB = parseInt(counts?.substring(3, 6) ?? '0');
                          if (Number.isFinite(nA)) {
                            for (let i = 0; i < nA; i++) {
                              const ln = lines[4 + i];
                              if (ln && ln.length >= 34)
                                candAtomEl.set(i, ln.substring(31, 34).trim());
                            }
                          }
                          if (Number.isFinite(nA) && Number.isFinite(nB)) {
                            for (let j = 0; j < nB; j++) {
                              const ln = lines[4 + nA + j];
                              if (!ln || ln.length < 9) continue;
                              const bt = parseInt(ln.substring(6, 9));
                              const a = parseInt(ln.substring(0, 3)) - 1;
                              const b = parseInt(ln.substring(3, 6)) - 1;
                              if (bt === 2 || bt === 4)
                                candDoubleBonds.add(candBondKey(a, b));
                            }
                          }
                        }
                        const isRingAromatic = (ring: Set<number>): boolean => {
                          let dbl = 0;
                          const arr = [...ring];
                          for (let i = 0; i < arr.length; i++) {
                            for (let j = i + 1; j < arr.length; j++) {
                              if (candDoubleBonds.has(candBondKey(arr[i], arr[j])))
                                dbl++;
                            }
                          }
                          // Benzene: 3 dbl / 6 atoms. Pyrrole: 2/5.
                          // Heuristic: at least floor(size/3) doubles
                          // catches all standard aromatics.
                          return dbl >= Math.floor(ring.size / 3);
                        };
                        // `markedRingSystemCount` was computed
                        // earlier (above the BFS) and drives:
                        //   (a) BFS depth (5 vs 3) so multi-ring
                        //       marked regions reach further.
                        //   (b) skipping the "novel >= 4" check on
                        //       adjacent rings for 2+ ring marked
                        //       regions (so sorafenib's pyridyl
                        //       and ponatinib's left pyridyl still
                        //       count as role-equivalent rings).
                        //   (c) linker-atom inclusion (below).
                        const selectedRingAtoms = new Set<number>();
                        for (const ring of rings) {
                          // ANCHOR ring: always include (contains
                          // directAnchors — the user explicitly
                          // marked these atoms via the full or per-
                          // ring substruct match, FMCS pairing, or
                          // the no-anchor fallback). Add ALL ring
                          // atoms, ignoring conservedByMatch — the
                          // pure-benzene ref ring substruct happens
                          // to match the marked region's tolyl too,
                          // but we shouldn't drop user-marked atoms.
                          let containsAnchor = false;
                          for (const a of ring)
                            if (directAnchors.has(a)) {containsAnchor = true; break;}

                          if (containsAnchor) {
                            for (const a of ring) selectedRingAtoms.add(a);
                            continue;
                          }
                          // MARKED-REGION-PRESERVED gate: when the
                          // full marked-region SMILES substruct-
                          // matched the candidate (close-analog
                          // cases like lapatinib / afatinib which
                          // share gefitinib's anilinoquinazoline),
                          // DON'T include adjacent novel rings at
                          // substituent positions. The candidate
                          // already has the user's marked region
                          // verbatim; substituents like lapatinib's
                          // furan-at-position-6 or benzyl-fluoro-
                          // benzyl arm are decoration, not the
                          // pharmacophore the user marked.
                          if (markedRegionPreserved) continue;
                          // NON-ANCHOR rings: include the WHOLE ring
                          // if it (a) has a reachable atom OR an
                          // atom adjacent to a reachable atom (so
                          // fused rings whose conserved half is
                          // BFS-blocked still get detected via the
                          // novel half's neighbour), (b) has >= 4
                          // novel atoms (atoms not in
                          // conservedByMatch), and (c) is
                          // "significant" — >= 5 atoms with a
                          // heteroatom OR aromatic. The whole-ring
                          // inclusion (vs. just novel-atom subset)
                          // keeps the SMILES output a valid ring;
                          // a partial ring is chemical nonsense.
                          let isReachableOrAdjacent = false;
                          for (const a of ring) {
                            if (reachable.has(a)) {isReachableOrAdjacent = true; break;}
                            for (const nbr of candAdj.get(a) ?? [])
                              if (reachable.has(nbr)) {isReachableOrAdjacent = true; break;}

                            if (isReachableOrAdjacent) break;
                          }
                          if (!isReachableOrAdjacent) continue;
                          if (ring.size < 5) continue;
                          let hasHetero = false;
                          for (const a of ring) {
                            const el = candAtomEl.get(a) ?? 'C';
                            if (el === 'N' || el === 'O' || el === 'S') hasHetero = true;
                          }
                          const isAromatic = isRingAromatic(ring);
                          if (!hasHetero && !isAromatic) continue;
                          // Require >= 4 novel atoms — non-anchor
                          // rings have to bring something new to
                          // earn inclusion. Close-analog adjacent
                          // rings (e.g. nilotinib's pyridyl
                          // substituent, identical to imatinib's
                          // non-marked pyridyl ref ring) score
                          // novel=0 and get correctly excluded.
                          // Bosutinib's quinoline has 4 benzo-half
                          // novel atoms so it passes; sorafenib's
                          // chloro-CF3-phenyl is wholly conserved
                          // (matched by ref benzene) so it doesn't.
                          let novelCount = 0;
                          for (const a of ring) if (!conservedByMatch.has(a)) novelCount++;
                          if (novelCount < 4) continue;
                          // Whole-ring inclusion (no conservation
                          // filter) so the SMILES stays a closed
                          // ring even when half the atoms happen to
                          // be in conservedByMatch (e.g. bosutinib's
                          // quinoline pyridyl half).
                          for (const a of ring) selectedRingAtoms.add(a);
                        }
                        const out = new Set<number>(selectedRingAtoms);
                        // Build per-atom ring membership for the
                        // linker-inclusion check below.
                        const atomInAnyRing = new Set<number>();
                        for (const ring of rings)
                          for (const a of ring) atomInAnyRing.add(a);
                        // Always include FMCS-paired NON-RING atoms
                        // in directAnchors. The marked region may
                        // include linker atoms (e.g. the amide
                        // N-C(=O) atoms in `c1ccc(C(=O)N)cc1`); when
                        // FMCS pairs them with cand atoms, those
                        // pairings carry positional information the
                        // user explicitly marked. Without this, the
                        // output drops the amide and shows just the
                        // ring (benzamide_phenyl → dasatinib loses
                        // the amide; aminotolyl → dasatinib loses
                        // the NH + methyl). Restricted to non-ring
                        // atoms — ring atoms are handled by the
                        // anchor-containing-ring rule above.
                        for (const a of directAnchors)
                          if (!atomInAnyRing.has(a)) out.add(a);

                        // Decide whether to add LINKER ATOMS based on
                        // the structure of the MARKED REGION in ref.
                        // If the user marked a multi-ring region (e.g.
                        // pyrimidine + NH + aminotolyl, 2 ring
                        // systems), we want the candidate's output to
                        // mirror that shape — multiple rings connected
                        // by linker atoms (dasatinib's thiazole–NH–
                        // pyrimidine). If the marked region is a
                        // single ring (e.g. DLK's pyrimidine alone),
                        // the candidate should show ONLY the
                        // replacement ring (DLK's pyrazole), without
                        // pendant N-isopropyl / cyclopropyl side
                        // chains. `markedRingSystemCount` was
                        // computed above.
                        if (markedRingSystemCount >= 2) {
                          // Add linker atoms only on SHORTEST
                          // PATHS between selected ring systems —
                          // not all reachable non-ring atoms.
                          // Otherwise substituents of anchor rings
                          // (e.g. nilotinib's amide C=O, which is
                          // a substituent of the central tolyl but
                          // NOT between two selected rings) get
                          // swept in. The shortest-path
                          // restriction keeps the cell tight:
                          // dasatinib's NH between pyrimidine and
                          // thiazole is included; nilotinib's
                          // amide outward from aminotolyl isn't.
                          const ringSystemIds: number[][] = [];
                          for (const ring of rings) {
                            let hasSelected = false;
                            for (const a of ring)
                              if (selectedRingAtoms.has(a)) {hasSelected = true; break;}

                            if (hasSelected) ringSystemIds.push([...ring]);
                          }
                          // BFS from each selected ring to find
                          // shortest paths to every other selected
                          // ring; include atoms on those paths.
                          const addPathAtoms = (from: number[], target: Set<number>) => {
                            const parent = new Map<number, number>();
                            const visited = new Set<number>(from);
                            let frontier = new Set<number>(from);
                            let found: number | null = null;
                            while (frontier.size > 0 && found === null) {
                              const next = new Set<number>();
                              for (const a of frontier) {
                                if (target.has(a) && !from.includes(a)) {
                                  found = a;
                                  break;
                                }
                                for (const nbr of candAdj.get(a) ?? []) {
                                  if (visited.has(nbr)) continue;
                                  visited.add(nbr);
                                  parent.set(nbr, a);
                                  next.add(nbr);
                                }
                              }
                              if (found !== null) break;
                              frontier = next;
                            }
                            if (found !== null) {
                              let cur: number | undefined = found;
                              while (cur !== undefined && !from.includes(cur)) {
                                out.add(cur);
                                cur = parent.get(cur);
                              }
                            }
                          };
                          for (let i = 0; i < ringSystemIds.length; i++) {
                            const fromAtoms = ringSystemIds[i];
                            for (let j = i + 1; j < ringSystemIds.length; j++) {
                              const targetSet = new Set<number>(ringSystemIds[j]);
                              addPathAtoms(fromAtoms, targetSet);
                            }
                          }
                        }
                        // MARKED-REGION-PRESERVED override: when
                        // the candidate contains the full marked
                        // region verbatim (e.g. afatinib's
                        // anilinoquinazoline = gefitinib's marked
                        // region), the MCS-direct path may have
                        // ALREADY put substituent atoms into
                        // `replacedAtoms` (afatinib's oxetane was
                        // tagged "replaced" because it doesn't
                        // pair with gefitinib's morpholine).
                        // OVERRIDE replacedAtoms with the
                        // substruct-fallback output (= just the
                        // marked-region analog), dropping any
                        // substituent-at-the-wrong-position
                        // contributions from MCS-direct.
                        if (markedRegionPreserved)
                          replacedAtoms = new Set<number>(out);
                        else
                          for (const a of out) replacedAtoms.add(a);

                        // POST-CATS GLOBAL PRUNING: if CATS picked
                        // any rings, also filter the MCS-direct
                        // contributions to replacedAtoms (which
                        // were added BEFORE the substruct fallback
                        // ran). The MCS-direct path may have placed
                        // atoms on the OPPOSITE side of the
                        // molecule from where CATS identified the
                        // marked-region analog (e.g. sorafenib:
                        // MCS-direct pairs imatinib's marked atoms
                        // with sorafenib's urea+right-phenyl half,
                        // but CATS endorses the pyridyl+middle-
                        // phenyl half as the actual pharmacophore
                        // analog). Filter to the CATS-endorsed
                        // cluster.
                        if (catsPickedAtoms.size > 0) {
                          const catsReach2 = new Set<number>(catsPickedAtoms);
                          let cfront2 = new Set<number>(catsPickedAtoms);
                          for (let d = 0; d < 3; d++) {
                            const nxt = new Set<number>();
                            for (const a of cfront2) {
                              for (const nbr of candAdj.get(a) ?? []) {
                                if (catsReach2.has(nbr)) continue;
                                catsReach2.add(nbr);
                                nxt.add(nbr);
                              }
                            }
                            cfront2 = nxt;
                            if (cfront2.size === 0) break;
                          }
                          const candRings2 = findRingSystems(candAdj);
                          const validAtoms2 = new Set<number>(catsReach2);
                          for (const ring of candRings2) {
                            let touches = false;
                            for (const a of ring)
                              if (catsReach2.has(a)) {touches = true; break;}

                            if (touches)
                              for (const a of ring) validAtoms2.add(a);
                          }
                          const keptRA = new Set<number>();
                          for (const a of replacedAtoms)
                            if (validAtoms2.has(a)) keptRA.add(a);

                          replacedAtoms = keptRA;
                        }
                      }
                    }
                  } catch {/* substruct fallback failed — leave atoms set from ring-sub */}
                }
                // DISCONNECTED-OUTPUT GUARD: when `replacedAtoms` spans
                // multiple connected components in the candidate graph
                // (FMCS pairs the marked region's atoms onto two disjoint
                // clusters in cand — e.g. aminotolyl → dasatinib's
                // thiazole AND chlorotolyl), keep only the SMALLEST
                // component. The smallest fragment is usually the
                // chemotype-swap analog at the central / hinge position
                // (dasatinib: thiazole 5 atoms ≪ chlorotolyl-amide 7+;
                // bosutinib: chloroaniline 6 ≪ quinoline 10) — that
                // matches the user's "show the chemotype-swap analog,
                // not the structurally-identical decoration" intuition.
                if (replacedAtoms.size > 1 && candAdj) {
                  // Compute connected components within replacedAtoms.
                  const comps: Set<number>[] = [];
                  const seenComp = new Set<number>();
                  for (const a of replacedAtoms) {
                    if (seenComp.has(a)) continue;
                    const comp = new Set<number>([a]);
                    let frontier = [a];
                    seenComp.add(a);
                    while (frontier.length > 0) {
                      const next: number[] = [];
                      for (const x of frontier) {
                        for (const nbr of candAdj.get(x) ?? []) {
                          if (!replacedAtoms.has(nbr) || seenComp.has(nbr)) continue;
                          seenComp.add(nbr);
                          comp.add(nbr);
                          next.push(nbr);
                        }
                      }
                      frontier = next;
                    }
                    comps.push(comp);
                  }
                  if (comps.length > 1) {
                    comps.sort((a, b) => a.size - b.size);
                    replacedAtoms = comps[0];
                  }
                }
                // SINGLE-RING-MARKED RESTRICTION: if the user marked
                // exactly ONE ring system in ref (e.g. morpholine alone,
                // chlorofluoroaniline alone, pyrrolopyrimidine alone),
                // the output should also be a single ring system. The
                // multi-ring outputs that sometimes leak through
                // (lapatinib chloroaniline+phenyl, osimertinib pyrimidine
                // +indolyl, bosutinib chloroaniline+quinoline) come from
                // FMCS pairings spanning multiple rings. Restrict to the
                // ring system containing the MOST replacedAtoms + its
                // immediately attached non-ring anchors.
                if (replacedAtoms.size > 0 && candAdj && markedRefAtoms.size > 0) {
                  const refRingsCount = findRingSystems(refAdjMap);
                  let markedRingSystems = 0;
                  for (const refRing of refRingsCount) {
                    for (const a of refRing)
                      if (markedRefAtoms.has(a)) {markedRingSystems++; break;}
                  }
                  if (markedRingSystems === 1) {
                    const candRings = findRingSystems(candAdj);
                    // Parse cand atom elements for halogen check.
                    const candEl: Map<number, string> = new Map();
                    if (candMolblock) {
                      const lines2 = candMolblock.split('\n');
                      const nA2 = parseInt(lines2[3]?.substring(0, 3) ?? '0');
                      if (Number.isFinite(nA2)) {
                        for (let i = 0; i < nA2; i++) {
                          const ln = lines2[4 + i];
                          if (ln && ln.length >= 34)
                            candEl.set(i, ln.substring(31, 34).trim());
                        }
                      }
                    }
                    // Parse ref atom elements.
                    const refEl: Map<number, string> = new Map();
                    if (refMolblockForLocal) {
                      const lines2 = refMolblockForLocal.split('\n');
                      const nA2 = parseInt(lines2[3]?.substring(0, 3) ?? '0');
                      if (Number.isFinite(nA2)) {
                        for (let i = 0; i < nA2; i++) {
                          const ln = lines2[4 + i];
                          if (ln && ln.length >= 34)
                            refEl.set(i, ln.substring(31, 34).trim());
                        }
                      }
                    }
                    // Count "substantial" substituents on a ring system:
                    // each connected branch off the ring with >= 3 non-
                    // halogen heavy atoms. Captures chains-to-rings,
                    // amide arms, NH-aryl substituents; rejects F, Cl,
                    // Br, I, methyl, methoxy, CF3 (the "small
                    // decorations" the user wants to ignore).
                    const countSubstantial = (
                      ringSys: Set<number>,
                      adj: Map<number, number[]>,
                      el: Map<number, string>,
                    ): number => {
                      const visited = new Set<number>();
                      const subs: Set<number>[] = [];
                      for (const ringAtom of ringSys) {
                        for (const nbr of adj.get(ringAtom) ?? []) {
                          if (ringSys.has(nbr) || visited.has(nbr)) continue;
                          const comp = new Set<number>([nbr]);
                          visited.add(nbr);
                          let frontier = [nbr];
                          while (frontier.length > 0) {
                            const next: number[] = [];
                            for (const a of frontier) {
                              for (const n of adj.get(a) ?? []) {
                                if (ringSys.has(n) || visited.has(n)) continue;
                                visited.add(n);
                                comp.add(n);
                                next.push(n);
                              }
                            }
                            frontier = next;
                          }
                          subs.push(comp);
                        }
                      }
                      let substantial = 0;
                      for (const sub of subs) {
                        let nonHal = 0;
                        for (const a of sub) {
                          const e = el.get(a) ?? 'C';
                          if (e !== 'F' && e !== 'Cl' && e !== 'Br' && e !== 'I') nonHal++;
                        }
                        if (nonHal >= 3) substantial++;
                      }
                      return substantial;
                    };
                    // Compute marked ref ring's substantial count.
                    let markedRefRing: Set<number> | null = null;
                    for (const refRing of refRingsCount) {
                      for (const a of refRing)
                        if (markedRefAtoms.has(a)) {markedRefRing = refRing; break;}

                      if (markedRefRing) break;
                    }
                    const refSubstantial = markedRefRing ?
                      countSubstantial(markedRefRing, refAdjMap, refEl) : 0;
                    // Compute graph distance from each ring system to the
                    // LARGEST ring system in the molecule (used as a
                    // proxy for "molecular centroid"). Marked ref ring
                    // and best cand ring should be at similar depths —
                    // morpholine in gefitinib is ~5 bonds from the
                    // quinazoline core, so the analog should also be far
                    // from the cand's core ring (fluorobenzyl in
                    // lapatinib, not chloroaniline at depth 1).
                    const depthToLargestRing = (
                      ringSys: Set<number>, allRings: Set<number>[],
                      adj: Map<number, number[]>,
                    ): number => {
                      // Largest ring system = the one with the most atoms.
                      let largest = allRings[0]; let maxSz = 0;
                      for (const r of allRings) {
                        if (r === ringSys) continue;
                        if (r.size > maxSz) {maxSz = r.size; largest = r;}
                      }
                      if (!largest || largest === ringSys) return 0;
                      // BFS from ringSys atoms to find shortest path to
                      // largest ring system.
                      const visited = new Set<number>(ringSys);
                      let frontier = [...ringSys];
                      let depth = 0;
                      while (frontier.length > 0) {
                        const next: number[] = [];
                        for (const a of frontier) {
                          if (largest.has(a) && !ringSys.has(a)) return depth;
                          for (const n of adj.get(a) ?? []) {
                            if (visited.has(n)) continue;
                            visited.add(n);
                            next.push(n);
                          }
                        }
                        frontier = next;
                        depth++;
                        if (depth > 20) return depth;
                      }
                      return depth;
                    };
                    const refDepth = markedRefRing ?
                      depthToLargestRing(markedRefRing, refRingsCount, refAdjMap) : 0;
                    // Pick cand ring system whose (substantialCount,
                    // depth) best matches (refSubstantial, refDepth).
                    // Considers ALL cand rings (not just FMCS-paired
                    // ones) — when the user marks an OUTER decoration
                    // ring (morpholine, chlorofluoroaniline), the
                    // position-analog cand ring may not be FMCS-paired
                    // at all (FMCS focuses on the central scaffold).
                    // Two-pass search: track BEST FMCS-paired ring and
                    // BEST non-FMCS ring separately.
                    let bestFmcsIdx = -1;
                    let bestFmcsScore: [number, number, number] = [Infinity, Infinity, 0];
                    let bestNonFmcsIdx = -1;
                    let bestNonFmcsScore: [number, number, number] = [Infinity, Infinity, 0];
                    for (let r = 0; r < candRings.length; r++) {
                      let inRA = 0;
                      for (const a of candRings[r]) if (replacedAtoms.has(a)) inRA++;
                      const candSub = countSubstantial(candRings[r], candAdj, candEl);
                      const subDiff = Math.abs(candSub - refSubstantial);
                      const candDepth = depthToLargestRing(candRings[r], candRings, candAdj);
                      const depthDiff = Math.abs(candDepth - refDepth);
                      const score: [number, number, number] = [subDiff, depthDiff, -inRA];
                      const beats = (s: [number, number, number], best: [number, number, number]) =>
                        s[0] < best[0] ||
                        (s[0] === best[0] && s[1] < best[1]) ||
                        (s[0] === best[0] && s[1] === best[1] && s[2] < best[2]);
                      if (inRA > 0) {
                        if (bestFmcsIdx < 0 || beats(score, bestFmcsScore)) {
                          bestFmcsScore = score;
                          bestFmcsIdx = r;
                        }
                      } else {
                        if (bestNonFmcsIdx < 0 || beats(score, bestNonFmcsScore)) {
                          bestNonFmcsScore = score;
                          bestNonFmcsIdx = r;
                        }
                      }
                    }
                    // Narrow FMCS-preference relaxation: allow a non-
                    // FMCS-paired ring to win IF (a) it has exactly
                    // diff=0 substantial-count match, (b) the best
                    // FMCS-paired ring has diff ≥ 1, (c) the non-FMCS
                    // ring has ≥ 5 atoms (rules out tiny aliphatic
                    // carbocycle noise). This fixes lapatinib morpholine
                    // → fluorobenzyl (diff=0, non-FMCS) over
                    // chloroaniline (diff=1, FMCS) without affecting
                    // BCR-ABL / DLK tests where the FMCS-paired ring
                    // typically also has diff=0.
                    let bestRingIdx = bestFmcsIdx;
                    if (bestFmcsIdx >= 0 && bestNonFmcsIdx >= 0 &&
                        bestNonFmcsScore[0] === 0 && bestFmcsScore[0] >= 1 &&
                        candRings[bestNonFmcsIdx].size >= 5)
                      bestRingIdx = bestNonFmcsIdx;

                    // Safety: only apply restriction if the picked ring
                    // makes sense; otherwise keep replacedAtoms as-is.
                    if (bestRingIdx < 0) bestRingIdx = -1;
                    if (bestRingIdx >= 0) {
                      const keep = new Set<number>(candRings[bestRingIdx]);
                      // Add non-ring replacedAtoms ADJACENT to kept ring.
                      const atomInAnyRing = new Set<number>();
                      for (const r of candRings) for (const a of r) atomInAnyRing.add(a);
                      for (const a of replacedAtoms) {
                        if (atomInAnyRing.has(a)) continue;
                        for (const nbr of candAdj.get(a) ?? [])
                          if (keep.has(nbr)) {keep.add(a); break;}
                      }
                      // OVERRIDE replacedAtoms with the chosen ring +
                      // adjacent substituents. The substantial-count
                      // matching is the authoritative signal — when it
                      // picks a non-FMCS-paired ring, we use that ring's
                      // atoms as the analog (e.g. lapatinib morpholine
                      // → fluorobenzyl, pyridyl_substituent → dasatinib
                      // chlorotolyl: both correctly chosen as 1-
                      // substantial outer-phenyl analogs of 1-
                      // substantial outer-ring marked regions).
                      replacedAtoms = keep;
                    }
                  }
                }
                if (replacedAtoms.size > 0 && candMolblock) {
                  // Compute bonds inside the selected ring system.
                  const lines = candMolblock.split('\n');
                  const counts = lines[3];
                  const nA = parseInt(counts.substring(0, 3));
                  const nB = parseInt(counts.substring(3, 6));
                  if (Number.isFinite(nA) && Number.isFinite(nB)) {
                    for (let j = 0; j < nB; j++) {
                      const bondLine = lines[4 + nA + j];
                      if (!bondLine || bondLine.length < 6) continue;
                      const a = parseInt(bondLine.substring(0, 3)) - 1;
                      const b = parseInt(bondLine.substring(3, 6)) - 1;
                      if (replacedAtoms.has(a) && replacedAtoms.has(b))
                        replacedBonds.push(j);
                    }
                  }
                }
              }

              // Fallback: ErG pharmacophore matching when ring-substitution
              // found nothing (acyclic marked region, no rings near the
              // boundary, or no refToCand available).
              if (replacedAtoms.size === 0 && familyQmols) {
                const ergResult = computeErgSharedAtoms(
                  refMolForErg, candMolForErg, markedRefAtoms, familyQmols,
                  rdKitModule);
                if (ergResult.candAtoms.size > 0) {
                  replacedAtoms = ergResult.candAtoms;
                  replacedBonds = ergResult.candBonds;
                }
              }

              if (replacedAtoms.size > 0) {
                // Compute bonds inside the chosen atom set (may have changed
                // since the first pass if the substruct fallback extended
                // the atom set).
                replacedBonds = [];
                if (candMolblock) {
                  const lines = candMolblock.split('\n');
                  const counts = lines[3];
                  const nA = parseInt(counts.substring(0, 3));
                  const nB = parseInt(counts.substring(3, 6));
                  if (Number.isFinite(nA) && Number.isFinite(nB)) {
                    for (let j = 0; j < nB; j++) {
                      const bondLine = lines[4 + nA + j];
                      if (!bondLine || bondLine.length < 6) continue;
                      const a = parseInt(bondLine.substring(0, 3)) - 1;
                      const b = parseInt(bondLine.substring(3, 6)) - 1;
                      if (replacedAtoms.has(a) && replacedAtoms.has(b))
                        replacedBonds.push(j);
                    }
                  }
                }
                const highlightAtomColors: {[k: number]: number[]} = {};
                for (const a of replacedAtoms)
                  highlightAtomColors[a] = [...sharedColorRgb];
                const highlightBondColors: {[k: number]: number[]} = {};
                for (const b of replacedBonds)
                  highlightBondColors[b] = [...sharedColorRgb];
                sharedSubstructs.set(rowIdx, {
                  atoms: [...replacedAtoms],
                  bonds: replacedBonds,
                  highlightAtomColors,
                  highlightBondColors,
                });
                if (candMolblock) {
                  // Use the with-RGroups molblock extractor so the cell
                  // renders the ring + R1/R2/R3 attachment markers at
                  // each boundary bond — matching the Local-mode output
                  // and giving the chemist a visual cue of where the
                  // substituents attach. Plain SMILES extraction loses
                  // this information (an N-substituted pyrazole renders
                  // as "c1cn[nH]c1" instead of "c1cn(R1)nc1R2-R3-style"
                  // — the user sees just the ring atoms with implicit
                  // hydrogens, not the substitution pattern). Falls
                  // back to plain SMILES on extraction failure (e.g.
                  // valence issues at the boundary).
                  let frag = extractFragmentMolblockWithRGroups(
                    candMolblock, replacedAtoms, rdKitModule);
                  if (!frag) {
                    frag = extractFragmentSmiles(
                      candMolblock, replacedAtoms, rdKitModule);
                  }
                  if (frag) sharedFragments.set(rowIdx, frag);
                }
              }
            }
          }
        } catch {/* extraction failure on this row — leave both views empty */} finally {
          candMolForErg?.delete();
          refMolForErg?.delete();
        }
      }
    } catch (e: any) {
      // Outer FMCS row-handler safety net. Before this catch existed the
      // outer try at the top of this loop body had only a `finally`,
      // meaning a single bad row (malformed SMARTS from get_qmol, an
      // RDKit panic deeper in the chain, a JSON.parse failure on a
      // substruct-match output) propagated past `finally` and aborted
      // the entire sweep — `writeOutputColumns` was never reached and
      // the user got NO Replaced Region / Replacement columns at all.
      //
      // The cooperative cancel check is at line 350, BEFORE this try
      // block opens, so any throw caught here is genuinely a per-row
      // failure (not user cancellation). Log + skip is the right
      // response: the sweep keeps going, the failed row stays at
      // FNULL ratio with empty fragment cells.
      console.warn(`[scaffold-hopping] FMCS row ${rowIdx} threw: ` +
        `${e?.message ?? e}; skipping this row, sweep continues.`);
    } finally {
      mcsMol?.delete();
    }

    if ((k % 20) === 0) {
      const pct = 60 + Math.floor(35 * (k / Math.max(1, topRows.length)));
      progress.update(pct, `MCS ${k + 1}/${topRows.length}...`);
    }
  }
  if (familyQmols) disposeFamilyQmols(familyQmols);
  return {
    refMcsAtoms: out, sharedSubstructs, sharedFragments,
    refToCandByRow, candAdjByRow, refAdj: refAdjMap,
  };
}
