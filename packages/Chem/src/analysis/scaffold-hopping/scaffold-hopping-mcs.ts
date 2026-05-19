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
} from './erg-scaffold-match';
import {
  extractFragmentSmiles, parseV2000Adjacency, extractFragmentMolblockWithRGroups,
  buildR1AnchorMolblock,
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
}> {
  const out = new Map<number, Set<number>>();
  const sharedSubstructs = new Map<number, ISubstruct>();
  const sharedFragments = new Map<number, string>();
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
  // adjacency for "edge anchor" detection. Fetched once outside the
  // per-row loop since the reference doesn't change between rows.
  let refMolblockForLocal = '';
  if (useRGroup && markedRefAtoms.size > 0) {
    let tmpRefMol: any = null;
    try {
      tmpRefMol = rdKitModule.get_mol(refSmiles);
      refMolblockForLocal = tmpRefMol?.get_molblock?.() ?? '';
    } catch {/* leave empty — Local-mode falls back to largest-component */}
    finally {tmpRefMol?.delete();}
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
  if (refAtomCount === 0) return {refMcsAtoms: out, sharedSubstructs, sharedFragments};

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
  } else try {
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
      try {
        refMol = rdKitModule.get_mol(refSmiles);
        const refMatchJson = refMol.get_substruct_match(mcsMol);
        if (refMatchJson && refMatchJson !== '{}') {
          const parsed = JSON.parse(refMatchJson);
          const refMatchAtoms: number[] = parsed?.atoms ?? [];
          out.set(rowIdx, new Set(refMatchAtoms));
        }
      } catch {/* ignore */} finally {
        refMol?.delete();
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
          const candForLocal = molecules.get(rowIdx);
          if (candForLocal && sharedColorRgb) {
            candMolForLocal = rdKitModule.get_mol(candForLocal);
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
              const refAdj = parseV2000Adjacency(refMolblockForLocal);
              const refBorderAtoms = new Set<number>();
              for (const m of markedRefAtoms) {
                for (const nbr of refAdj.get(m) ?? [])
                  if (!markedRefAtoms.has(nbr)) refBorderAtoms.add(nbr);
              }
              const edgeAnchors = new Set<number>();
              for (let i = 0; i < pairLen; i++) {
                if (refBorderAtoms.has(refMcsAtomsArr[i]))
                  edgeAnchors.add(candMcsAtomsArr[i]);
              }

              // BFS components on non-preserved candidate atoms.
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

              // Pick components that either (a) contain a direct anchor
              // or (b) have an atom adjacent (in candidate graph) to an
              // edge anchor. Union them — multi-region marked-area
              // replacements still come through whole.
              const includedComponents: Set<number>[] = [];
              for (const c of components) {
                let touches = false;
                for (const a of c) {
                  if (directAnchors.has(a)) {touches = true; break;}
                  let adjacentEdge = false;
                  for (const nbr of candAdj.get(a) ?? []) {
                    if (edgeAnchors.has(nbr)) {adjacentEdge = true; break;}
                  }
                  if (adjacentEdge) {touches = true; break;}
                }
                if (touches) includedComponents.push(c);
              }

              // If no component touched an anchor at all (pathological:
              // marked region fully outside any MCS pairing AND not
              // adjacent to any preserved boundary atom), fall back to
              // the largest non-preserved component so the column at
              // least shows SOMETHING rather than empty.
              const chosen = new Set<number>();
              if (includedComponents.length > 0) {
                for (const c of includedComponents) for (const a of c) chosen.add(a);
              } else if (components.length > 0) {
                let largest = components[0];
                for (const c of components) if (c.size > largest.size) largest = c;
                for (const a of largest) chosen.add(a);
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
        } catch {/* Local-mode extraction failure — leave the row empty */}
        finally {candMolForLocal?.delete(); refMolForLocal?.delete();}
      } else {
        // ErG path. Pharmacophore-equivalent matching with path-completion.
        // Used for Easy / Middle / Hard presets where the "kept" region
        // may shift chemotype (pyrimidine ↔ thiazole etc.).
        let candMolForErg: any = null;
        let refMolForErg: any = null;
        try {
          const candForErg = molecules.get(rowIdx);
          if (candForErg && familyQmols && sharedColorRgb) {
            candMolForErg = rdKitModule.get_mol(candForErg);
            refMolForErg = rdKitModule.get_mol(refSmiles);
            if (candMolForErg && refMolForErg) {
              const ergResult = computeErgSharedAtoms(
                refMolForErg, candMolForErg, markedRefAtoms, familyQmols, rdKitModule);
              if (ergResult.candAtoms.size > 0) {
                const highlightAtomColors: {[k: number]: number[]} = {};
                for (const a of ergResult.candAtoms)
                  highlightAtomColors[a] = [...sharedColorRgb];
                const highlightBondColors: {[k: number]: number[]} = {};
                for (const b of ergResult.candBonds)
                  highlightBondColors[b] = [...sharedColorRgb];
                sharedSubstructs.set(rowIdx, {
                  atoms: [...ergResult.candAtoms],
                  bonds: ergResult.candBonds,
                  highlightAtomColors,
                  highlightBondColors,
                });
                const candMolblock = candMolForErg.get_molblock?.() ?? '';
                if (candMolblock) {
                  const frag = extractFragmentSmiles(
                    candMolblock, ergResult.candAtoms, rdKitModule);
                  if (frag) sharedFragments.set(rowIdx, frag);
                }
              }
            }
          }
        } catch {/* ErG failure on this row — leave both views empty */} finally {
          candMolForErg?.delete();
          refMolForErg?.delete();
        }
      }
    } finally {
      mcsMol?.delete();
    }

    if ((k % 20) === 0) {
      const pct = 60 + Math.floor(35 * (k / Math.max(1, topRows.length)));
      progress.update(pct, `MCS ${k + 1}/${topRows.length}...`);
    }
  }
  if (familyQmols) disposeFamilyQmols(familyQmols);
  return {refMcsAtoms: out, sharedSubstructs, sharedFragments};
}
