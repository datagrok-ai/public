/** V2000 molblock manipulation helpers used by scaffold hopping's Local-mode
 *  R-group decomposition and per-row alignment code. RDKit-WASM doesn't
 *  expose direct bond / atom mutation APIs on `RDMol`, so we parse the
 *  molblock text and emit a modified one — cheap (microseconds for
 *  drug-sized molecules) and the V2000 format is fixed-width.
 *
 *  Extracted from `scaffold-hopping.ts` as part of the multi-module split;
 *  no behaviour change. */

import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/** Extracts a sub-fragment of a molecule consisting of the atoms in
 *  `atomsToKeep` (0-indexed) and the bonds between them, returning the
 *  result as a canonical SMILES string. Boundary-atom valences are filled
 *  in by RDKit's implicit-hydrogen perception at parse time, so cuts read
 *  as plain valence completion rather than `[*]` attachment markers —
 *  cleaner for visual rendering, less informative for fragment-set
 *  analyses (callers that need attachment markers can post-process).
 *
 *  Implementation: parse the parent's V2000 molblock as text, drop atom
 *  records for atoms not in `atomsToKeep`, drop bond records that touch
 *  any removed atom, renumber surviving bond endpoints (V2000 is 1-
 *  indexed), re-emit the molblock, and let RDKit parse the result.
 *  Disconnected fragments come back as multi-component SMILES separated
 *  by `.` — informative ("the marked region maps to two pieces in this
 *  candidate") rather than an error.
 *
 *  Returns `null` on parse failure (invalid valence, malformed molblock,
 *  empty keep-set), so callers get a clean miss signal rather than a
 *  bogus partial structure. */
export function extractFragmentSmiles(
  parentMolblock: string, atomsToKeep: Set<number>, rdkit: RDModule,
): string | null {
  if (atomsToKeep.size === 0 || !parentMolblock) return null;
  const lines = parentMolblock.split('\n');
  if (lines.length < 5) return null;

  // V2000 counts line: "%3d%3d ..." — first 3 cols nAtoms, next 3 nBonds.
  const countsLine = lines[3];
  const nAtoms = parseInt(countsLine.substring(0, 3));
  const nBonds = parseInt(countsLine.substring(3, 6));
  if (!Number.isFinite(nAtoms) || !Number.isFinite(nBonds) || nAtoms === 0) return null;

  // Build old (0-indexed) → new (1-indexed) atom-number map for survivors.
  const oldToNew = new Map<number, number>();
  let newIdx = 1;
  for (let i = 0; i < nAtoms; i++) {
    if (atomsToKeep.has(i)) {
      oldToNew.set(i, newIdx);
      newIdx++;
    }
  }
  if (oldToNew.size === 0) return null;

  // Atom block: lines[4 .. 4+nAtoms-1]. Preserve original record verbatim
  // for kept atoms — that retains coords, charges, isotopes, etc.
  const atomLines: string[] = [];
  for (let i = 0; i < nAtoms; i++)
    if (atomsToKeep.has(i)) atomLines.push(lines[4 + i]);

  // Bond block: lines[4+nAtoms .. 4+nAtoms+nBonds-1]. Format
  // `%3d%3d%3d%3d ...` with 1-indexed atom positions; renumber the two
  // endpoints to the new indices, drop bonds that cross the keep boundary.
  const bondLines: string[] = [];
  for (let j = 0; j < nBonds; j++) {
    const bondLine = lines[4 + nAtoms + j];
    if (bondLine.length < 6) continue;
    const a = parseInt(bondLine.substring(0, 3)) - 1;
    const b = parseInt(bondLine.substring(3, 6)) - 1;
    if (!atomsToKeep.has(a) || !atomsToKeep.has(b)) continue;
    const newA = oldToNew.get(a)!;
    const newB = oldToNew.get(b)!;
    bondLines.push(
      newA.toString().padStart(3) + newB.toString().padStart(3) + bondLine.substring(6));
  }

  // Re-emit counts line with updated atom/bond counts; preserve trailing
  // fields (chiralFlag, V2000 marker etc.).
  const newCountsLine =
    atomLines.length.toString().padStart(3) +
    bondLines.length.toString().padStart(3) +
    countsLine.substring(6);

  const newMolblock = [
    lines[0] || '', // title (often empty)
    lines[1] || '  scaffold-hop fragment',
    lines[2] || '', // comment (often empty)
    newCountsLine,
    ...atomLines,
    ...bondLines,
    'M  END',
  ].join('\n');

  let mol: any = null;
  try {
    mol = rdkit.get_mol(newMolblock);
    if (!mol) return null;
    const smi = mol.get_smiles();
    return smi || null;
  } catch {
    return null;
  } finally {
    mol?.delete();
  }
}

/** Returns ring atoms AND the non-bridge edges between them. A bond is
 *  a "ring edge" iff it is part of at least one cycle (= removing it
 *  doesn't disconnect its endpoints). Tested per-bond by BFS from one
 *  endpoint avoiding the bond — if the other endpoint is reached, the
 *  bond is in a cycle.
 *
 *  Both outputs are returned together because `findRingSystems` needs
 *  the ring-edge subgraph (not the full adjacency restricted to ring
 *  atoms) to group atoms correctly: pyrazole-CH2-piperidine has both
 *  pyrazole AND piperidine atoms as ring atoms, but the bond between
 *  them is a BRIDGE — so they're two distinct ring systems. Walking
 *  the full adjacency restricted to ring atoms would incorrectly merge
 *  them.
 *
 *  O(B × (A + B)) where A = atom count, B = bond count. For drug-sized
 *  molecules (~40 atoms / ~45 bonds) this is microseconds — fine for
 *  our use case. A linear-time Tarjan bridge-finding implementation
 *  would be cleaner but the per-bond BFS is dead-simple and avoids
 *  recursion-depth concerns on the WASM-bridged worker. */
function findRingAtomsAndEdges(
  adj: Map<number, number[]>,
): {ringAtoms: Set<number>; ringEdges: Map<number, Set<number>>} {
  const ringAtoms = new Set<number>();
  const ringEdges = new Map<number, Set<number>>();
  for (const a of adj.keys()) ringEdges.set(a, new Set());
  const seen = new Set<string>();
  for (const [u, nbrs] of adj.entries()) {
    for (const v of nbrs) {
      const key = u < v ? `${u}|${v}` : `${v}|${u}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const visited = new Set<number>([u]);
      const queue: number[] = [u];
      let reached = false;
      while (queue.length > 0 && !reached) {
        const x = queue.shift()!;
        for (const y of adj.get(x) ?? []) {
          if ((x === u && y === v) || (x === v && y === u)) continue;
          if (visited.has(y)) continue;
          if (y === v) {reached = true; break;}
          visited.add(y);
          queue.push(y);
        }
      }
      if (reached) {
        ringAtoms.add(u);
        ringAtoms.add(v);
        ringEdges.get(u)?.add(v);
        ringEdges.get(v)?.add(u);
      }
    }
  }
  return {ringAtoms, ringEdges};
}

/** Returns the set of atoms that participate in at least one ring.
 *  Thin wrapper around `findRingAtomsAndEdges` for callers that only
 *  need the atom set. */
export function findRingAtoms(adj: Map<number, number[]>): Set<number> {
  return findRingAtomsAndEdges(adj).ringAtoms;
}

/** Selects the candidate's ring system that occupies the position of
 *  the user's marked region in the reference. Used by the Replaced
 *  Region extraction in BOTH the Local mode and the ErG (Easy/Middle/
 *  Hard) modes so the column shows just the heterocycle that replaced
 *  the marked ring (pyrimidine → pyrazole / thiazole / pyridine), not
 *  the union of every non-MCS-preserved component.
 *
 *  Criterion: a cand ring system is the "replacement" iff it is NOT
 *  entirely inside `candPreserved` AND at least one ring atom has a
 *  neighbour OUTSIDE the ring that is in `edgeAnchors` (= cand atoms
 *  paired to ref atoms that bordered the marked region in ref).
 *
 *  Returns the cand atoms NOT in candPreserved that belong to the
 *  selected ring system(s). Empty set if no qualifying ring exists
 *  (caller falls back to BFS / ErG / etc.). */
export function selectReplacementRingAtoms(
  candAdj: Map<number, number[]>,
  candPreserved: Set<number>,
  edgeAnchors: Set<number>,
): Set<number> {
  const ringSystems = findRingSystems(candAdj);
  const selected = new Set<number>();
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
    if (bondedToEdge) {
      for (const a of ring)
        if (!candPreserved.has(a)) selected.add(a);
    }
  }
  return selected;
}

/** Groups ring atoms into maximal connected components of the
 *  RING-EDGE subgraph — i.e. each returned Set is one ring system (a
 *  single ring OR a set of fused rings sharing atoms, like
 *  pyrrolopyridine's 9 atoms forming one system). Atoms NOT in any
 *  ring are absent from every returned Set. Rings connected only by
 *  bridge bonds (e.g. biphenyl's two phenyls, or DLK compound 3's
 *  pyrazole–piperidine pair) are returned as SEPARATE systems because
 *  the connecting bond is a bridge — by design, since chemically those
 *  are two distinct ring systems linked by a rotatable bond.
 *
 *  Critically, we walk over `ringEdges` (non-bridge bonds) NOT the
 *  full adjacency. Walking over full adjacency would incorrectly merge
 *  bridge-connected ring systems into a single component because the
 *  bridge bond's endpoints are BOTH ring atoms (each individually in
 *  a cycle, just not in the SAME cycle through the bridge).
 *
 *  Used by the Local-mode "Replaced Region" extraction to surface the
 *  candidate's ring system that occupies the position of the user's
 *  marked ring in the reference — bypassing unrelated substituents the
 *  BFS-on-non-preserved-atoms approach would also pick up. */
export function findRingSystems(adj: Map<number, number[]>): Set<number>[] {
  const {ringAtoms, ringEdges} = findRingAtomsAndEdges(adj);
  const systems: Set<number>[] = [];
  const visited = new Set<number>();
  for (const a of ringAtoms) {
    if (visited.has(a)) continue;
    const system = new Set<number>();
    const stack: number[] = [a];
    while (stack.length > 0) {
      const x = stack.pop()!;
      if (visited.has(x)) continue;
      visited.add(x);
      system.add(x);
      for (const y of ringEdges.get(x) ?? [])
        if (!visited.has(y)) stack.push(y);
    }
    systems.push(system);
  }
  return systems;
}

/** Parses a V2000 molblock's bond block into an undirected atom-atom
 *  adjacency map (atom indices are 0-based). RDKit-WASM doesn't expose
 *  a direct bond iterator on `RDMol`, so we parse the molblock text
 *  ourselves — cheap (~microseconds for drug-sized molecules) and the
 *  V2000 bond record format is fixed-width and easy to slice. */
export function parseV2000Adjacency(molblock: string): Map<number, number[]> {
  const adj = new Map<number, number[]>();
  if (!molblock) return adj;
  const lines = molblock.split('\n');
  if (lines.length < 5) return adj;
  const counts = lines[3];
  if (!counts) return adj;
  const nA = parseInt(counts.substring(0, 3));
  const nB = parseInt(counts.substring(3, 6));
  if (!Number.isFinite(nA) || !Number.isFinite(nB)) return adj;
  for (let i = 0; i < nA; i++) adj.set(i, []);
  for (let j = 0; j < nB; j++) {
    const bondLine = lines[4 + nA + j];
    if (!bondLine || bondLine.length < 6) continue;
    const a = parseInt(bondLine.substring(0, 3)) - 1;
    const b = parseInt(bondLine.substring(3, 6)) - 1;
    if (a >= 0 && b >= 0) {
      adj.get(a)?.push(b);
      adj.get(b)?.push(a);
    }
  }
  return adj;
}

/** Sibling to `extractFragmentSmiles` that ALSO adds explicit R-group
 *  labels (R1, R2, …) at every boundary bond (one endpoint kept, the
 *  other dropped). Returns a V2000 molblock string with `M RGP` records
 *  marking each dummy atom as an R-group, so RDKit's renderer displays
 *  them as `R1`, `R2`, … instead of `*:1` / `*:2`.
 *
 *  Why molblock and not SMILES: SMILES has no R-group concept — the
 *  closest analog is `[*:1]` (atom-map numbering) which renders as
 *  literal `*:1` in the grid. V2000's `M RGP` record is the standard
 *  way to say "this dummy atom is R-group N" and the Chem cell renderer
 *  picks it up automatically (same mechanism the existing MMP "from" /
 *  "to" columns use).
 *
 *  Without the attachment markers, single-atom replacements render as
 *  "CH4" or "H2O" — chemically meaningless. With them, the same atom
 *  becomes e.g. `C-R1` — a methyl substituent at R-group position 1 —
 *  which reads correctly as "a methyl that connects to the rest of the
 *  candidate at this position." Used by the Local-mode replacement
 *  path AND by the reference-row population so the column reads
 *  top-to-bottom as "this is what I marked → here's each candidate's
 *  version of it." */
export function extractFragmentMolblockWithRGroups(
  parentMolblock: string, atomsToKeep: Set<number>, rdkit: RDModule,
): string | null {
  if (atomsToKeep.size === 0 || !parentMolblock) return null;
  const lines = parentMolblock.split('\n');
  if (lines.length < 5) return null;
  const countsLine = lines[3];
  const nAtoms = parseInt(countsLine.substring(0, 3));
  const nBonds = parseInt(countsLine.substring(3, 6));
  if (!Number.isFinite(nAtoms) || !Number.isFinite(nBonds) || nAtoms === 0) return null;

  // 1-indexed renumbering for the kept atoms.
  const oldToNew = new Map<number, number>();
  let newIdx = 1;
  for (let i = 0; i < nAtoms; i++)
    if (atomsToKeep.has(i)) {oldToNew.set(i, newIdx); newIdx++;}

  if (oldToNew.size === 0) return null;

  const realAtomLines: string[] = [];
  for (let i = 0; i < nAtoms; i++) if (atomsToKeep.has(i)) realAtomLines.push(lines[4 + i]);
  const dummyAtomLines: string[] = [];

  // dummyAtomIndices[i] = 1-indexed molblock position of the i-th dummy
  // atom we add. Used to build the M RGP record at the end.
  const dummyAtomIndices: number[] = [];

  const bondLines: string[] = [];
  for (let j = 0; j < nBonds; j++) {
    const bondLine = lines[4 + nAtoms + j];
    if (!bondLine || bondLine.length < 6) continue;
    const a = parseInt(bondLine.substring(0, 3)) - 1;
    const b = parseInt(bondLine.substring(3, 6)) - 1;
    const aKept = atomsToKeep.has(a);
    const bKept = atomsToKeep.has(b);
    if (aKept && bKept) {
      const newA = oldToNew.get(a)!;
      const newB = oldToNew.get(b)!;
      bondLines.push(
        newA.toString().padStart(3) + newB.toString().padStart(3) + bondLine.substring(6));
    } else if (aKept || bKept) {
      const keptEnd = aKept ? a : b;
      const newKeptEnd = oldToNew.get(keptEnd)!;
      // V2000 atom record — `R#` symbol marks this as an R-group
      // placeholder. The actual R-group label is set via the M RGP
      // record at the bottom. Using `R#` plus an RGP record (rather
      // than just `*`) is what makes the renderer show `R1` / `R2`.
      dummyAtomLines.push('    0.0000    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0');
      const newDummyIdx = realAtomLines.length + dummyAtomLines.length;
      dummyAtomIndices.push(newDummyIdx);
      bondLines.push(
        newKeptEnd.toString().padStart(3) + newDummyIdx.toString().padStart(3) + '  1  0  0  0  0');
    }
  }

  const atomLines = realAtomLines.concat(dummyAtomLines);
  const newCountsLine =
    atomLines.length.toString().padStart(3) +
    bondLines.length.toString().padStart(3) +
    countsLine.substring(6);

  // Build the M RGP line. V2000 packs up to 8 (atomIdx, rgrpLabel)
  // pairs per line, each field 4 chars wide, right-justified.
  const rgpLines: string[] = [];
  if (dummyAtomIndices.length > 0) {
    for (let i = 0; i < dummyAtomIndices.length; i += 8) {
      const chunk = dummyAtomIndices.slice(i, i + 8);
      let line = `M  RGP${chunk.length.toString().padStart(3)}`;
      for (let j = 0; j < chunk.length; j++) {
        const idx = chunk[j];
        const label = i + j + 1; // R1, R2, …
        line += idx.toString().padStart(4) + label.toString().padStart(4);
      }
      rgpLines.push(line);
    }
  }

  const newMolblock = [
    lines[0] || '',
    lines[1] || '  scaffold-hop fragment',
    lines[2] || '',
    newCountsLine,
    ...atomLines,
    ...bondLines,
    ...rgpLines,
    'M  END',
  ].join('\n');

  // Round-trip through RDKit so aromaticity / valences get re-perceived
  // properly, then regenerate 2D coordinates so the R# atoms sit at a
  // sensible distance from their bonded neighbours. Without the
  // coord-regen step, the dummy atoms keep their (0, 0, 0) placeholder
  // coords and the renderer draws very long bonds from the ring to R1.
  // `useCoordGen: true` picks the CoordGen layout engine, which gives
  // cleaner aromatic-ring + substituent placement than the default
  // RDKit 2D depictor.
  let mol: any = null;
  try {
    mol = rdkit.get_mol(newMolblock);
    if (!mol) return newMolblock;
    try {mol.set_new_coords?.(true);} catch {/* layout engine missing — keep raw coords */}
    const canon = mol.get_molblock?.() ?? '';
    return canon || newMolblock;
  } catch {
    return newMolblock;
  } finally {
    mol?.delete();
  }
}


/** Builds a CLEAN fragment molblock (no R-group markers, just the kept
 *  atoms + bonds-between-them) for use as an alignment scaffold. Datagrok's
 *  cell renderer reads the `.%chem-scaffold-align` column tag, parses
 *  the value as a molblock, and calls `mol.generate_aligned_coords(scaffold)`
 *  on every cell to align them to a common orientation. Boundary atoms
 *  are H-filled by RDKit on parse — that's fine for substructure-matching
 *  during alignment, the H's don't interfere. */
export function buildCleanFragmentMolblock(
  parentMolblock: string, atomsToKeep: Set<number>, rdkit: RDModule,
): string | null {
  const smi = extractFragmentSmiles(parentMolblock, atomsToKeep, rdkit);
  if (!smi) return null;
  let mol: any = null;
  try {
    mol = rdkit.get_mol(smi);
    if (!mol) return null;
    try {mol.set_new_coords?.(true);} catch {/* layout engine missing */}
    const molblock = mol.get_molblock?.() ?? '';
    return molblock || null;
  } catch {
    return null;
  } finally {mol?.delete();}
}

/** Builds a 3-atom alignment anchor from a fragment molblock that has
 *  R# dummy atoms: R# + the atom it's bonded to + ONE of that atom's
 *  other neighbours. Three atoms define a triangle, which fixes both
 *  position AND rotation — a 2-atom anchor only fixes direction, so
 *  every fragment could still flip 180° around the anchor bond. The
 *  third atom locks the rotation: rendered fragments come out with R1
 *  in the same direction AND the ring extending the same way.
 *
 *  RDKit's substructure search with `allowRGroups: true` treats R# as
 *  a wildcard, so this anchor matches any fragment with an `R-X-Y`
 *  pattern where X is the right element (typically N for piperazine /
 *  piperidine / morpholine variants). The third atom Y is the X's
 *  other ring neighbour, so it pretty much always matches a C in any
 *  candidate fragment. */
export function buildR1AnchorMolblock(
  fragmentMolblock: string,
): string | null {
  if (!fragmentMolblock) return null;
  const lines = fragmentMolblock.split('\n');
  if (lines.length < 5) return null;
  const counts = lines[3];
  if (!counts) return null;
  const nA = parseInt(counts.substring(0, 3));
  const nB = parseInt(counts.substring(3, 6));
  if (!Number.isFinite(nA) || !Number.isFinite(nB)) return null;

  // Find the R# atom (V2000 atom block, cols 32-34 hold the symbol).
  let rIdx = -1;
  for (let i = 0; i < nA; i++) {
    const line = lines[4 + i];
    if (!line || line.length < 34) continue;
    const sym = line.substring(31, 34).trim();
    if (sym === 'R#' || sym === '*' || sym.startsWith('R')) {rIdx = i; break;}
  }
  if (rIdx < 0) return null;

  // Build adjacency (atom → list of (neighbour, bondLine)) so we can
  // find R#'s neighbour and that neighbour's other neighbours.
  const adj = new Map<number, Array<{neighbour: number; bondLine: string}>>();
  for (let i = 0; i < nA; i++) adj.set(i, []);
  for (let j = 0; j < nB; j++) {
    const line = lines[4 + nA + j];
    if (!line || line.length < 6) continue;
    const a = parseInt(line.substring(0, 3)) - 1;
    const b = parseInt(line.substring(3, 6)) - 1;
    adj.get(a)?.push({neighbour: b, bondLine: line});
    adj.get(b)?.push({neighbour: a, bondLine: line});
  }

  const rNeighbours = adj.get(rIdx) ?? [];
  if (rNeighbours.length === 0) return null;
  const nA1 = rNeighbours[0].neighbour; // R#'s sole neighbour
  const rToN1Bond = rNeighbours[0].bondLine;

  // Pick A1's OTHER neighbours (not the R atom). For a typical ring
  // attachment point, A1 has two or three connections — the R atom
  // itself, ring-C-left, and (if present) ring-C-right.
  // We pick the first non-R neighbour.
  const a1Neighbours = (adj.get(nA1) ?? []).filter((x) => x.neighbour !== rIdx);
  if (a1Neighbours.length === 0) {
    // Pathological — only 2 atoms total. Fall back to the 2-atom anchor.
    const aOrig = parseInt(rToN1Bond.substring(0, 3)) - 1;
    const newA = aOrig === rIdx ? 1 : 2;
    const newB = aOrig === rIdx ? 2 : 1;
    const newBondLine = newA.toString().padStart(3) +
      newB.toString().padStart(3) + rToN1Bond.substring(6);
    const newCounts = '  2  1' + counts.substring(6);
    return ['', '  scaffold-hop anchor', '', newCounts,
      lines[4 + rIdx], lines[4 + nA1], newBondLine,
      'M  RGP  1   1   1', 'M  END'].join('\n');
  }
  const nA2 = a1Neighbours[0].neighbour;
  const n1ToN2Bond = a1Neighbours[0].bondLine;

  // Emit a 3-atom molblock. Renumber the three atoms 1..3:
  //   1 = R#, 2 = A1 (R#'s neighbour), 3 = A2 (A1's other neighbour).
  // Bond 1: R# (atom 1) ↔ A1 (atom 2). Bond 2: A1 (atom 2) ↔ A2 (atom 3).
  // The M RGP record marks atom 1 as R1.
  const r1AtomLine = lines[4 + rIdx];
  const a1AtomLine = lines[4 + nA1];
  const a2AtomLine = lines[4 + nA2];
  const rToN1Orig = parseInt(rToN1Bond.substring(0, 3)) - 1;
  const newBond1 = (rToN1Orig === rIdx ? '  1  2' : '  2  1') + rToN1Bond.substring(6);
  const n1ToN2Orig = parseInt(n1ToN2Bond.substring(0, 3)) - 1;
  const newBond2 = (n1ToN2Orig === nA1 ? '  2  3' : '  3  2') + n1ToN2Bond.substring(6);
  const newCounts = '  3  2' + counts.substring(6);
  return [
    '',
    '  scaffold-hop anchor',
    '',
    newCounts,
    r1AtomLine,
    a1AtomLine,
    a2AtomLine,
    newBond1,
    newBond2,
    'M  RGP  1   1   1',
    'M  END',
  ].join('\n');
}
