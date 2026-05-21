import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {findRingSystems as findBridgeAwareRingSystems} from './scaffold-hopping-molblock';

/** ErG (Stiefl 2006) reduced-graph scaffold matching for the scaffold-
 *  hopping result table.
 *
 *  Why this exists: strict MCS (atom-element identity, connected match)
 *  treats pyrimidine ↔ thiazole as different rings and refuses to span
 *  ring matches across non-matching linker atoms. ErG-style reduced
 *  graphs collapse each ring system to a single abstract node labelled
 *  by its pharmacophore types — so pyrimidine [aromatic, acceptor×2]
 *  and thiazole [aromatic, acceptor] both match imatinib's pyridyl
 *  [aromatic, acceptor], capturing the kinase-hinge-binder equivalence
 *  that strict MCS misses.
 *
 *  Reuses the same 7-family SMARTS the rest of Chem uses (Pharmacophore
 *  Features panel + CATS descriptor), so atom typing stays consistent.
 *
 *  Reference: Stiefl, Watson, Baumann, Mika 2006. "ErG: 2D Pharmacophore
 *  Descriptions for Scaffold Hopping." JCIM 46:208-220.
 *
 *  Implementation chooses the *disconnected* relaxation of ErG-MCS: each
 *  reference node is greedily matched with the best-scoring candidate
 *  node by pharmacophore-label overlap, with no connectedness constraint.
 *  This matters because in real scaffold hops the equivalent rings are
 *  often separated by different linker chains (imatinib's tolyl-NH-
 *  pyrimidine vs. dasatinib's chloroaryl-NH-amide-thiazole-NH-pyrimidine
 *  for instance) — strict connected matching would miss the pyrimidine
 *  match. We trade some specificity for correct chemistry. */

/** Pharmacophore family letter → human-readable name. Mirrors the table
 *  in `panels/pharmacophore-features.ts` and the CATS Python script —
 *  single source of truth for atom typing across the Chem package. */
export const ERG_FAMILY_LETTER_TO_NAME: Record<string, string> = {
  'D': 'Donor',
  'A': 'Acceptor',
  'H': 'Hydrophobic',
  'a': 'Aromatic',
  'P': 'Positive',
  'N': 'Negative',
  'X': 'HalogenBond',
};

/** A node in the ErG reduced graph. Either a ring system (multiple atoms)
 *  or a single chain atom. */
export interface ErgNode {
  /** Original-mol atom indices contained in this node. For rings: all
   *  atoms in the ring system. For chain atoms: a singleton. */
  atoms: number[];
  /** Pharmacophore labels (family names) aggregated over `atoms`. */
  pharma: Set<string>;
  /** True if this node abstracts a ring system. */
  isRing: boolean;
  /** Number of atoms in this node — used as a tiebreaker for ring-size
   *  preferences in matching (e.g. prefer 6-membered ↔ 6-membered when
   *  multiple candidates have matching pharma labels). */
  size: number;
}

/** Compiled query mols for the 7 pharmacophore families. Multiple SMARTS
 *  per family are unioned at match time (any match → atom belongs). */
export type FamilyQmols = Map<string, any[]>;

/** Builds the map of compiled family-SMARTS query mols from the
 *  pharmacophore-features dataframe (the same `pharmacophore-features.csv`
 *  the CATS step loads). Caller is responsible for `.delete()`-ing the
 *  returned qmols when done — `disposeFamilyQmols` is the symmetric
 *  cleanup. */
export function buildFamilyQmols(
  featuresDf: any, rdkit: RDModule,
): FamilyQmols {
  const out: FamilyQmols = new Map();
  for (const name of Object.values(ERG_FAMILY_LETTER_TO_NAME))
    out.set(name, []);
  if (!featuresDf?.rowCount) return out;
  const familyCol = featuresDf.col('family');
  const smartsCol = featuresDf.col('smarts');
  if (!familyCol || !smartsCol) return out;
  for (let i = 0; i < featuresDf.rowCount; i++) {
    const letter: string = familyCol.get(i);
    const smarts: string = smartsCol.get(i);
    if (!letter || !smarts) continue;
    const familyName = ERG_FAMILY_LETTER_TO_NAME[letter];
    if (!familyName) continue;
    try {
      const qmol = rdkit.get_qmol(smarts);
      if (qmol) out.get(familyName)!.push(qmol);
    } catch {/* skip unparseable SMARTS */}
  }
  return out;
}

export function disposeFamilyQmols(qmols: FamilyQmols): void {
  for (const arr of qmols.values()) for (const q of arr) q?.delete?.();
}

/** Tags every atom of `mol` with the set of pharmacophore family names
 *  that it belongs to. Returns `Map<atomIdx, Set<familyName>>`; atoms
 *  matching no family get an empty set (still present in the map for
 *  uniform iteration). */
export function labelAtomsByPharmacophore(
  mol: any, familyQmols: FamilyQmols,
): Map<number, Set<string>> {
  const labels = new Map<number, Set<string>>();
  const nAtoms = mol.get_num_atoms();
  for (let i = 0; i < nAtoms; i++) labels.set(i, new Set());
  for (const [familyName, qmols] of familyQmols.entries()) {
    for (const q of qmols) {
      try {
        const matchesJson = mol.get_substruct_matches(q);
        if (!matchesJson) continue;
        const matches = JSON.parse(matchesJson);
        if (!Array.isArray(matches)) continue;
        for (const match of matches) {
          const atomIdxs: number[] = match?.atoms ?? [];
          for (const a of atomIdxs) labels.get(a)!.add(familyName);
        }
      } catch {/* SMARTS match failure on individual family — skip */}
    }
  }
  return labels;
}

/** Identifies ring systems in `mol`. A ring system is a connected
 *  component in the **ring-edge** subgraph — fused rings (sharing an
 *  edge) collapse into one system, spiro systems likewise. Rings
 *  connected only by a **bridge bond** (e.g. a biaryl, or DLK
 *  compound 3's pyrazole-piperidine where the C-C bond between the
 *  two ring atoms is a rotatable single bond, NOT in any ring) are
 *  returned as SEPARATE systems — chemically correct since they're
 *  two distinct rings linked by a rotatable bond.
 *
 *  Delegates to `findBridgeAwareRingSystems` from `scaffold-hopping-
 *  molblock.ts`, which walks ring-edges only (not full adjacency
 *  restricted to ring atoms). Walking the full restricted adjacency
 *  is the older approach this function used to do — it INCORRECTLY
 *  merges bridge-connected rings into a single component because
 *  both endpoints of the bridge ARE ring atoms (each individually in
 *  a cycle, just not in the same cycle through the bridge bond).
 *  That over-merging propagated into ErG's reduced-graph: pyrazole +
 *  piperidine became one "ring system" node, the ErG match treated
 *  the combined system as the candidate's image of ref's marked
 *  pyrimidine, and the Replaced Region surfaced both rings instead
 *  of just the heterocycle. The bridge-aware version fixes this.
 *
 *  Returns an array of atom-index arrays, one per ring system. */
export function findRingSystems(mol: any, _rdkit: RDModule): number[][] {
  const molblock = mol.get_molblock?.() ?? '';
  if (!molblock) return [];
  const bonds = parseV2000Bonds(molblock);
  if (bonds.length === 0) return [];
  // Build an undirected adjacency map covering EVERY atom, not just ring
  // atoms — the bridge-detection inside `findBridgeAwareRingSystems`
  // needs the full graph to correctly classify each bond.
  const nAtoms = mol.get_num_atoms?.() ?? 0;
  const adj = new Map<number, number[]>();
  for (let i = 0; i < nAtoms; i++) adj.set(i, []);
  for (const [a, b] of bonds) {
    adj.get(a)?.push(b);
    adj.get(b)?.push(a);
  }
  const systems = findBridgeAwareRingSystems(adj);
  return systems.map((s) => Array.from(s));
}

/** Parses the V2000 molblock bond block into [atom_a (0-indexed),
 *  atom_b (0-indexed)] pairs. Robust to leading/trailing whitespace
 *  on lines; returns empty list if the molblock is malformed. */
function parseV2000Bonds(molblock: string): Array<[number, number]> {
  if (!molblock) return [];
  const lines = molblock.split('\n');
  if (lines.length < 5) return [];
  const counts = lines[3];
  const nAtoms = parseInt(counts.substring(0, 3));
  const nBonds = parseInt(counts.substring(3, 6));
  if (!Number.isFinite(nAtoms) || !Number.isFinite(nBonds)) return [];
  const out: Array<[number, number]> = [];
  for (let j = 0; j < nBonds; j++) {
    const line = lines[4 + nAtoms + j];
    if (!line || line.length < 6) continue;
    const a = parseInt(line.substring(0, 3));
    const b = parseInt(line.substring(3, 6));
    if (Number.isFinite(a) && Number.isFinite(b))
      out.push([a - 1, b - 1]);
  }
  return out;
}

/** Builds the ErG reduced graph for `mol`. */
export function buildReducedGraph(
  mol: any, ringSystems: number[][], atomLabels: Map<number, Set<string>>,
  _rdkit: RDModule,
): {nodes: ErgNode[]; adj: Set<number>[]} {
  const nAtoms = mol.get_num_atoms();
  const ringAtomToSystem = new Map<number, number>();
  ringSystems.forEach((sys, sysIdx) => {
    for (const a of sys) ringAtomToSystem.set(a, sysIdx);
  });

  // Each ring system → one node; each non-ring atom → one node.
  const nodes: ErgNode[] = [];
  const atomToNode = new Map<number, number>();
  for (let s = 0; s < ringSystems.length; s++) {
    const sysAtoms = ringSystems[s];
    const pharma = new Set<string>();
    for (const a of sysAtoms) for (const p of atomLabels.get(a) ?? []) pharma.add(p);
    const nodeIdx = nodes.length;
    nodes.push({atoms: sysAtoms.slice(), pharma, isRing: true, size: sysAtoms.length});
    for (const a of sysAtoms) atomToNode.set(a, nodeIdx);
  }
  for (let a = 0; a < nAtoms; a++) {
    if (atomToNode.has(a)) continue;
    const nodeIdx = nodes.length;
    nodes.push({atoms: [a], pharma: new Set(atomLabels.get(a) ?? []), isRing: false, size: 1});
    atomToNode.set(a, nodeIdx);
  }

  // Build adjacency from the molblock bonds: each bond becomes an edge
  // between the two nodes its endpoints belong to (skip self-loops on a
  // ring-internal bond).
  const adj: Set<number>[] = nodes.map(() => new Set<number>());
  const molblock = mol.get_molblock?.() ?? '';
  for (const [a, b] of parseV2000Bonds(molblock)) {
    const na = atomToNode.get(a);
    const nb = atomToNode.get(b);
    if (na == null || nb == null || na === nb) continue;
    adj[na].add(nb);
    adj[nb].add(na);
  }
  return {nodes, adj};
}

/** Compatibility score between two ErG nodes — Jaccard over pharmacophore
 *  label sets, with `isRing` parity and size proximity as tiebreakers.
 *  Returns 0 when `isRing` flags differ (a ring node never matches a
 *  chain atom). Higher = better. */
export function nodeCompatibility(a: ErgNode, b: ErgNode): number {
  if (a.isRing !== b.isRing) return 0;
  // Jaccard on pharma labels. Empty ∩ empty → return 0 (no match): two
  // chain atoms with no pharmacophore tags (typical CH2/CH3 linker
  // atoms) carry no positive evidence that they should be considered
  // equivalent. Earlier this branch returned 0.05 to allow "neutral
  // linker atoms" to pair with each other, but that 0.05 happens to
  // equal the match threshold in `matchReducedGraphs`, so EVERY pair
  // of unlabeled chain atoms got matched — inflating the Replacement-
  // column highlights with chemically meaningless chain-atom matches
  // on any drug-like molecule with a propyl/butyl linker.
  let inter = 0;
  for (const p of a.pharma) if (b.pharma.has(p)) inter++;
  const union = a.pharma.size + b.pharma.size - inter;
  let score = union === 0 ? 0 : inter / union;
  // Ring-size proximity bonus (small, just a tiebreaker).
  if (a.isRing && b.isRing && a.size > 0 && b.size > 0)
    score += 0.05 * (1 - Math.abs(a.size - b.size) / Math.max(a.size, b.size));
  return score;
}

/** Disconnected greedy match between two reduced graphs. For each ref
 *  node, picks the best-scoring still-unmatched cand node (compatibility
 *  ≥ threshold). Ref nodes with no compatible candidate stay unmatched.
 *
 *  Two-tier ordering when `markedRefAtoms` is provided: ref nodes that
 *  do NOT contain marked atoms (= the "conserved core") are matched
 *  FIRST, in score-desc order. Then ref nodes that DO contain marked
 *  atoms (= the user's region of interest) match against whatever cand
 *  nodes are left. This is load-bearing for the Replaced Region
 *  extraction: without it, a marked ref ring can greedy-steal the
 *  candidate node that should have paired with a non-marked ref ring,
 *  surfacing the wrong cand ring as the "replacement". E.g. DLK ref
 *  pyrimidine (marked, a + 2A, size 6) and cand aminopyridyl (a + 1A,
 *  size 6) have the SAME size and similar pharmacophore signature, so
 *  the unconstrained matcher pairs them — but cand aminopyridyl is the
 *  candidate's image of ref aminopyridyl, NOT of ref pyrimidine. The
 *  two-tier order ensures ref aminopyridyl claims cand aminopyridyl
 *  first, leaving ref pyrimidine to match the actual replacement ring
 *  (cand pyrazole / thiazole / etc.).
 *
 *  When `markedRefAtoms` is empty (global-hop run), the order
 *  collapses to pure score-desc — same behaviour as before.
 *
 *  Returns `Map<refNodeIdx, candNodeIdx>`. */
export function matchReducedGraphs(
  refRG: {nodes: ErgNode[]; adj: Set<number>[]},
  candRG: {nodes: ErgNode[]; adj: Set<number>[]},
  threshold: number = 0.05,
  markedRefAtoms: Set<number> = new Set(),
): Map<number, number> {
  type Triple = {r: number; c: number; s: number; markedRef: boolean};
  const isRefMarked = (node: ErgNode): boolean => {
    if (markedRefAtoms.size === 0) return false;
    for (const a of node.atoms) if (markedRefAtoms.has(a)) return true;
    return false;
  };
  const triples: Triple[] = [];
  for (let r = 0; r < refRG.nodes.length; r++) {
    const markedRef = isRefMarked(refRG.nodes[r]);
    for (let c = 0; c < candRG.nodes.length; c++) {
      const s = nodeCompatibility(refRG.nodes[r], candRG.nodes[c]);
      if (s >= threshold) triples.push({r, c, s, markedRef});
    }
  }
  // Sort: non-marked ref nodes first (so they claim cand nodes first),
  // then marked ref nodes. Within each tier, score-desc.
  triples.sort((x, y) => {
    if (x.markedRef !== y.markedRef) return x.markedRef ? 1 : -1;
    return y.s - x.s;
  });
  const refUsed = new Set<number>();
  const candUsed = new Set<number>();
  const out = new Map<number, number>();
  for (const t of triples) {
    if (refUsed.has(t.r) || candUsed.has(t.c)) continue;
    refUsed.add(t.r);
    candUsed.add(t.c);
    out.set(t.r, t.c);
  }
  return out;
}

/** Builds an undirected adjacency list from a list of bond pairs.
 *  `nAtoms` upper-bounds the index space; bonds referencing higher
 *  indices (malformed molblock) are silently dropped. */
function buildAdjacency(
  bondPairs: Array<[number, number]>, nAtoms: number,
): Set<number>[] {
  const adj: Set<number>[] = [];
  for (let i = 0; i < nAtoms; i++) adj.push(new Set());
  for (const [a, b] of bondPairs) {
    if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms) continue;
    adj[a].add(b);
    adj[b].add(a);
  }
  return adj;
}

/** Connected components of an atom set, where two atoms are in the same
 *  component iff there's a path between them using only bonds whose
 *  endpoints are *both* in `atoms`. */
function findConnectedComponents(
  atoms: Set<number>, adj: Set<number>[],
): number[][] {
  const visited = new Set<number>();
  const components: number[][] = [];
  for (const start of atoms) {
    if (visited.has(start)) continue;
    const stack: number[] = [start];
    const component: number[] = [];
    while (stack.length > 0) {
      const x = stack.pop()!;
      if (visited.has(x)) continue;
      visited.add(x);
      component.push(x);
      for (const y of adj[x] ?? [])
        if (atoms.has(y) && !visited.has(y)) stack.push(y);
    }
    components.push(component);
  }
  return components;
}

/** BFS shortest path from any atom in `fromSet` to any atom in `toSet`.
 *  Returns the path as an array of atom indices (including both
 *  endpoints), or `null` if `toSet` is unreachable from `fromSet` —
 *  which is the correct signal for genuinely-disconnected molecule
 *  pieces (salts, co-crystals). All atoms in `fromSet` are seeded as
 *  distance-0; the BFS expands outward and stops at the first atom in
 *  `toSet` that's not already in `fromSet`. */
function shortestPathBetweenSets(
  fromSet: Set<number>, toSet: Set<number>, adj: Set<number>[],
): number[] | null {
  if (fromSet.size === 0 || toSet.size === 0) return null;
  // Already overlapping → single-atom "path" (no extension needed).
  for (const a of fromSet) if (toSet.has(a)) return [a];

  const parent = new Map<number, number>();
  const queue: number[] = [];
  for (const a of fromSet) {
    parent.set(a, -1);
    queue.push(a);
  }

  let target = -1;
  let head = 0;
  while (head < queue.length) {
    const x = queue[head++];
    if (toSet.has(x) && !fromSet.has(x)) {target = x; break;}
    for (const y of adj[x] ?? []) {
      if (parent.has(y)) continue;
      parent.set(y, x);
      queue.push(y);
    }
  }
  if (target < 0) return null;

  // Reconstruct path back to the source frontier.
  const path: number[] = [];
  let cur: number | undefined = target;
  while (cur !== undefined && cur !== -1) {
    path.push(cur);
    cur = parent.get(cur);
  }
  return path.reverse();
}

/** Path-completion pass: iteratively merges connected components within
 *  `seedAtoms` by adding the atoms on the shortest path between the
 *  closest pair of components. Stops when only one component remains
 *  OR when no path connects any pair (genuinely disconnected — e.g.
 *  salt / co-crystal in the candidate; correct fallback to keep the
 *  fragment SMILES multi-component). Mirrors Bemis-Murcko scaffold
 *  semantics: ring systems plus the linker atoms that connect them. */
function connectMatchedAtomsViaShortestPaths(
  seedAtoms: Set<number>, adj: Set<number>[],
): Set<number> {
  const expanded = new Set(seedAtoms);
  // Cap the loop at |seedAtoms| iterations — each iteration merges at
  // least two components, so the loop is bounded even in degenerate cases.
  for (let iter = 0; iter < seedAtoms.size; iter++) {
    const components = findConnectedComponents(expanded, adj);
    if (components.length <= 1) break;

    // Find the closest pair of components.
    let bestPath: number[] | null = null;
    for (let i = 0; i < components.length; i++) {
      const fromSet = new Set(components[i]);
      for (let j = i + 1; j < components.length; j++) {
        const toSet = new Set(components[j]);
        const path = shortestPathBetweenSets(fromSet, toSet, adj);
        if (path && (bestPath === null || path.length < bestPath.length))
          bestPath = path;
      }
    }
    if (!bestPath) break; // genuinely disconnected — leave as separate fragments
    for (const a of bestPath) expanded.add(a);
  }
  return expanded;
}

/** Top-level ErG scaffold-match: builds reduced graphs for both molecules,
 *  runs the disconnected matcher, and resolves the matched node-pairs back
 *  to original-mol atom indices + the bonds between those atoms (used by
 *  the cell renderer to colour both atoms and bonds in the highlight).
 *  The marked-region mode restricts the return to ref atoms whose
 *  containing reduced-graph node has a marked atom in it (and the
 *  candidate-side image of those nodes); the no-marks mode returns the
 *  union of all matched nodes' atoms.
 *
 *  Returns `{refAtoms, candAtoms, candBonds, refRingMatches}`:
 *  - `refAtoms` / `candAtoms`: the matched atom sets, suitable for
 *    rendering the Shared Region fragment + the highlighted-on-full-
 *    molecule view.
 *  - `candBonds`: candidate-mol bond indices (V2000 order) where BOTH
 *    endpoints are in `candAtoms` — feeds the cell renderer's
 *    `highlightBondColors` so the green tint covers entire ring/chain
 *    segments rather than just isolated atoms.
 *  - `refRingMatches`: count of matched RING node pairs (used by the
 *    caller as a "did the ErG match find anything substantial" check). */
export function computeErgSharedAtoms(
  refMol: any, candMol: any, markedRefAtoms: Set<number>,
  familyQmols: FamilyQmols, rdkit: RDModule,
): {
  refAtoms: Set<number>; candAtoms: Set<number>;
  candBonds: number[]; refRingMatches: number;
} {
  const refRings = findRingSystems(refMol, rdkit);
  const candRings = findRingSystems(candMol, rdkit);
  const refLabels = labelAtomsByPharmacophore(refMol, familyQmols);
  const candLabels = labelAtomsByPharmacophore(candMol, familyQmols);
  const refRG = buildReducedGraph(refMol, refRings, refLabels, rdkit);
  const candRG = buildReducedGraph(candMol, candRings, candLabels, rdkit);
  // Pass `markedRefAtoms` so the matcher orders non-marked ref nodes
  // first — see the comment block on `matchReducedGraphs` for why this
  // matters for the Replaced Region extraction.
  const matching = matchReducedGraphs(refRG, candRG, 0.05, markedRefAtoms);

  const refAtoms = new Set<number>();
  const seedCandAtoms = new Set<number>();
  let refRingMatches = 0;
  for (const [r, c] of matching.entries()) {
    const refNode = refRG.nodes[r];
    const candNode = candRG.nodes[c];
    let include = false;
    if (markedRefAtoms.size === 0)
      include = true;
    else {
      for (const a of refNode.atoms)
        if (markedRefAtoms.has(a)) {include = true; break;}
    }
    if (!include) continue;
    for (const a of refNode.atoms) refAtoms.add(a);
    for (const a of candNode.atoms) seedCandAtoms.add(a);
    if (refNode.isRing && candNode.isRing) refRingMatches++;
  }

  // Path-completion: extend the candidate-side atom set with the shortest
  // paths through the candidate's molecular graph that connect each pair
  // of matched ring/chain components. Mirrors Bemis-Murcko scaffold
  // semantics — the matched scaffold + the linker atoms between the
  // matched rings, rendered as a single connected fragment. If two
  // matched rings live in genuinely disconnected components of the
  // candidate (salt, co-crystal), no path exists and the algorithm
  // leaves them disconnected — correct behaviour, the fragment SMILES
  // ends up `.`-separated. Same for the unmatched-NH-bridge case where
  // a path *does* exist: the linker atoms come along even though they
  // didn't have a direct correspondent in the reference, because the
  // visual answer "what does the corresponding scaffold piece look like"
  // requires connectivity, not literal atom-by-atom mapping.
  const candMolblock = candMol.get_molblock?.() ?? '';
  const candNumAtoms = candMol.get_num_atoms() ?? 0;
  const bondPairs = candMolblock ? parseV2000Bonds(candMolblock) : [];
  const candAdj = buildAdjacency(bondPairs, candNumAtoms);
  const candAtoms = bondPairs.length > 0 && seedCandAtoms.size > 0 ?
    connectMatchedAtomsViaShortestPaths(seedCandAtoms, candAdj) :
    seedCandAtoms;

  // Bond indices for the highlight: every candidate bond whose endpoints
  // are both in the (path-completed) candAtoms. V2000-block 0-based
  // indices, lined up with what RDKit's `substruct_match.bonds` returns.
  const candBonds: number[] = [];
  if (candAtoms.size > 0) {
    for (let i = 0; i < bondPairs.length; i++) {
      const [a, b] = bondPairs[i];
      if (candAtoms.has(a) && candAtoms.has(b)) candBonds.push(i);
    }
  }

  return {refAtoms, candAtoms, candBonds, refRingMatches};
}
