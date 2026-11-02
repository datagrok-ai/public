import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

/* Engine-free chemistry checks on molblocks for the atomic-level tests (Bio on Crux). The tests no longer ask RDKit
 * (or the engine that built the molecule) whether a molecule is right: they read the molblock's own atom/bond graph
 * here and check the structural features directly. Each pattern function counts what the SMARTS in its comment
 * matches (distinct atom sets, as RDKit's uniquified matches), on heavy atoms; aromaticity is not perceived, so
 * aromatic SMARTS atoms are read as their ring and substituent pattern. Checked against RDKit's SMARTS counts on
 * nucleotides, cholesterol, biotin, GalNAc and LNA (V3000, V2000, explicit hydrogens). */

/** Canonical SMILES through Chem (the engine in effect), to compare a molecule with an expected SMILES: both sides
 * go through the same call. */
export async function canonicalSmiles(molecule: string): Promise<string> {
  return await grok.functions.call('Chem:convertMolNotation', {
    molecule, sourceNotation: DG.chem.Notation.Unknown, targetNotation: DG.chem.Notation.Smiles,
  });
}

export interface MolblockBond { a: number, b: number, order: number }

export interface MolblockGraph {
  /** Element symbol per atom (as written: 'C', 'O', 'R#', ...) */
  elements: string[];
  x: number[];
  y: number[];
  bonds: MolblockBond[];
  /** Heavy-atom neighbours per atom, with the bond order */
  nbrs: { atom: number, order: number }[][];
  /** SSSR rings over heavy atoms */
  rings: number[][];
  /** Number of SSSR rings per atom */
  ringCount: number[];
  /** Simple 5- and 6-membered cycles (atom lists in ring order), unique by atom set */
  smallCycles: number[][];
}

function joinV3KContinuations(lines: string[]): string[] {
  const res: string[] = [];
  let pending = '';
  for (const raw of lines) {
    const line = raw.replace(/\r$/, '');
    if (pending) {
      pending += line.replace(/^M {2}V30 /, '');
      if (!pending.endsWith('-')) {res.push(pending); pending = '';} else pending = pending.slice(0, -1);
    } else if (line.startsWith('M  V30 ') && line.endsWith('-'))
      pending = line.slice(0, -1);
    else
      res.push(line);
  }
  if (pending) res.push(pending);
  return res;
}

function parseAtomsAndBonds(molfile: string):
  { elements: string[], x: number[], y: number[], bonds: MolblockBond[] } {
  const elements: string[] = [];
  const x: number[] = [];
  const y: number[] = [];
  const bonds: MolblockBond[] = [];
  const lines = molfile.split('\n');
  if (molfile.includes('V3000')) {
    let block: 'none' | 'atom' | 'bond' = 'none';
    for (const line of joinV3KContinuations(lines)) {
      if (line.startsWith('M  V30 BEGIN ATOM')) block = 'atom';
      else if (line.startsWith('M  V30 BEGIN BOND')) block = 'bond';
      else if (line.startsWith('M  V30 END')) block = 'none';
      else if (block !== 'none' && line.startsWith('M  V30 ')) {
        const t = line.substring(7).trim().split(/\s+/);
        if (block === 'atom') {
          elements.push(t[1]);
          x.push(parseFloat(t[2]));
          y.push(parseFloat(t[3]));
        } else
          bonds.push({order: parseInt(t[1]), a: parseInt(t[2]) - 1, b: parseInt(t[3]) - 1});
      }
    }
  } else if (lines.length > 3) {
    const atomCount = parseInt(lines[3].substring(0, 3));
    const bondCount = parseInt(lines[3].substring(3, 6));
    for (let i = 0; i < atomCount; ++i) {
      const l = lines[4 + i];
      x.push(parseFloat(l.substring(0, 10)));
      y.push(parseFloat(l.substring(10, 20)));
      elements.push(l.substring(31, 34).trim());
    }
    for (let i = 0; i < bondCount; ++i) {
      const l = lines[4 + atomCount + i];
      bonds.push({a: parseInt(l.substring(0, 3)) - 1, b: parseInt(l.substring(3, 6)) - 1,
        order: parseInt(l.substring(6, 9))});
    }
  }
  return {elements, x, y, bonds};
}

const isHydrogen = (el: string): boolean => el === 'H' || el === 'D' || el === 'T';

/** Edge bitset helpers for the cycle basis */
function xorInto(dst: Uint32Array, src: Uint32Array): void {
  for (let i = 0; i < dst.length; ++i) dst[i] ^= src[i];
}

/** SSSR over heavy atoms: the shortest cycle through each edge as candidates, kept smallest first while independent
 * over GF(2), until the cyclomatic number is reached (enough for the fused and bridged rings of these tests). */
function findSssr(atomCount: number, edges: [number, number][], nbrs: { atom: number, order: number }[][]): number[][] {
  const edgeIdx = new Map<string, number>();
  edges.forEach(([a, b], i) => edgeIdx.set(a < b ? `${a},${b}` : `${b},${a}`, i));
  // cyclomatic number = E - V + components (heavy atoms with any bond, plus isolated ones)
  const parent = Array.from({length: atomCount}, (_, i) => i);
  const find = (v: number): number => {while (parent[v] !== v) v = parent[v] = parent[parent[v]]; return v;};
  for (const [a, b] of edges) parent[find(a)] = find(b);
  let components = 0;
  for (let i = 0; i < atomCount; ++i) if (find(i) === i) ++components;
  const cyclomatic = edges.length - atomCount + components;
  if (cyclomatic <= 0) return [];

  const candidates: number[][] = [];
  for (let e = 0; e < edges.length; ++e) {
    const [src, dst] = edges[e];
    // BFS from src to dst without the edge itself
    const prev = new Int32Array(atomCount).fill(-1);
    prev[src] = src;
    const queue = [src];
    for (let qi = 0; qi < queue.length && prev[dst] < 0; ++qi) {
      const v = queue[qi];
      for (const {atom: w} of nbrs[v]) {
        if (prev[w] >= 0 || (v === src && w === dst)) continue;
        prev[w] = v;
        queue.push(w);
      }
    }
    if (prev[dst] < 0) continue;
    const cycle: number[] = [];
    for (let v = dst; v !== src; v = prev[v]) cycle.push(v);
    cycle.push(src);
    candidates.push(cycle);
  }
  candidates.sort((p, q) => p.length - q.length);

  const words = Math.ceil(edges.length / 32);
  const basis: { vec: Uint32Array, pivot: number }[] = [];
  const rings: number[][] = [];
  const seen = new Set<string>();
  for (const cycle of candidates) {
    const key = [...cycle].sort((p, q) => p - q).join(',');
    if (seen.has(key)) continue;
    seen.add(key);
    const vec = new Uint32Array(words);
    for (let i = 0; i < cycle.length; ++i) {
      const a = cycle[i];
      const b = cycle[(i + 1) % cycle.length];
      const ei = edgeIdx.get(a < b ? `${a},${b}` : `${b},${a}`)!;
      vec[ei >>> 5] ^= 1 << (ei & 31);
    }
    for (const {vec: bv, pivot} of basis)
      if (vec[pivot >>> 5] & (1 << (pivot & 31))) xorInto(vec, bv);
    let pivot = -1;
    for (let w = 0; w < words && pivot < 0; ++w)
      if (vec[w]) pivot = w * 32 + (31 - Math.clz32(vec[w] & -vec[w])); // the lowest set bit
    if (pivot < 0) continue;
    // keep the basis reduced on the new pivot
    for (const bvItem of basis)
      if (bvItem.vec[pivot >>> 5] & (1 << (pivot & 31))) xorInto(bvItem.vec, vec);
    basis.push({vec, pivot});
    rings.push(cycle);
    if (rings.length === cyclomatic) break;
  }
  return rings;
}

/** Simple cycles of 5 and 6 heavy atoms, unique by atom set. */
function findSmallCycles(atomCount: number, nbrs: { atom: number, order: number }[][]): number[][] {
  const res: number[][] = [];
  const seen = new Set<string>();
  const path: number[] = [];
  const onPath = new Uint8Array(atomCount);
  const dfs = (start: number, v: number): void => {
    for (const {atom: w} of nbrs[v]) {
      if (w === start && path.length >= 5) {
        const key = [...path].sort((p, q) => p - q).join(',');
        if (!seen.has(key)) {seen.add(key); res.push([...path]);}
        continue;
      }
      if (onPath[w] || w < start || path.length >= 6) continue;
      onPath[w] = 1; path.push(w);
      dfs(start, w);
      onPath[w] = 0; path.pop();
    }
  };
  for (let s = 0; s < atomCount; ++s) {
    if (nbrs[s].length < 2) continue;
    onPath[s] = 1; path.push(s);
    dfs(s, s);
    onPath[s] = 0; path.pop();
  }
  return res;
}

export function parseMolblockGraph(molfile: string): MolblockGraph {
  const {elements, x, y, bonds} = parseAtomsAndBonds(molfile);
  const n = elements.length;
  const nbrs: { atom: number, order: number }[][] = Array.from({length: n}, () => []);
  const edges: [number, number][] = [];
  for (const bond of bonds) {
    if (bond.a < 0 || bond.b < 0 || bond.a >= n || bond.b >= n) continue;
    if (isHydrogen(elements[bond.a]) || isHydrogen(elements[bond.b])) continue;
    nbrs[bond.a].push({atom: bond.b, order: bond.order});
    nbrs[bond.b].push({atom: bond.a, order: bond.order});
    edges.push([bond.a, bond.b]);
  }
  const heavy = Array.from({length: n}, (_, i) => !isHydrogen(elements[i]));
  const rings = findSssr(n, edges, nbrs).filter((r) => r.every((i) => heavy[i]));
  const ringCount = new Array<number>(n).fill(0);
  for (const r of rings) for (const i of r) ++ringCount[i];
  return {elements, x, y, bonds, nbrs, rings, ringCount, smallCycles: findSmallCycles(n, nbrs)};
}

// -- atom-level predicates --

const deg = (g: MolblockGraph, i: number): number => g.nbrs[i].length;
const isEl = (g: MolblockGraph, i: number, el: string): boolean => g.elements[i] === el;
/** Only single (or aromatic) bonds: an sp3-like atom for CX4 / NX3 */
const allSingle = (g: MolblockGraph, i: number): boolean => g.nbrs[i].every((nb) => nb.order === 1 || nb.order === 4);
/** A methyl: carbon with one heavy neighbour over a single bond */
const isMethyl = (g: MolblockGraph, i: number): boolean =>
  isEl(g, i, 'C') && deg(g, i) === 1 && g.nbrs[i][0].order === 1;
const inRing = (g: MolblockGraph, i: number): boolean => g.ringCount[i] > 0;
/** Exocyclic neighbours of a ring atom */
const exoNbrs = (g: MolblockGraph, i: number, ring: number[]) => g.nbrs[i].filter((nb) => !ring.includes(nb.atom));

export function countElement(g: MolblockGraph, el: string): number {
  return g.elements.filter((e) => e === el).length;
}

export function heavyAtomCount(g: MolblockGraph): number {
  return g.elements.filter((e) => !isHydrogen(e)).length;
}

export function heavyBondCount(g: MolblockGraph): number {
  return g.bonds.filter((b) => !isHydrogen(g.elements[b.a]) && !isHydrogen(g.elements[b.b])).length;
}

/** Heavy-atom and heavy-bond counts (RDKit's get_num_atoms / get_num_bonds of a molblock read without Hs). */
export function heavyAtomBondCounts(molfile: string): { atoms: number, bonds: number } {
  const g = parseMolblockGraph(molfile);
  return {atoms: heavyAtomCount(g), bonds: heavyBondCount(g)};
}

/** Connected components over all atoms and bonds. */
export function fragmentCount(g: MolblockGraph): number {
  const n = g.elements.length;
  const parent = Array.from({length: n}, (_, i) => i);
  const find = (v: number): number => {while (parent[v] !== v) v = parent[v] = parent[parent[v]]; return v;};
  for (const b of g.bonds) if (b.a >= 0 && b.b >= 0 && b.a < n && b.b < n) parent[find(b.a)] = find(b.b);
  let res = 0;
  for (let i = 0; i < n; ++i) if (find(i) === i) ++res;
  return res;
}

// -- patterns --

/** Cycles of the given size with exactly these element counts, all ring bonds single or aromatic. */
function cyclesOf(g: MolblockGraph, size: number, counts: { [el: string]: number }): number[][] {
  return g.smallCycles.filter((c) => {
    if (c.length !== size) return false;
    for (let i = 0; i < c.length; ++i) {
      const a = c[i];
      const b = c[(i + 1) % c.length];
      const order = g.nbrs[a].find((nb) => nb.atom === b)?.order;
      if (order !== 1 && order !== 4) return false;
    }
    const got: { [el: string]: number } = {};
    for (const i of c) got[g.elements[i]] = (got[g.elements[i]] ?? 0) + 1;
    return Object.keys(counts).every((el) => got[el] === counts[el]) &&
      Object.keys(got).every((el) => counts[el] !== undefined);
  });
}

/** `[#6;R]1[#6;R][#6;R][#6;R][O;R]1` */
export const furanoseRings = (g: MolblockGraph): number[][] => cyclesOf(g, 5, {C: 4, O: 1});
/** `[#6]1[#6][#6][#6][#6][O]1` */
export const pyranoseRings = (g: MolblockGraph): number[][] => cyclesOf(g, 6, {C: 5, O: 1});
/** `[#6]1[#6][#6][#6][#6]1` */
export const cyclopentaneRings = (g: MolblockGraph): number[][] => cyclesOf(g, 5, {C: 5});

/** Atoms of 5- and 6-membered rings that contain nitrogen: the nucleobases (aromatic in RDKit, `[a]`). */
export function nucleobaseRingAtoms(g: MolblockGraph): Set<number> {
  const res = new Set<number>();
  for (const c of g.smallCycles)
    if (c.some((i) => isEl(g, i, 'N'))) for (const i of c) res.add(i);
  return res;
}

/** Phosphates `[#6][OX2][PX4](=[OX1])(X)[OX2][#6]`: P with a P=O, two O-C arms and a fourth substituent that
 * passes {@link isX}; one match per such P. */
function phosphateDiesters(g: MolblockGraph, isX: (atom: number, order: number) => boolean): number {
  let res = 0;
  for (let p = 0; p < g.elements.length; ++p) {
    if (!isEl(g, p, 'P') || deg(g, p) !== 4) continue;
    const nb = g.nbrs[p];
    const oxo = nb.filter((x) => x.order === 2 && isEl(g, x.atom, 'O') && deg(g, x.atom) === 1);
    if (oxo.length === 0) continue;
    const arms = nb.filter((x) => x.order === 1 && isEl(g, x.atom, 'O') && deg(g, x.atom) === 2 &&
      g.nbrs[x.atom].some((y) => y.atom !== p && isEl(g, y.atom, 'C')));
    if (arms.length < 2) continue;
    const rest = nb.filter((x) => x !== oxo[0] && !arms.slice(0, 2).includes(x));
    if (rest.length === 1 && isX(rest[0].atom, rest[0].order)) ++res;
  }
  return res;
}

/** `[#6][OX2][PX4](=[OX1])([OX2,SX2,OX1H,SX1H,OX1-,SX1-])[OX2][#6]` */
export function phosphodiesters(g: MolblockGraph): number {
  return phosphateDiesters(g, (a, order) => order === 1 && (isEl(g, a, 'O') || isEl(g, a, 'S')) && deg(g, a) <= 2);
}

/** `[#6][OX2][PX4](=[OX1])([SX2,SX1H,SX1-])[OX2][#6]` */
export function phosphorothioateDiesters(g: MolblockGraph): number {
  return phosphateDiesters(g, (a, order) => order === 1 && isEl(g, a, 'S') && deg(g, a) <= 2);
}

/** `[#6][OX2][PX4](=[OX1])([CH3])[OX2][#6]` */
export function methylphosphonateDiesters(g: MolblockGraph): number {
  return phosphateDiesters(g, (a, order) => order === 1 && isMethyl(g, a));
}

/** `[CX4][PX4]`: bonds between an sp3 carbon and a four-connected phosphorus */
export function directCPBonds(g: MolblockGraph): number {
  return g.bonds.filter(({a, b}) => {
    const [c, p] = isEl(g, a, 'C') ? [a, b] : [b, a];
    return isEl(g, c, 'C') && isEl(g, p, 'P') && allSingle(g, c) && deg(g, p) === 4;
  }).length;
}

/** `[PX4][OX2][PX4]`: oxygens bridging two phosphorus atoms */
export function popBridges(g: MolblockGraph): number {
  let res = 0;
  for (let o = 0; o < g.elements.length; ++o) {
    if (isEl(g, o, 'O') && deg(g, o) === 2 && allSingle(g, o) &&
      g.nbrs[o].every((nb) => isEl(g, nb.atom, 'P') && deg(g, nb.atom) === 4)) ++res;
  }
  return res;
}

/** `[PX4][SX2,SX1H,SX1-]`: single P-S bonds to a sulfur with at most two heavy neighbours */
export function phosphorusSulfurSingle(g: MolblockGraph): number {
  return g.bonds.filter(({a, b, order}) => {
    const [p, s] = isEl(g, a, 'P') ? [a, b] : [b, a];
    return order === 1 && isEl(g, p, 'P') && isEl(g, s, 'S') && deg(g, s) <= 2;
  }).length;
}

/** `[PX4]=S` */
export function phosphorusSulfurDouble(g: MolblockGraph): number {
  return g.bonds.filter(({a, b, order}) => order === 2 &&
    ((isEl(g, a, 'P') && isEl(g, b, 'S')) || (isEl(g, a, 'S') && isEl(g, b, 'P')))).length;
}

/** `[CX4][OX2][PX4](=[OX1])[SX2,SX1H,SX1-]`: C-O-P bridges on phosphorothioates, one per bridging oxygen */
export function thiophosphateCOPBridges(g: MolblockGraph): number {
  let res = 0;
  for (let o = 0; o < g.elements.length; ++o) {
    if (!isEl(g, o, 'O') || deg(g, o) !== 2 || !allSingle(g, o)) continue;
    const c = g.nbrs[o].find((nb) => isEl(g, nb.atom, 'C') && allSingle(g, nb.atom));
    const p = g.nbrs[o].find((nb) => isEl(g, nb.atom, 'P') && deg(g, nb.atom) === 4);
    if (!c || !p) continue;
    const pn = g.nbrs[p.atom];
    const hasOxo = pn.some((x) => x.order === 2 && isEl(g, x.atom, 'O') && deg(g, x.atom) === 1);
    const hasS = pn.some((x) => x.order === 1 && isEl(g, x.atom, 'S') && deg(g, x.atom) <= 2);
    if (hasOxo && hasS) ++res;
  }
  return res;
}

/** `[F][CX4;R]` */
export function fluorineOnRingCarbon(g: MolblockGraph): number {
  return g.bonds.filter(({a, b}) => {
    const [f, c] = isEl(g, a, 'F') ? [a, b] : [b, a];
    return isEl(g, f, 'F') && isEl(g, c, 'C') && allSingle(g, c) && inRing(g, c);
  }).length;
}

/** `[NX3]C(=O)[CH3]` */
export function nAcetylGroups(g: MolblockGraph): number {
  let res = 0;
  for (let c = 0; c < g.elements.length; ++c) {
    if (!isEl(g, c, 'C')) continue;
    const nb = g.nbrs[c];
    const oxo = nb.some((x) => x.order === 2 && isEl(g, x.atom, 'O'));
    const methyl = nb.some((x) => x.order === 1 && isMethyl(g, x.atom));
    const amideN = nb.filter((x) => x.order === 1 && isEl(g, x.atom, 'N') && allSingle(g, x.atom));
    if (oxo && methyl) res += amideN.length;
  }
  return res;
}

/** `[CH3][OX2][#6;R]` */
export function methoxyOnRingCarbon(g: MolblockGraph): number {
  let res = 0;
  for (let o = 0; o < g.elements.length; ++o) {
    if (!isEl(g, o, 'O') || deg(g, o) !== 2 || !allSingle(g, o)) continue;
    const [p, q] = g.nbrs[o].map((nb) => nb.atom);
    if ((isMethyl(g, p) && isEl(g, q, 'C') && inRing(g, q)) || (isMethyl(g, q) && isEl(g, p, 'C') && inRing(g, p)))
      ++res;
  }
  return res;
}

/** `[#6;R2]`: carbons in exactly two SSSR rings */
export function carbonsInTwoRings(g: MolblockGraph): number {
  return g.elements.filter((el, i) => el === 'C' && g.ringCount[i] === 2).length;
}

/** Bonds between atoms of two elements (any order), e.g. S-S bridges */
export function bondsBetween(g: MolblockGraph, el1: string, el2: string): number {
  return g.bonds.filter(({a, b}) =>
    (isEl(g, a, el1) && isEl(g, b, el2)) || (isEl(g, a, el2) && isEl(g, b, el1))).length;
}

/** Double bonds outside rings (C=O, C=C, C=N of side chains and backbones), independent of aromaticity */
export function chainDoubleBonds(g: MolblockGraph): number {
  return g.bonds.filter(({a, b, order}) => order === 2 &&
    !g.rings.some((r) => r.includes(a) && r.includes(b))).length;
}

/** `[#16;R]` */
export function ringSulfurs(g: MolblockGraph): number {
  return g.elements.filter((el, i) => el === 'S' && inRing(g, i)).length;
}

/** Six-membered rings with two nitrogens and four carbons (pyrimidines, the six-ring of purines) */
function diazineRings(g: MolblockGraph): number[][] {
  return g.smallCycles.filter((c) => c.length === 6 &&
    c.filter((i) => isEl(g, i, 'N')).length === 2 && c.filter((i) => isEl(g, i, 'C')).length === 4);
}

function ringSubstituents(g: MolblockGraph, ring: number[]) {
  const carbonsWith = (pred: (x: { atom: number, order: number }) => boolean): number[] =>
    ring.filter((i) => isEl(g, i, 'C') && exoNbrs(g, i, ring).some(pred));
  const oxo = carbonsWith((x) => x.order === 2 && isEl(g, x.atom, 'O'));
  const amino = carbonsWith((x) => x.order === 1 && isEl(g, x.atom, 'N'));
  const methyl = carbonsWith((x) => isMethyl(g, x.atom));
  const nSubst = ring.filter((i) => isEl(g, i, 'N') && exoNbrs(g, i, ring).length > 0);
  return {oxo, amino, methyl, nSubst};
}

/** `[CH3]c1cn([!#1])c(=O)nc1N`: 5-methylcytosine rings */
export function methylCytosineRings(g: MolblockGraph): number {
  return diazineRings(g).filter((r) => {
    const s = ringSubstituents(g, r);
    return s.oxo.length === 1 && s.amino.length === 1 && s.methyl.length === 1 && s.nSubst.length >= 1;
  }).length;
}

/** `[CH3][#6]1=[#6][#7]([!#1])[#6](=O)[#7][#6]1=O` / `Cc1cn([!#1])c(=O)[nH]c1=O`: thymine rings */
export function thymineRings(g: MolblockGraph): number {
  return diazineRings(g).filter((r) => {
    const s = ringSubstituents(g, r);
    return s.oxo.length === 2 && s.methyl.length === 1 && s.nSubst.length >= 1;
  }).length;
}

/** `n1cnc2c1ncnc2N`: adenine (a purine whose six-ring carbon next to the fusion bears an amino group) */
export function adenineRings(g: MolblockGraph): number {
  let res = 0;
  for (const six of diazineRings(g)) {
    for (const five of g.smallCycles) {
      if (five.length !== 5) continue;
      const shared = five.filter((i) => six.includes(i));
      if (shared.length !== 2 || five.filter((i) => isEl(g, i, 'N')).length !== 2) continue;
      // six-ring carbons bonded to a fusion atom, not fusion atoms themselves
      const nextToFusion = six.filter((i) => !shared.includes(i) && isEl(g, i, 'C') &&
        g.nbrs[i].some((nb) => shared.includes(nb.atom)));
      if (nextToFusion.some((i) => exoNbrs(g, i, six).some((x) => x.order === 1 && isEl(g, x.atom, 'N'))))
        ++res;
    }
  }
  return res;
}

/** `[#7;R]1[#6;R](=[OX1])[#7;R][#6;R][#6;R]1`: cyclic ureas (biotin's ureido ring) */
export function cyclicUreaRings(g: MolblockGraph): number {
  return g.smallCycles.filter((c) => {
    if (c.length !== 5 || c.filter((i) => isEl(g, i, 'N')).length !== 2) return false;
    return c.some((i) => isEl(g, i, 'C') && g.nbrs[i].filter((nb) => c.includes(nb.atom))
      .every((nb) => isEl(g, nb.atom, 'N')) &&
      exoNbrs(g, i, c).some((x) => x.order === 2 && isEl(g, x.atom, 'O') && deg(g, x.atom) === 1));
  }).length;
}
