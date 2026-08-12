import * as DG from 'datagrok-api/dg';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {tanimotoSimilarity} from '@datagrok-libraries/ml/src/distance-metrics-methods';

import {Fingerprint, rdKitFingerprintToBitArray} from '../../utils/chem-common';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {bitbirchWorker} from '../bit-birch/bitbirch-clustering';
import {getMmpFrags} from '../molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {MmpFragments} from '../molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {CoreCluster, MatchedSeries, SeriesMember} from './sar-matrix-types';

/**
 * Step 3 — build matched series directly from the MMP fragmentation.
 *
 * Groups molecules that share the same core fragment; each such group with at
 * least two members is a matched series, and the varying substituent per member
 * is the value fragment of that cut. No `MMPA.init` is needed — this reads the
 * raw `getMmpFrags` output and applies the fragment-size cutoff itself (mirroring
 * `getMmpRulesCPU`).
 *
 * @param frags Output of `getMmpFrags` for the molecule set.
 * @param fragmentCutoff Maximum substituent-to-core size ratio to keep a cut.
 * @returns One `MatchedSeries` per core with at least two members.
 */
export function buildMatchedSeries(frags: MmpFragments, fragmentCutoff: number): MatchedSeries[] {
  const {fragCodes, idToName, sizes} = frags;
  const membersByCore = new Map<number, Map<number, number>>(); // coreId -> (molIdx -> substId)

  for (let molIdx = 0; molIdx < fragCodes.length; molIdx++) {
    for (const [coreId, substId] of fragCodes[molIdx]) {
      if (idToName[coreId] === '' || sizes[coreId] === 0)
        continue;
      if (sizes[substId] / sizes[coreId] > fragmentCutoff)
        continue;
      let mols = membersByCore.get(coreId);
      if (!mols)
        membersByCore.set(coreId, mols = new Map<number, number>());
      // First cut that yields this core wins for a given molecule.
      if (!mols.has(molIdx))
        mols.set(molIdx, substId);
    }
  }

  const series: MatchedSeries[] = [];
  for (const [coreFragId, mols] of membersByCore) {
    if (mols.size < 2)
      continue;
    const members: SeriesMember[] = [];
    for (const [molIdx, substId] of mols)
      members.push({molIdx, substSmiles: idToName[substId]});
    series.push({coreSmiles: idToName[coreFragId], members});
  }
  return series;
}

/** Isotope standing in for a series' own attachment point while its core is cut a second time. Without
 *  it both the old and the newly opened point come back as `[*:1]`, so cores varying at two unrelated
 *  sites collapse onto one group key and land in a single block-diagonal matrix. */
const SITE_MARKER = '[3*]';

/** A second-cut value holding no atoms at all — the cut severed the marker dummy itself. Its key is
 *  just the core with its own site reopened, which no other core can produce, so the cut can never
 *  form a matrix and only costs an index entry. */
const ATOMLESS_VALUE = /^(?:\[3\*\]|\[\*:\d+\])+$/;

/** True when every member of `a` is also in `b`. */
function isSubset(a: Set<number>, b: Set<number>): boolean {
  if (a.size > b.size)
    return false;
  for (const x of a) {
    if (!b.has(x))
      return false;
  }
  return true;
}

/**
 * Step 4 — group series into the rows of one matrix by cutting their cores a second time.
 *
 * Two cores that differ at exactly one point leave the same remainder once that point is cut away, so
 * grouping on the remainder ("key of key") is exact: no fingerprint, no similarity threshold and
 * nothing to tune. Membership is derived rather than chosen — a core is a row of every group its own
 * cuts produce — so the result cannot depend on evaluation order and two runs are identical.
 *
 * The series' own attachment point is marked before the cut so it stays distinguishable from the one
 * the cut opens; without that, both come back as `[*:1]` and two unrelated sites collapse onto one key.
 * Whether the marker survives into the key decides only how a matrix is *labelled*, not whether it is
 * built: a key that kept it groups cores varying away from their substituent, so a column compares one
 * position across rows; a key that lost it groups cores whose substituent rides on the very part that
 * differs, so a column compares the same substituent on different scaffolds. Both are matrices under
 * the method, which keys columns by substituent identity and not by positional equivalence — the
 * second kind is what puts a ring hop in a row, and dropping it would forbid that outright.
 *
 * A core can be cut at many bonds, so it joins one group per cuttable site — an overlapping cover
 * rather than a partition, and a compound appearing in several matrices is the intended behaviour,
 * each matrix asking a different question of it. The count stays linear in the number of series, since
 * a core has only a handful of cuttable bonds. Groups whose series are all contained in a larger group
 * are dropped, since they only restate it over fewer rows; groups that merely overlap are kept, being
 * genuinely different comparisons.
 *
 * No size ratio is applied here. The larger fragment is already taken as the key, which bounds the
 * removed piece at half the core; the level-1 cutoff answers "substituent or half a molecule?", a
 * question that has no meaning once the thing being cut is itself a core. The only cut worth
 * discarding is the one that carries no chemical change at all — see `ATOMLESS_VALUE`.
 *
 * @param series Matched series from `buildMatchedSeries`.
 * @returns One `CoreCluster` per site, plus a singleton per series no site grouped.
 */
export async function groupSeriesBySite(series: MatchedSeries[]): Promise<CoreCluster[]> {
  if (series.length < 2)
    return series.map((s, i) => ({id: `c${i}`, series: [s], siteKey: '', level: 2}));

  const [frags] = await getMmpFrags(series.map((s) => s.coreSmiles.replace(/\[\*:\d+\]/g, SITE_MARKER)));
  const {fragCodes, idToName, sizes} = frags;

  const bySite = new Map<string, Set<number>>();
  for (let i = 0; i < series.length; i++) {
    for (const [keyId, valueId] of fragCodes[i] ?? []) {
      const key = idToName[keyId];
      if (!key || sizes[keyId] === 0 || ATOMLESS_VALUE.test(idToName[valueId] ?? ''))
        continue;
      let group = bySite.get(key);
      if (!group)
        bySite.set(key, group = new Set<number>());
      group.add(i);
    }
  }

  // Largest first, so a group survives only when it covers series the ones already kept do not. Two
  // groups over the identical series set are subsets of each other, so exactly one of them survives.
  const groups = [...bySite.entries()].filter(([, g]) => g.size > 1).sort((a, b) => b[1].size - a[1].size);
  const kept: {key: string, members: Set<number>}[] = [];
  for (const [key, members] of groups) {
    if (!kept.some((k) => isSubset(members, k.members)))
      kept.push({key, members});
  }

  const clusters: CoreCluster[] = [];
  const placed = new Set<number>();
  for (const {key, members} of kept) {
    // Sorted so a cluster's row order follows the series order and does not depend on which cut of
    // which core happened to reach the group first.
    const ordered = [...members].sort((a, b) => a - b);
    for (const i of ordered)
      placed.add(i);
    clusters.push({id: `c${clusters.length}`, series: ordered.map((i) => series[i]),
      siteKey: key, level: 2});
  }
  // A lone series has no second row to be aligned against, so nothing about it is nominal.
  for (let i = 0; i < series.length; i++) {
    if (!placed.has(i))
      clusters.push({id: `c${clusters.length}`, series: [series[i]], siteKey: '', level: 2});
  }
  return clusters;
}

/** Any site marked by an earlier level. Each level marks the points it inherits with its own isotope
 *  before cutting, so a key carries one dummy per level that opened a site still present in it. */
const ANY_MARKER = /\[\d+\*\]/;

/**
 * One further round of fragmentation over a set of keys: mark whatever points each key already
 * carries, cut, and group keys leaving the same remainder.
 *
 * A key must retain at least one marked site to be grouped. Without that a family can form on a
 * fragment that merely occurs in several keys by coincidence — a bare 4-substituted phenyl gathers
 * chemotypes sharing nothing else — and every level makes the surviving fragment smaller and so more
 * prone to it.
 *
 * @param keys One key per item; an empty entry never joins a group.
 * @param marker Isotope dummy standing in for this level's inherited attachment points.
 * @returns The group key per item, empty where nothing grouped it. An item reachable from several
 *   groups takes the largest, so the result partitions its input.
 */
async function groupKeysOnce(keys: string[], marker: string): Promise<string[]> {
  const grouped = new Array<string>(keys.length).fill('');
  const indexed = keys.map((k, i) => ({k, i})).filter((e) => e.k !== '');
  if (indexed.length < 2)
    return grouped;

  const [frags] = await getMmpFrags(indexed.map((e) => e.k.replace(/\[\*:\d+\]/g, marker)));
  const {fragCodes, idToName, sizes} = frags;

  const byKey = new Map<string, Set<number>>();
  for (let n = 0; n < indexed.length; n++) {
    for (const [keyId, valueId] of fragCodes[n] ?? []) {
      const key = idToName[keyId];
      if (!key || sizes[keyId] === 0 || ATOMLESS_VALUE.test(idToName[valueId] ?? '') || !ANY_MARKER.test(key))
        continue;
      let group = byKey.get(key);
      if (!group)
        byKey.set(key, group = new Set<number>());
      group.add(indexed[n].i);
    }
  }

  // Largest group first; the key tie-break keeps equally-sized groups from swapping run to run.
  for (const [key, members] of [...byKey.entries()].filter(([, g]) => g.size > 1)
    .sort((a, b) => b[1].size - a[1].size || (a[0] < b[0] ? -1 : 1))) {
    for (const i of members) {
      if (grouped[i] === '')
        grouped[i] = key;
    }
  }
  return grouped;
}

/**
 * Step 4b — build the coarser levels above the matrices, as many as asked for.
 *
 * Two matrices whose cores become the same core one level deeper belong together: a benzanilide
 * matrix and an isoxazole-carboxamide one differ by more than a single site, so no level-2 key can
 * hold both, yet strip the ring from each key and the same remainder falls out. Those siblings are
 * folded into one coarser group holding every series they contain between them, which assembles into
 * a matrix of its own — the parent covers its children's compounds over fewer, broader rows.
 *
 * Each further level removes more of the molecule, so groups get broader as the level rises. This is
 * the dial between many small matrices and few large ones, and it is a chemical statement — how much
 * of the molecule is held constant — rather than a threshold to tune.
 *
 * A group that gathers nothing at a level simply stops there, so branches end at different depths.
 *
 * @param clusters The level-2 groups, from `groupSeriesBySite`.
 * @param extraLevels How many coarser levels to build above them; 0 adds none.
 * @returns Every group, the input included, each carrying its level and the id of its parent.
 */
export async function buildCoarserLevels(clusters: CoreCluster[], extraLevels: number): Promise<CoreCluster[]> {
  const all = [...clusters];
  let current = clusters;
  for (let step = 0; step < extraLevels && current.length > 1; step++) {
    // The marker index has to clear SITE_MARKER's isotope and every level's before it, or a fresh
    // mark would be indistinguishable from an inherited one and the two sites would collapse.
    const grouped = await groupKeysOnce(current.map((c) => c.siteKey), `[${4 + step}*]`);
    const byKey = new Map<string, CoreCluster[]>();
    grouped.forEach((key, i) => {
      if (key === '')
        return;
      byKey.set(key, [...(byKey.get(key) ?? []), current[i]]);
    });
    if (byKey.size === 0)
      break;

    const parents: CoreCluster[] = [];
    for (const [key, members] of byKey) {
      if (members.length < 2)
        continue;
      const id = `L${3 + step}-${parents.length}`;
      // Series are unioned by core, so a series reached through two children is not made a row twice.
      const series = new Map<string, MatchedSeries>();
      for (const member of members) {
        for (const s of member.series)
          series.set(s.coreSmiles, s);
        member.parentId = id;
      }
      parents.push({id, series: [...series.values()], siteKey: key, level: 3 + step});
    }
    if (parents.length === 0)
      break;
    all.push(...parents);
    current = parents;
  }
  return all;
}


/** Past this many cores the all-pairs comparison stops being worth its cost. */
const MAX_BUTINA_CORES = 3000;

/**
 * Taylor-Butina sphere-exclusion clustering over core fingerprints, returning a cluster id per input
 * (-1 for a structure whose fingerprint could not be built).
 *
 * Each core's neighbours within `cutoff` Tanimoto are collected, the core with the most neighbours
 * claims them all as a cluster, and the process repeats over what is left. Nothing forces a lone core
 * into a cluster it barely resembles: it simply ends up alone, which matters because a core absorbed
 * on a weak match brings its whole series into a matrix it does not belong in.
 *
 * The cost of that rule is that membership is decided against cluster centres rather than against
 * whole clusters, so a core can sit above the cutoff from another core and still land elsewhere when
 * neither of them becomes a centre. Raising `threshold` tightens the centres and reduces it.
 *
 * Neighbour counts are taken once, before any assignment, so the result does not depend on input order.
 */
async function butinaCluster(structures: string[], cutoff: number,
  fingerprintType: Fingerprint): Promise<number[]> {
  const service = await getRdKitService();
  const raw = (await service.getFingerprints(fingerprintType, structures, false)).fps;
  const fps: (BitArray | null)[] = raw.map((fp) => fp ? rdKitFingerprintToBitArray(fp) : null);

  const n = fps.length;
  const neighbours: number[][] = Array.from({length: n}, () => []);
  for (let i = 0; i < n; i++) {
    if (!fps[i])
      continue;
    for (let j = i + 1; j < n; j++) {
      if (!fps[j] || tanimotoSimilarity(fps[i]!, fps[j]!) < cutoff)
        continue;
      neighbours[i].push(j);
      neighbours[j].push(i);
    }
  }

  const assigned = new Array<number>(n).fill(-1);
  // Densest first; the index tie-break keeps equally-connected cores from swapping roles run to run.
  const order = [...neighbours.keys()]
    .sort((a, b) => neighbours[b].length - neighbours[a].length || a - b);
  let cluster = 0;
  for (const centre of order) {
    if (assigned[centre] !== -1 || !fps[centre])
      continue;
    assigned[centre] = cluster;
    for (const neighbour of neighbours[centre]) {
      if (assigned[neighbour] === -1)
        assigned[neighbour] = cluster;
    }
    cluster++;
  }
  return assigned;
}

/**
 * Step 4, by similarity instead of by site — cluster structurally related series into one matrix.
 *
 * Unlike `groupSeriesBySite` this reaches across a change inside the core, but it gives no guarantee
 * that the rows of a matrix vary at a common site, which is what the reference substituents and the
 * folded row identities in assembly exist to work around.
 *
 * Series are grouped by the fingerprint similarity of their own cores, so similarity (not identity)
 * still puts ring hops such as pyrazole and thiazole in one matrix, which is what SAR transfer
 * between rows depends on. The Murcko scaffold is deliberately not the key — it is still the right
 * anchor for the downstream decomposition, but as a similarity key it carries every ring in the
 * molecule, including a substituent's, and that swamps the shared core.
 *
 * Morgan circular fingerprints compare exact atom environments, so cores differing by a ring atom
 * separate instead of merging. MACCS groups such ring hops together, which on diverse sets collapses
 * hundreds of unrelated cores into one matrix. Note that a threshold does not carry across
 * fingerprint types — Morgan Tanimoto runs lower than MACCS for the same pair, so the same number is a
 * tighter cut here and produces smaller clusters.
 *
 * @param series Matched series from `buildMatchedSeries`.
 * @param threshold Similarity cutoff; higher keeps clusters tighter.
 * @param fingerprintType Fingerprint used for core similarity.
 * @returns One `CoreCluster` per group of related cores.
 */
export async function clusterRelatedCores(series: MatchedSeries[],
  threshold: number = 0.5, fingerprintType: Fingerprint = Fingerprint.Morgan): Promise<CoreCluster[]> {
  if (series.length === 0)
    return [];
  if (series.length === 1)
    return [{id: 'c0', series: [series[0]], siteKey: '', level: 2}];

  // A series is compared by its own core, which is exactly what its members share and is identical
  // for every one of them. A whole-molecule scaffold is the wrong key here: it carries every ring in
  // the molecule, so a ring-bearing substituent joins the comparison and dominates it — two series
  // over the same core score far apart when one varies an aryl and the other an alkyl. Reading the
  // core also removes the dependence on which member happened to be listed first.
  const keys = series.map((s) => s.coreSmiles);

  let assignment: number[];
  if (keys.length <= MAX_BUTINA_CORES)
    assignment = await butinaCluster(keys, threshold, fingerprintType);
  else {
    // Butina compares every pair, so past this size the quadratic pass costs more than the sharper
    // clusters are worth; BitBIRCH builds its tree in one incremental sweep instead.
    const coreCol = DG.Column.fromStrings('core', keys);
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    const clusterCol = await bitbirchWorker(coreCol, threshold, fingerprintType);
    assignment = series.map((_s, i) => clusterCol.isNone(i) ? -1 : clusterCol.get(i));
  }

  // Group series by cluster id; anything left unplaced becomes its own singleton.
  const seriesByCluster = new Map<number, MatchedSeries[]>();
  let singletonKey = -1;
  for (let i = 0; i < series.length; i++) {
    const cluster = assignment[i] < 0 ? singletonKey-- : assignment[i];
    let group = seriesByCluster.get(cluster);
    if (!group)
      seriesByCluster.set(cluster, group = []);
    group.push(series[i]);
  }

  let n = 0;
  const clusters: CoreCluster[] = [];
  // Every row here is decomposed against one shared anchor scaffold, which is what makes this path's
  // columns comparable in the first place, so its positions are aligned by construction.
  for (const group of seriesByCluster.values())
    clusters.push({id: `c${n++}`, series: group, siteKey: '', level: 2});
  return clusters;
}
