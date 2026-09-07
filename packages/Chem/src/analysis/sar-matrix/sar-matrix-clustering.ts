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
 * Step 3 — build matched series from MMP fragmentation: molecules sharing a core
 * form a series, with the fragment-size cutoff applied here (mirrors `getMmpRulesCPU`).
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
  // Canonical order, keyed on the chemistry rather than on fragment ids. Those ids are minted as the
  // workers discover fragments, so this map's order varies between runs even though its contents do
  // not. Ordering once here settles every downstream tie at the same time.
  for (const matched of series)
    matched.members.sort((a, b) => a.molIdx - b.molIdx);
  series.sort((a, b) => a.coreSmiles < b.coreSmiles ? -1 : a.coreSmiles > b.coreSmiles ? 1 : 0);
  return series;
}

/** Marks a series' own attachment point during a second cut so it stays distinct from the newly
 *  opened one; without it both return as `[*:1]` and unrelated sites collapse onto one group key. */
const SITE_MARKER = '[3*]';

/** A second-cut value with no atoms (the cut severed the marker dummy itself); it can never form a
 *  matrix, so it is skipped. */
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
 * Group the series by a value the user supplies per compound: everything carrying the same value
 * becomes one matrix, labelled with that value. A core whose compounds carry different values is split
 * across them, so a compound never lands in a matrix whose name contradicts its own column. Compounds
 * with no value sit out rather than forming a "blank" series the user never named.
 */
export function groupSeriesByColumn(series: MatchedSeries[], values: (string | null)[]): CoreCluster[] {
  const byValue = new Map<string, Map<string, MatchedSeries>>();
  for (const matched of series) {
    for (const member of matched.members) {
      const value = values[member.molIdx];
      if (value === null || value === '')
        continue;
      let cores = byValue.get(value);
      if (cores === undefined)
        byValue.set(value, cores = new Map<string, MatchedSeries>());
      const core = cores.get(matched.coreSmiles);
      if (core === undefined)
        cores.set(matched.coreSmiles, {coreSmiles: matched.coreSmiles, members: [member]});
      else
        core.members.push(member);
    }
  }
  return [...byValue.entries()].map(([value, cores], i) => ({
    id: `u${i}`, series: [...cores.values()], siteKey: '', level: 2, label: value,
  }));
}

/**
 * Series indices per site: cutting each core once more, the remainder is the site its substituents
 * hang off, and two cores leaving the same remainder vary at the same place. A core reachable by
 * several cuts appears under each of them, so the sets overlap.
 */
async function seriesBySite(series: MatchedSeries[]): Promise<Map<string, Set<number>>> {
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
  return bySite;
}

export async function groupSeriesBySite(series: MatchedSeries[]): Promise<CoreCluster[]> {
  if (series.length < 2)
    return series.map((s, i) => ({id: `c${i}`, series: [s], siteKey: '', level: 2}));

  const bySite = await seriesBySite(series);

  // Largest first, so a group survives only when it covers series the kept ones do not; which of two
  // same-sized groups is kept decides the whole set below it, hence the key tiebreak.
  const groups = [...bySite.entries()].filter(([, g]) => g.size > 1)
    .sort((a, b) => b[1].size - a[1].size || (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0));
  const kept: {key: string, members: Set<number>}[] = [];
  for (const [key, members] of groups) {
    if (!kept.some((k) => isSubset(members, k.members)))
      kept.push({key, members});
  }

  const clusters: CoreCluster[] = [];
  const placed = new Set<number>();
  for (const {key, members} of kept) {
    const ordered = [...members].sort((a, b) => a - b);
    for (const i of ordered)
      placed.add(i);
    clusters.push({id: `c${clusters.length}`, series: ordered.map((i) => series[i]),
      siteKey: key, level: 2});
  }
  for (let i = 0; i < series.length; i++) {
    if (!placed.has(i))
      clusters.push({id: `c${clusters.length}`, series: [series[i]], siteKey: '', level: 2});
  }
  return clusters;
}

/** Matches any earlier level's attachment-point isotope. A key must retain a marked site to be
 *  grouped, otherwise a family forms on a fragment shared only by coincidence. */
const ANY_MARKER = /\[\d+\*\]/;

/**
 * One further round of fragmentation over a set of keys: mark each key's existing points, cut, and
 * group keys leaving the same remainder. An item reachable from several groups takes the largest.
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

  // Largest group first; the key tie-break keeps equally-sized groups stable run to run.
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
 * Step 4b — build coarser levels above the matrices. Matrices whose keys leave the same remainder one
 * level deeper are folded into a parent group holding every series they contain between them. Each
 * level removes more of the molecule; a group that gathers nothing simply stops there.
 */
export async function buildCoarserLevels(clusters: CoreCluster[], extraLevels: number): Promise<CoreCluster[]> {
  const all = [...clusters];
  let current = clusters;
  for (let step = 0; step < extraLevels && current.length > 1; step++) {
    // Marker isotope must clear SITE_MARKER's and every prior level's, or a fresh mark would be
    // indistinguishable from an inherited one and the two sites would collapse.
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
      // Union by core, so a series reached through two children is not made a row twice.
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
 * Taylor-Butina sphere-exclusion clustering over core fingerprints (cluster id per input, -1 when the
 * fingerprint could not be built). The densest core claims its within-`cutoff` neighbours as a cluster,
 * repeating over the rest; a lone core stays alone rather than being absorbed on a weak match.
 * Neighbour counts are taken once, before assignment, so the result is order-independent.
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
  // Densest first; the index tie-break keeps equally-connected cores stable run to run.
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
 * Step 4, by similarity instead of by site — cluster related series into one matrix by the fingerprint
 * similarity of their own cores. This reaches across a change inside the core (ring hops such as
 * pyrazole/thiazole land together) but gives no guarantee that rows vary at a common site. The core,
 * not the Murcko scaffold, is the key: the scaffold carries a substituent's rings and swamps the core.
 * A threshold does not carry across fingerprint types — Morgan Tanimoto runs lower than MACCS.
 */
export async function clusterRelatedCores(series: MatchedSeries[],
  threshold: number = 0.5, fingerprintType: Fingerprint = Fingerprint.Morgan): Promise<CoreCluster[]> {
  if (series.length === 0)
    return [];
  if (series.length === 1)
    return [{id: 'c0', series: [series[0]], siteKey: '', level: 2}];

  const keys = series.map((s) => s.coreSmiles);

  const assignment = await similarityAssignment(keys, threshold, fingerprintType);

  // Group series by cluster id; anything unplaced becomes its own singleton.
  const seriesByCluster = new Map<number, number[]>();
  let singletonKey = -1;
  for (let i = 0; i < series.length; i++) {
    const cluster = assignment[i] < 0 ? singletonKey-- : assignment[i];
    let group = seriesByCluster.get(cluster);
    if (!group)
      seriesByCluster.set(cluster, group = []);
    group.push(i);
  }

  // Similarity alone lets a cluster hold cores whose substituents hang off different places, and the
  // matrix pools those into one axis — a column that means one position on some rows and another on
  // the rest, comparable across neither. Split so the rows of a cluster share a site.
  const bySite = await seriesBySite(series);
  const clusters: CoreCluster[] = [];
  const emit = (members: number[], siteKey: string): void => {
    clusters.push({id: `c${clusters.length}`, series: members.map((i) => series[i]), siteKey, level: 2});
  };
  for (const group of seriesByCluster.values()) {
    if (group.length < 2) {
      emit(group, '');
      continue;
    }
    let remaining = new Set(group);
    while (remaining.size > 1) {
      let bestKey = '';
      let best: number[] = [];
      for (const [key, members] of bySite) {
        const shared = [...remaining].filter((i) => members.has(i));
        // Key breaks a tie so equally large shares always split the same way.
        if (shared.length > best.length || (shared.length === best.length && shared.length > 1 &&
          key < bestKey))
          [best, bestKey] = [shared, key];
      }
      if (best.length < 2)
        break;
      emit(best.sort((a, b) => a - b), bestKey);
      remaining = new Set([...remaining].filter((i) => !best.includes(i)));
    }
    for (const i of [...remaining].sort((a, b) => a - b))
      emit([i], '');
  }
  return clusters;
}

/**
 * Cluster id per structure, `-1` for anything left unplaced. Butina compares all pairs, which is
 * quadratic, so past {@link MAX_BUTINA_CORES} BitBIRCH builds its tree in one incremental sweep.
 */
async function similarityAssignment(structures: string[], threshold: number,
  fingerprintType: Fingerprint): Promise<number[]> {
  if (structures.length <= MAX_BUTINA_CORES)
    return butinaCluster(structures, threshold, fingerprintType);
  const col = DG.Column.fromStrings('structure', structures);
  col.semType = DG.SEMTYPE.MOLECULE;
  const clusterCol = await bitbirchWorker(col, threshold, fingerprintType);
  return structures.map((_s, i) => clusterCol.isNone(i) ? -1 : clusterCol.get(i));
}

/**
 * Similarity groups over `structures`, as index groups ordered by the structure each starts with — so
 * which group is which does not follow the order the workers returned them. Singleton groups are
 * kept; the caller decides whether one structure on its own is worth a cluster.
 */
async function similarityGroups(structures: string[], threshold: number,
  fingerprintType: Fingerprint): Promise<number[][]> {
  const assignment = await similarityAssignment(structures, threshold, fingerprintType);
  const byCluster = new Map<number, number[]>();
  let singletonKey = -1;
  for (let i = 0; i < structures.length; i++) {
    const key = assignment[i] < 0 ? singletonKey-- : assignment[i];
    // Appended, not rebuilt: the pooled set runs to thousands of structures, and copying the group on
    // every insert makes filling one cluster quadratic in its own size.
    const group = byCluster.get(key);
    if (group === undefined)
      byCluster.set(key, [i]);
    else
      group.push(i);
  }
  return [...byCluster.values()]
    // Equal structures must compare 0: duplicates are routine here, and a comparator that never
    // reports equality leaves their order up to the sort implementation, which is the one thing this
    // ordering exists to prevent.
    .map((group) => [...group].sort((a, b) =>
      structures[a] < structures[b] ? -1 : structures[a] > structures[b] ? 1 : a - b))
    .sort((a, b) => structures[a[0]] < structures[b[0]] ? -1 :
      structures[a[0]] > structures[b[0]] ? 1 : a[0] - b[0]);
}

/**
 * Pool the series no core grouping could place, so an MCS has something to anchor.
 *
 * A series left on its own holds one core, and one core is one row: it can never become a matrix, so
 * its compounds simply do not appear. Grouping several by core similarity gives them rows to compare
 * against each other, and since they share no cut site — that is why they were left over — an MCS is
 * the only thing that can find their common core. Clusters the grouping placed are untouched.
 */
export async function poolUngroupedSeries(clusters: CoreCluster[], threshold: number,
  fingerprintType: Fingerprint = Fingerprint.Morgan): Promise<CoreCluster[]> {
  const isLone = (c: CoreCluster): boolean => !c.siteKey && c.series.length === 1;
  const lone = clusters.filter(isLone);
  if (lone.length < 2)
    return clusters;

  const groups = await similarityGroups(lone.map((c) => c.series[0].coreSmiles), threshold, fingerprintType);
  const out = clusters.filter((c) => !isLone(c));
  for (const group of groups)
    out.push({id: `p${out.length}`, series: group.map((i) => lone[i].series[0]), siteKey: '', level: 2});
  return out;
}

/**
 * Group the compounds that never formed a matched series at all.
 *
 * {@link poolUngroupedSeries} can only rescue a series that exists; a compound whose every cut left it
 * without a partner has no core, so nothing upstream ever offers it to a matrix. Singletons are
 * dropped — one compound is no use here.
 *
 * Marked `requiresDecomposition`: the series below are placeholders carrying no real core, so without
 * a decomposition the single-position fallback would read each compound as its own bare core.
 */
export async function poolUngroupedMolecules(clusters: CoreCluster[], molecules: string[],
  threshold: number, fingerprintType: Fingerprint = Fingerprint.Morgan): Promise<CoreCluster[]> {
  const covered = new Set<number>();
  for (const cluster of clusters) {
    for (const matched of cluster.series) {
      for (const member of matched.members)
        covered.add(member.molIdx);
    }
  }
  const leftover = molecules.map((_m, i) => i).filter((i) => !covered.has(i) && molecules[i]);
  if (leftover.length < 2)
    return clusters;

  const groups = await similarityGroups(leftover.map((i) => molecules[i]), threshold, fingerprintType);
  const out = [...clusters];
  for (const group of groups.filter((g) => g.length > 1)) {
    out.push({
      id: `m${out.length}`,
      series: group.map((i) => ({
        coreSmiles: molecules[leftover[i]],
        members: [{molIdx: leftover[i], substSmiles: ''}],
      })),
      siteKey: '', level: 2, requiresDecomposition: true,
    });
  }
  return out;
}
