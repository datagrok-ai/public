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
 * Step 4 — group series into matrix rows by cutting their cores a second time and grouping on the
 * remainder ("key of key"). This is exact and order-independent: cores differing at one point leave
 * the same remainder, with no fingerprint or threshold. A core joins one group per cuttable site (an
 * overlapping cover), and groups fully contained in a larger one are dropped.
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

  // Largest first, so a group survives only when it covers series the kept ones do not.
  const groups = [...bySite.entries()].filter(([, g]) => g.size > 1).sort((a, b) => b[1].size - a[1].size);
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

  let assignment: number[];
  if (keys.length <= MAX_BUTINA_CORES)
    assignment = await butinaCluster(keys, threshold, fingerprintType);
  else {
    // Past MAX_BUTINA_CORES the quadratic all-pairs pass costs too much; BitBIRCH builds its tree in
    // one incremental sweep instead.
    const coreCol = DG.Column.fromStrings('core', keys);
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    const clusterCol = await bitbirchWorker(coreCol, threshold, fingerprintType);
    assignment = series.map((_s, i) => clusterCol.isNone(i) ? -1 : clusterCol.get(i));
  }

  // Group series by cluster id; anything unplaced becomes its own singleton.
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
  for (const group of seriesByCluster.values())
    clusters.push({id: `c${n++}`, series: group, siteKey: '', level: 2});
  return clusters;
}
