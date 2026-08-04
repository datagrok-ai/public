import * as DG from 'datagrok-api/dg';

import {Fingerprint} from '../../utils/chem-common';
import {bitbirchWorker} from '../bit-birch/bitbirch-clustering';
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
      members.push({molIdx, substFragId: substId, substSmiles: idToName[substId]});
    series.push({coreFragId, coreSmiles: idToName[coreFragId], members});
  }
  return series;
}

/**
 * Step 4 — cluster structurally related series into the rows of one matrix.
 *
 * Series are grouped by the fingerprint similarity of their **Bemis-Murcko scaffold**
 * (ring systems + linkers, computed server-side), reusing BitBIRCH. Clustering on the
 * scaffold rather than the raw single-cut core gives homogeneous clusters — every
 * molecule in a cluster shares a real, decomposable core — so the downstream R-group
 * decomposition stays clean instead of producing garbage fragments. Similarity (not
 * identity) still groups ring/scaffold hops (pyrazole vs thiazole) into one matrix,
 * which is what SAR transfer between rows depends on.
 *
 * When `scaffolds` is empty (the Murcko script/compute environment is unavailable),
 * this falls back to clustering on the single-cut core SMILES.
 *
 * MACCS keys are the default because their substructure bits group heteroaromatic ring
 * hops better than Morgan circular fingerprints. The threshold is looser than the
 * BitBIRCH default so moderately similar scaffolds still merge; both are exposed for tuning.
 *
 * @param series Matched series from `buildMatchedSeries`.
 * @param scaffolds Murcko scaffold SMILES per molecule (by global index); empty to fall back.
 * @param threshold BitBIRCH merge threshold; lower groups more distant scaffolds.
 * @param fingerprintType Fingerprint used for scaffold similarity.
 * @returns One `CoreCluster` per group of related cores.
 */
export async function clusterRelatedCores(series: MatchedSeries[], scaffolds: string[] = [],
  threshold: number = 0.5, fingerprintType: Fingerprint = Fingerprint.MACCS): Promise<CoreCluster[]> {
  if (series.length === 0)
    return [];
  if (series.length === 1)
    return [{id: 'c0', series: [series[0]]}];

  // Representative structure per series: its Murcko scaffold when available (all members of a
  // single-cut series share the same ring system), else the single-cut core.
  const keyOf = (s: MatchedSeries): string =>
    (scaffolds.length ? scaffolds[s.members[0].molIdx] : '') || s.coreSmiles;
  const coreCol = DG.Column.fromStrings('core', series.map(keyOf));
  coreCol.semType = DG.SEMTYPE.MOLECULE;
  const clusterCol = await bitbirchWorker(coreCol, threshold, fingerprintType);

  // Group series by cluster id; cores BitBIRCH could not place get their own singleton.
  const seriesByCluster = new Map<number, MatchedSeries[]>();
  let singletonKey = -1;
  for (let i = 0; i < series.length; i++) {
    const cluster = clusterCol.isNone(i) ? singletonKey-- : clusterCol.get(i);
    let group = seriesByCluster.get(cluster);
    if (!group)
      seriesByCluster.set(cluster, group = []);
    group.push(series[i]);
  }

  let n = 0;
  const clusters: CoreCluster[] = [];
  for (const group of seriesByCluster.values())
    clusters.push({id: `c${n++}`, series: group});
  return clusters;
}
