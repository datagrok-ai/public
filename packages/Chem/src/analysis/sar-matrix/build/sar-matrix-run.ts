import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../../package';
import {getRdKitModule} from '../../../utils/chem-common-rdkit';
import {isCruxRequested} from '../../../engine/chem-engine';
import {convertMolNotationByEngine} from '../../../engine/worker-jobs';
import {getMmpFrags} from '../../molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {MmpFragments} from '../../molecular-matched-pairs/mmp-analysis/mmpa-misc';
import {SCALING_METHODS} from '../../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {assembleMultiPositionMatrix, linkStaged} from './sar-matrix-assemble';
import {canonicalSmiles, checkAgainstMolecules, cutWarning, decomposeByColumns, SarFragmentColumns, SeriesCut}
  from './sar-matrix-columns';
import {cachedAtomCount, ClusterDecomposition, decomposeClusters} from './sar-matrix-decompose';
import {computeMatrixConfidence} from './sar-matrix-confidence';
import {buildMatchedSeries, buildCoarserLevels, clusterRelatedCores, groupSeriesByColumn, groupSeriesBySite,
  poolUngroupedSeries, poolUngroupedMolecules}
  from './sar-matrix-clustering';
import {rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {CoreCluster, logSarTime, observedMolecules, SarMatrix, SarMatrixCell} from '../sar-matrix-types';

/** Joins a core to one substituent at its single attachment point. */
const SINGLE_CUT = [{R1: [1]}];

const STEREO_MARK = /[@/\\]/;

/**
 * Link each virtual analog's row core to its column substituent in one batched worker call.
 *
 * Single-cut matrices only (empty `refValues`): decomposed matrices are linked during assembly.
 */
export async function linkVirtualStructures(matrices: SarMatrix[]): Promise<void> {
  const cores: string[] = [];
  const fragments: string[] = [];
  const targets: SarMatrixCell[] = [];

  for (const matrix of matrices) {
    if (Object.keys(matrix.refValues).length > 0)
      continue;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      for (let ci = 0; ci < matrix.columns.length; ci++) {
        const cell = matrix.cells[ri][ci];
        if (cell.kind === 'virtual' && cell.smiles === null) {
          cores.push(matrix.rows[ri].coreSmiles);
          fragments.push(matrix.columns[ci].substSmiles);
          targets.push(cell);
        }
      }
    }
  }
  if (targets.length === 0)
    return;

  const linked = await linkStaged(cores, cores.map(() => SINGLE_CUT), (i) => fragments[i]);
  for (let i = 0; i < targets.length; i++)
    targets[i].smiles = linked[i];
}

/** The component with the most heavy atoms, canonical: fragmenting a salt carries the counter-ion on a
 *  substituent. Counted in atoms, since a counter-ion such as citrate can be the longer string. */
function parentStructure(smiles: string): string {
  if (!smiles.includes('.'))
    return smiles;
  const rdkit = getRdKitModule();
  const size = (component: string): number => cachedAtomCount(component, rdkit);
  return canonicalSmiles(smiles.split('.')
    .reduce((a, b) => size(b) > size(a) || (size(b) === size(a) && b > a) ? b : a));
}

/**
 * Take back the predictions whose structure the table already holds under another cut, core spelling, salt
 * form or series: `unmeasured` when never assayed and in this series, else `excluded`. No measurement moves.
 * One without stereo that the table holds with stereo is `excluded` but names no compound: the cut lost the
 * configuration, so it may be any of the isomers.
 */
function claimHeldStructures(matrices: SarMatrix[], molecules: string[], table: string[], assayed: Uint8Array,
  assigned: (string | null)[] | null, predictUnmeasured: boolean): void {
  const exact = new Map<string, number>();
  const parents = new Map<string, number>();
  const stereoFree = new Set<string>();
  const hold = (held: Map<string, number>, key: string, i: number): void => {
    const prior = held.get(key);
    // A repeated compound stands for its assayed copy.
    if (prior === undefined || (assayed[i] && !assayed[prior]))
      held.set(key, i);
  };
  table.forEach((smiles, i) => {
    if (!smiles)
      return;
    const parent = parentStructure(smiles);
    hold(exact, smiles, i);
    hold(parents, parent, i);
    if (STEREO_MARK.test(parent))
      stereoFree.add(canonicalSmiles(parent.replace(/@+/g, '').replace(/[/\\]/g, '')));
  });
  for (const matrix of matrices) {
    for (const row of matrix.cells) {
      row.forEach((cell, ci) => {
        if (cell.kind !== 'virtual' || cell.smiles === null)
          return;
        const parent = parentStructure(cell.smiles);
        // The exact form first: of a compound held both bare and as a salt, a proposal is the form it spells.
        const molIdx = exact.get(cell.smiles) ?? parents.get(parent);
        if (molIdx === undefined) {
          if (!STEREO_MARK.test(parent) && stereoFree.has(parent)) {
            matrix.virtualCount--;
            row[ci] = {kind: 'excluded', value: null, molIdx: null, smiles: null};
          }
          return;
        }
        matrix.virtualCount--;
        const own = assigned === null || assigned[molIdx] === matrix.label;
        if (own && !assayed[molIdx]) {
          row[ci] = predictUnmeasured ? {...cell, kind: 'unmeasured', molIdx, smiles: molecules[molIdx]} :
            {kind: 'unmeasured', value: null, molIdx, smiles: molecules[molIdx]};
        } else {
          row[ci] = {kind: 'excluded', value: null, molIdx: own ? molIdx : null,
            smiles: own ? molecules[molIdx] : null};
        }
      });
    }
  }
}

/** Each molecule column's canonical structures at the version read, since switching the column component
 *  rebuilds over the same molecules. Keyed by the Dart column: every read from the table is a new wrapper. By the
 *  engine in effect, the one that writes the linked structures they are compared with. */
const canonicalCache = new WeakMap<object, {version: number, crux: boolean, structures: Promise<string[]>}>();

function canonicalStructures(molecules: DG.Column, molList: string[]): Promise<string[]> {
  const crux = isCruxRequested();
  const cached = canonicalCache.get(molecules.dart);
  if (cached?.version === molecules.version && cached.crux === crux)
    return cached.structures;
  const structures = convertMolNotationByEngine(molList, DG.chem.Notation.Smiles, false);
  canonicalCache.set(molecules.dart, {version: molecules.version, crux, structures});
  structures.catch(() => {
    if (canonicalCache.get(molecules.dart)?.structures === structures)
      canonicalCache.delete(molecules.dart);
  });
  return structures;
}

/** Drops the cuts whose two pieces do not rebuild their own compound: the fragmenter can write a stereocentre
 *  at the cut mirrored. A cut that cannot be joined at all is kept, as its analogs fail anyway. */
async function dropMirroredCuts(frags: MmpFragments, parents: string[]): Promise<number> {
  const cuts: {row: number, k: number}[] = [];
  const cores: string[] = [];
  const values: string[] = [];
  frags.fragCodes.forEach((codes, row) => codes.forEach(([core, value], k) => {
    const [a, b] = [frags.idToName[core], frags.idToName[value]];
    if (STEREO_MARK.test(a) || STEREO_MARK.test(b)) {
      cuts.push({row, k});
      cores.push(a);
      values.push(b);
    }
  }));
  if (cuts.length === 0)
    return 0;
  const rebuilt = await linkStaged(cores, cores.map(() => SINGLE_CUT), (i) => values[i]);
  const dropped = new Map<number, Set<number>>();
  rebuilt.forEach((smiles, i) => {
    const parent = parents[cuts[i].row];
    if (smiles === null || smiles === parent || canonicalSmiles(smiles) === parent)
      return;
    if (!dropped.has(cuts[i].row))
      dropped.set(cuts[i].row, new Set());
    dropped.get(cuts[i].row)!.add(cuts[i].k);
  });
  let count = 0;
  for (const [row, ks] of dropped) {
    frags.fragCodes[row] = frags.fragCodes[row].filter((_code, k) => !ks.has(k));
    count += ks.size;
  }
  return count;
}

/** Row cores, column substituents and observed compounds of an assembled matrix, as comparable sets. */
function matrixContent(matrix: SarMatrix): {rows: Set<string>, cols: Set<string>, mols: Set<number>} {
  return {
    rows: new Set(matrix.rows.map((r) => r.keySmiles)),
    cols: new Set(matrix.columns.map((c) => `${c.position}${c.substSmiles}`)),
    mols: observedMolecules(matrix),
  };
}

/**
 * Drop matrices whose rows, columns and compounds are all subsets of a larger one. All three must
 * match: a different core/substituent split over the same compounds is a distinct view, not a crop.
 * Runs on assembled matrices, since two differently-keyed groups can assemble into the same table.
 */
function dropCroppedMatrices(matrices: SarMatrix[]): SarMatrix[] {
  const content = matrices.map(matrixContent);
  // Position-independent identity: which of two equal matrices survives follows what they hold, not
  // where they sit in the array.
  const rank = matrices.map((matrix, i) => [
    [...content[i].rows].sort().join(''),
    [...content[i].cols].sort().join(''),
    [...content[i].mols].sort((x, y) => x - y).join(','),
    matrix.level,
    matrix.siteKey,
  ].join(''));
  const subset = <T>(a: Set<T>, b: Set<T>): boolean => a.size <= b.size && [...a].every((x) => b.has(x));
  return matrices.filter((_matrix, i) => !matrices.some((_other, j) => {
    if (i === j)
      return false;
    const a = content[i];
    const b = content[j];
    if (matrices[i].level !== matrices[j].level) {
      // A coarser matrix containing finer ones is the hierarchy, not redundancy; drop it only when it
      // restates a finer matrix exactly, keeping the finer one.
      const identical = a.rows.size === b.rows.size && a.cols.size === b.cols.size &&
        a.mols.size === b.mols.size;
      if (!identical || matrices[i].level < matrices[j].level)
        return false;
    }
    const smaller = a.rows.size < b.rows.size || a.cols.size < b.cols.size || a.mols.size < b.mols.size;
    // Of two equal matrices exactly one must go, or they eliminate each other. The index only settles
    // pairs alike in content, level and key, which are interchangeable anyway.
    const tie = rank[j] < rank[i] || (rank[j] === rank[i] && j < i);
    if (!smaller && !(a.rows.size === b.rows.size && a.cols.size === b.cols.size && tie))
      return false;
    return subset(a.rows, b.rows) && subset(a.cols, b.cols) && subset(a.mols, b.mols);
  }));
}

/** Number the series by their place in the fold hierarchy: "Series 3_2" is the second matrix inside
 *  "Series 3". Walked from the roots so a parent is always labelled before its children; a matrix
 *  whose parent did not survive is itself a root. */
function assignSeriesLabels(matrices: SarMatrix[]): void {
  const present = new Set(matrices.map((matrix) => matrix.id));
  const children = new Map<string, SarMatrix[]>();
  const roots: SarMatrix[] = [];
  for (const matrix of matrices) {
    const parentId = matrix.parentId !== undefined && present.has(matrix.parentId) ? matrix.parentId : null;
    if (parentId === null)
      roots.push(matrix);
    else
      children.set(parentId, [...(children.get(parentId) ?? []), matrix]);
  }
  const label = (matrix: SarMatrix, path: string): void => {
    matrix.label = `Series ${path}`;
    (children.get(matrix.id) ?? []).forEach((child, i) => label(child, `${path}_${i + 1}`));
  };
  roots.forEach((matrix, i) => label(matrix, `${i + 1}`));
}

/** The user's series value per compound, or null when no column was given. Read through the column's
 *  own formatting so a numeric column groups as "1", "2", "3" rather than "1.0000" — the names must
 *  match the cells. Empty and missing values stay null and are left ungrouped. */
function readSeriesValues(column: DG.Column | null, rowCount: number): (string | null)[] | null {
  if (column === null)
    return null;
  const values: (string | null)[] = new Array(rowCount).fill(null);
  for (let i = 0; i < rowCount; i++) {
    if (column.isNone(i))
      continue;
    const text = column.getString(i).trim();
    values[i] = text === '' ? null : text;
  }
  return values;
}

/** How series become the rows of one matrix. */
export enum SarGrouping {
  /** Cut the cores a second time and group the ones that leave the same remainder. */
  Site = 'Site',
  /** Cluster the cores by fingerprint similarity. */
  Similarity = 'Similarity',
}

/** Ceiling on {@link SarMatrixParams.fragmentationLevels}: each level is another fragmentation pass,
 *  and deep tiers mostly restate the ones below (a key reduced to a bare ring cannot strip further). */
export const MAX_SERIES_LEVELS = 5;

export interface SarMatrixParams {
  scaling: SCALING_METHODS;
  fragmentCutoff: number;
  predictVirtual: boolean;
  grouping: SarGrouping;
  /** Tiers of matrices, counted as the UI numbers them (3 => navigator shows L1, L2, L3). Each tier
   *  past 1 relates whole matrices by cutting their keys again, giving fewer and broader ones. */
  fragmentationLevels: number;
  /** Whether higher scaled activity is more potent. Passed in rather than re-derived from `scaling`
   *  so ranking and the cards' best/mean text agree with the direction the user set. */
  higherIsBetter: boolean;
  /** Core-similarity clustering threshold (lower groups more distant cores). Similarity grouping only. */
  threshold: number;
  /** Fewest MEASURED compounds a matrix must hold to be kept. One measured compound gives no row
   *  effect, no column effect, no correlation and nothing to predict from, so such a matrix is a row
   *  of blanks around a single value. Matrices are grouped by structure before activity is attached,
   *  which is how they arise: the compounds are there, the measurements are not. */
  minCompounds: number;
  /** Whether to predict the activity of a compound the set already holds but has no measurement for.
   *  Same additive model as an unexplored cell, and only where its row and column both have
   *  measurements to fit — it just never becomes something to synthesize. */
  predictUnmeasured: boolean;
  /** Optional: series the user assigned themselves, one value per compound. Replaces `grouping`
   *  entirely and skips the automatic tiers, which would otherwise fold the user's series by
   *  chemistry into matrices they did not ask for. */
  seriesColumn?: DG.Column | null;
  /** Optional: core and R-group columns already in the table; they replace fragmentation and grouping. */
  fragmentColumns?: SarFragmentColumns | null;
  /** Whether to also cover the compounds no shared core could group, by pooling those series and
   *  searching for a common core with an MCS. No matrix is lost and no compound is dropped by turning
   *  it on, but a cluster whose site key is too small to anchor also falls through to the MCS, so a
   *  core-grouped matrix can come out multi-position instead of single-position. Off leaves the
   *  ungrouped compounds out entirely, which is also the only setting whose output cannot vary between
   *  runs — an MCS that overruns retires its worker, and the clusters queued behind it go
   *  undecomposed. */
  useMcsAnchors: boolean;
  rankScheme: SarRankScheme;
}


/**
 * Scaled activity per row, NaN where missing or unscalable: a null sentinel read as a number would poison
 * the Free-Wilson fit. `assayed` marks rows holding a value at all, usable or not: a log scale turns a real
 * zero or negative into NaN, and such a compound must not end up labelled untested.
 */
export function readActivities(activity: DG.Column<number>, scaling: SCALING_METHODS):
  {activities: Float32Array, assayed: Uint8Array, unscalable: number} {
  const scaledCol = scaling === SCALING_METHODS.NONE ? activity : scaleActivity(activity, scaling);
  // Bounded by row count, since column storage has spare capacity past the last row.
  const scaled = scaledCol.getRawData();
  const activities = new Float32Array(activity.length);
  const assayed = new Uint8Array(activity.length);
  let unscalable = 0;
  for (let i = 0; i < activities.length; i++) {
    if (activity.isNone(i)) {
      activities[i] = NaN;
      continue;
    }
    assayed[i] = 1;
    // The log of a zero or a negative is NaN, which the scaled column stores as its null: a finite number near 0.
    activities[i] = scaledCol !== activity && scaledCol.isNone(i) ? NaN : scaled[i];
    if (!Number.isFinite(activities[i]))
      unscalable++;
  }
  return {activities, assayed, unscalable};
}

/**
 * Run the full SAR Matrix pipeline: fragment molecules, build matched series, group related cores,
 * assemble a potency matrix per group (predicting virtual analogs), and return them ranked.
 */
export async function runSarMatrix(molecules: DG.Column, activity: DG.Column<number>,
  params: SarMatrixParams): Promise<SarMatrix[]> {
  const tTotal = performance.now();
  const molList = molecules.toList();
  const {activities, assayed, unscalable} = readActivities(activity, params.scaling);
  if (unscalable > 0) {
    // Without this the whole analysis can come back empty with nothing to explain it: an all-negative
    // column (ΔG, say) under the default −lg scaling leaves no observed cell anywhere.
    const message = `SAR Matrix: ${unscalable} of ${activity.length} "${activity.name}" values cannot be ` +
      `scaled by ${params.scaling}, which needs positive numbers. They are excluded from every matrix — ` +
      'set Scaling to "none" to use them as they are.';
    _package.logger.warning(message);
    grok.shell.warning(message);
  }

  // A column of user-assigned series wins over both automatic groupings.
  const assigned = readSeriesValues(params.seriesColumn ?? null, molecules.length);
  let t = performance.now();
  let clusters: CoreCluster[];
  let decomps: (ClusterDecomposition | null)[];
  let cuts: SeriesCut[] = [];
  // Fragmentation computes the table's canonical structures anyway; read from columns, they are fetched alongside.
  let canonical: Promise<string[]> | null = null;
  if (params.fragmentColumns) {
    canonical = canonicalStructures(molecules, molList);
    const read = decomposeByColumns(params.fragmentColumns, assigned, activities);
    ({clusters, decomps, cuts} = read);
    const {checked, mismatched} = await checkAgainstMolecules(read.decomps, molList);
    if (mismatched > 0 && mismatched * 10 >= checked) {
      const message = `SAR Matrix: ${mismatched} of ${checked} compounds checked are not what their core and ` +
        `R-groups build; the columns may not belong to "${molecules.name}".`;
      _package.logger.warning(message);
      grok.shell.warning(message);
    }
    logSarTime(`fragment columns (${clusters.length} groups, ${mismatched}/${checked} mismatched)`, t);
  } else {
    const tFrag = performance.now();
    // Read as the fragmenter reads them, since one it cannot read fails the whole call, and cut to the
    // parent, since a counter-ion would ride on a substituent. A row that cannot be read sits out.
    const structures: string[] = [];
    const parents: string[] = [];
    for (let i = 0; i < molList.length; i++) {
      // Read in slices, so a large table does not hold the page.
      if (i % 500 === 499)
        await new Promise((resolve) => setTimeout(resolve, 0));
      structures.push(canonicalSmiles(molList[i]));
      parents.push(parentStructure(structures[i]));
    }
    const readable = parents.flatMap((parent, i) => parent ? [i] : []);
    const [read] = await getMmpFrags(readable.map((i) => parents[i]));
    const fragCodes: [number, number][][] = parents.map(() => []);
    readable.forEach((row, k) => fragCodes[row] = read.fragCodes[k]);
    const frags = {...read, fragCodes};
    canonical = Promise.resolve(structures);
    logSarTime(`MMP fragmentation (${readable.length} molecules)`, tFrag);
    t = performance.now();
    const mirrored = await dropMirroredCuts(frags, parents);
    logSarTime(`cuts that do not rebuild their compound (${mirrored} dropped)`, t);
    t = performance.now();
    const series = buildMatchedSeries(frags, params.fragmentCutoff);
    logSarTime(`matched series (${series.length} series)`, t);
    t = performance.now();
    let base = assigned !== null ? await groupSeriesByColumn(series, assigned) :
      params.grouping === SarGrouping.Site ?
        await groupSeriesBySite(series) :
        await clusterRelatedCores(series, params.threshold);
    // A series no core grouping could place is alone in its cluster, so it holds one row and never
    // reaches a matrix. Pooling those is what the MCS is for, and the user asks for it: unasked,
    // compounds without a shared core are left out rather than grouped by a search.
    if (params.useMcsAnchors && assigned === null) {
      base = await poolUngroupedSeries(base, params.threshold);
      base = await poolUngroupedMolecules(base, parents, params.threshold);
    }
    logSarTime(`grouping by ${assigned !== null ? 'series column' : params.grouping.toLowerCase()} ` +
      `(${base.length} groups)`, t);
    // Clamped, not trusted: a programmatic caller can pass any number and each level costs a pass. The
    // user's own series are never folded further: the tiers above them would mix series they named.
    const levels = assigned !== null ? 1 : Math.min(MAX_SERIES_LEVELS, Math.max(1, params.fragmentationLevels));
    t = performance.now();
    clusters = await buildCoarserLevels(base, levels - 1);
    logSarTime(`coarser levels (${clusters.length} clusters)`, t);

    // Decompose all clusters in one batched pass before assembly so they run parallel across workers;
    // per-cluster calls under the Promise.all below would serialize on one worker. A null decomposition
    // falls back to single-position construction inside the assembler.
    const clusterMembers = clusters.map((c) => [...new Set(c.series.flatMap((s) => s.members.map((m) => m.molIdx)))]);
    // The grouping already knows each cluster's shared scaffold; handing it over is what keeps the MCS
    // to the clusters that have none.
    const clusterSiteKeys = clusters.map((c) => c.siteKey);
    t = performance.now();
    decomps = await decomposeClusters(clusterMembers, parents, params.useMcsAnchors, clusterSiteKeys);
    logSarTime(`decomposition total (${decomps.filter(Boolean).length}/${clusters.length} clusters decomposed)`, t);
  }

  // A matrix needs >=2 rows to compare cores; folding can collapse several series onto one row, so
  // this is checked after assembly. The floor counts MEASURED compounds, not clustered members.
  t = performance.now();
  const usable = (matrix: SarMatrix | null): matrix is SarMatrix =>
    matrix !== null && matrix.realCount >= params.minCompounds && matrix.columns.length > 0 &&
    matrix.rows.length >= 2;
  const axis = params.fragmentColumns?.column.name;
  const assembled = (await Promise.all(clusters.map((cluster, i) =>
    assembleMultiPositionMatrix(cluster, molList, activities, params.predictVirtual, decomps[i],
      params.predictUnmeasured, assayed, axis))))
    .filter(usable);

  // An anchor can decompose a cluster into fewer rows than it has series, and a matrix under two rows
  // is dropped — so a cluster that yields a matrix without an anchor could lose it to having one. The
  // MCS may only add, so those are rebuilt the way they would have been built without it; a cluster
  // whose series are placeholders returns null from that build and simply stays absent.
  const built = new Set(assembled.map((matrix) => matrix.id));
  // Columns mode has no anchor to cost anything, and its clusters carry no series to rebuild from.
  const rescued = params.fragmentColumns ? [] : (await Promise.all(clusters.map((cluster, i) =>
    decomps[i] !== null && !built.has(cluster.id) ?
      assembleMultiPositionMatrix(cluster, molList, activities, params.predictVirtual, null,
        params.predictUnmeasured, assayed) : null)))
    .filter(usable);
  if (rescued.length) {
    assembled.push(...rescued);
    _package.logger.debug(`SAR Matrix: rebuilt ${rescued.length} matrices the anchor would have cost`);
  }
  logSarTime(`assembly (${assembled.length} matrices)`, t);
  const matrices = dropCroppedMatrices(assembled);
  if (matrices.length < assembled.length)
    _package.logger.debug(`SAR Matrix: dropped ${assembled.length - matrices.length} cropped matrices`);

  // Re-link each survivor to its nearest surviving ancestor so dropping a matrix does not orphan the
  // ones below it. The chain is read from clusters, which hold links through groups that never assembled.
  const clusterParent = new Map(clusters.map((cluster) => [cluster.id, cluster.parentId]));
  const survived = new Set(matrices.map((matrix) => matrix.id));
  for (const matrix of matrices) {
    let parent = clusterParent.get(matrix.id);
    while (parent !== undefined && !survived.has(parent))
      parent = clusterParent.get(parent);
    matrix.parentId = parent;
  }

  if (params.fragmentColumns) {
    const linksOf = new Map(clusters.map((cluster, i) => [cluster.id, decomps[i]?.links]));
    let collapsed = 0;
    for (const matrix of matrices) {
      const structures = linksOf.get(matrix.id)?.structures;
      let named = false;
      matrix.rows.forEach((row, i) => {
        // Names are what tell such rows apart; structures are drawn.
        const names = Object.values(row.foldedValues).filter((v) => v !== '' && !structures?.has(v));
        named ||= names.length > 0;
        row.label = names.length > 0 ? names.join(' · ') : `Row ${i + 1}`;
      });
      if (!named && matrix.rows.length > 1 && new Set(matrix.rows.map((row) => row.keySmiles)).size === 1)
        collapsed++;
    }
    if (collapsed > 0) {
      const message = `SAR Matrix: the rows of ${collapsed} series look alike, because their R-groups do ` +
        'not join onto the core.';
      _package.logger.warning(message);
      grok.shell.warning(message);
    }
  }
  // Labels assigned before ranking so re-ranking never renames a matrix. Series the user grouped
  // themselves already carry their own name and must keep it.
  if (assigned === null)
    assignSeriesLabels(matrices);
  const cutMessage = cutWarning(cuts, (id) => matrices.find((matrix) => matrix.id === id)?.label);
  if (cutMessage !== null) {
    _package.logger.warning(cutMessage);
    grok.shell.warning(cutMessage);
  }
  // Started before the confidence pass, which is main-thread and has no await, while linking is a worker
  // round-trip. They touch disjoint fields: the fit writes `confidence`, linking `smiles` on virtual cells.
  t = performance.now();
  const linking = params.predictVirtual ? linkVirtualStructures(matrices) : null;
  matrices.forEach((matrix) => matrix.confidence = computeMatrixConfidence(matrix));
  logSarTime('confidence (leave-one-out)', t);
  if (linking !== null) {
    t = performance.now();
    await linking;
    logSarTime('virtual structure linking (overlapped)', t);
    t = performance.now();
    let table: string[] | null = null;
    try {
      table = await canonical!;
    } catch (e) {
      // The claim only corrects proposals, so the matrices are worth more than failing the run over it.
      _package.logger.warning(`SAR Matrix: the structures of "${molecules.name}" could not be canonicalized, ` +
        `so proposals the table already holds are not taken back: ${e instanceof Error ? e.message : e}`);
    }
    if (table !== null)
      claimHeldStructures(matrices, molList, table, assayed, assigned, params.predictUnmeasured);
    logSarTime('proposals the table already holds', t);
  }

  const ranked = rankMatrices(matrices, params.rankScheme, params.higherIsBetter);
  logSarTime(`TOTAL (${ranked.length} matrices)`, tTotal);
  return ranked;
}
