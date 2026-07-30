// Pure comparison logic: no DG/platform dependencies, unit-testable.
// See docs/run-comparison-tool.md for the rules implemented here.

export type MatchConfidence = 'exact' | 'normalized' | 'fuzzy';
export type UnitsCompatibility = 'match' | 'warn' | 'mismatch';

export const FUZZY_NAME_THRESHOLD = 0.7;

export interface ScalarNodeInfo {
  path: string;
  name: string;
  // human-readable path: step friendly names + IO caption; `path` stays the stable id
  friendlyPath?: string;
  valueType: string;
  units?: string;
  value: number | null;
}

export interface ColumnInfo {
  name: string;
  type: string;
  units?: string;
}

export interface TableNodeInfo {
  path: string;
  name: string;
  friendlyPath?: string;
  // nqName of the producing function; used to merge same-function tables in index selection
  nqName?: string;
  columns: ColumnInfo[];
  rowCount: number;
}

export interface ComparisonEntryNodes {
  entryId: string;
  entryName: string;
  scalars: ScalarNodeInfo[];
  tables: TableNodeInfo[];
}

export interface ScalarBinding {
  entryId: string;
  path: string;
  name: string;
  friendlyPath?: string;
  units?: string;
  value: number | null;
}

export interface ColumnBinding {
  entryId: string;
  tablePath: string;
  tableName: string;
  tableFriendlyPath?: string;
  columnName: string;
  units?: string;
  indexColumnName: string;
  splitColumnName?: string;
}

export interface TargetBase {
  key: string;
  displayName: string;
  confidence: MatchConfidence;
  unitsWarning: boolean;
  coverage: number;
  total: number;
}

export interface ScalarTarget extends TargetBase {
  kind: 'scalar';
  bindings: ScalarBinding[];
}

export interface ColumnTarget extends TargetBase {
  kind: 'column';
  bindings: ColumnBinding[];
}

export type ComparisonTarget = ScalarTarget | ColumnTarget;

const NUMERIC_TYPES = new Set(['int', 'double', 'float', 'number', 'bigint']);

export function isNumericType(type: string): boolean {
  return NUMERIC_TYPES.has(type);
}

export function normalizeName(name: string): string {
  return name.trim().toLowerCase().replace(/[\s_\-]+/g, ' ');
}

function bigrams(s: string): Map<string, number> {
  const res = new Map<string, number>();
  for (let i = 0; i < s.length - 1; i++) {
    const gram = s.slice(i, i + 2);
    res.set(gram, (res.get(gram) ?? 0) + 1);
  }
  return res;
}

/** Sørensen–Dice bigram similarity of normalized names, in [0, 1]. */
export function nameSimilarity(a: string, b: string): number {
  const na = normalizeName(a);
  const nb = normalizeName(b);
  if (na === nb)
    return 1;
  if (na.length < 2 || nb.length < 2)
    return 0;
  const ga = bigrams(na);
  const gb = bigrams(nb);
  let intersection = 0;
  let totalA = 0;
  let totalB = 0;
  for (const [gram, count] of ga) {
    totalA += count;
    intersection += Math.min(count, gb.get(gram) ?? 0);
  }
  for (const count of gb.values())
    totalB += count;
  return (2 * intersection) / (totalA + totalB);
}

export function nameMatchConfidence(a: string, b: string): MatchConfidence | null {
  if (a === b)
    return 'exact';
  if (normalizeName(a) === normalizeName(b))
    return 'normalized';
  if (nameSimilarity(a, b) >= FUZZY_NAME_THRESHOLD)
    return 'fuzzy';
  return null;
}

const CONFIDENCE_RANK: Record<MatchConfidence, number> = {exact: 0, normalized: 1, fuzzy: 2};

export function weakerConfidence(a: MatchConfidence, b: MatchConfidence): MatchConfidence {
  return CONFIDENCE_RANK[a] >= CONFIDENCE_RANK[b] ? a : b;
}

export function unitsCompatibility(a?: string, b?: string): UnitsCompatibility {
  const ua = a?.trim().toLowerCase() ?? '';
  const ub = b?.trim().toLowerCase() ?? '';
  if (ua === ub)
    return 'match';
  if (!ua || !ub)
    return 'warn';
  return 'mismatch';
}

interface ClusterItem<P> {
  entryId: string;
  name: string;
  units?: string;
  // used only to break ties between equally-confident candidates (e.g. table name for columns)
  secondaryName?: string;
  payload: P;
}

interface Cluster<P> {
  canonicalName: string;
  canonicalSecondary?: string;
  confidence: MatchConfidence;
  unitsWarning: boolean;
  items: ClusterItem<P>[];
  entryIds: Set<string>;
}

function clusterByName<P>(items: ClusterItem<P>[]): Cluster<P>[] {
  const clusters: Cluster<P>[] = [];
  for (const item of items) {
    let best: {cluster: Cluster<P>, confidence: MatchConfidence, score: number} | null = null;
    for (const cluster of clusters) {
      if (cluster.entryIds.has(item.entryId))
        continue;
      const confidence = nameMatchConfidence(cluster.canonicalName, item.name);
      if (!confidence)
        continue;
      if (unitsCompatibility(cluster.items[0].units, item.units) === 'mismatch')
        continue;
      const score = nameSimilarity(cluster.canonicalName, item.name) +
        (cluster.canonicalSecondary && item.secondaryName ?
          nameSimilarity(cluster.canonicalSecondary, item.secondaryName) * 0.5 : 0);
      if (!best || CONFIDENCE_RANK[confidence] < CONFIDENCE_RANK[best.confidence] ||
        (confidence === best.confidence && score > best.score))
        best = {cluster, confidence, score};
    }
    if (best) {
      const cluster = best.cluster;
      cluster.items.push(item);
      cluster.entryIds.add(item.entryId);
      cluster.confidence = weakerConfidence(cluster.confidence, best.confidence);
      if (unitsCompatibility(cluster.items[0].units, item.units) === 'warn')
        cluster.unitsWarning = true;
    } else {
      clusters.push({
        canonicalName: item.name,
        canonicalSecondary: item.secondaryName,
        confidence: 'exact',
        unitsWarning: false,
        items: [item],
        entryIds: new Set([item.entryId]),
      });
    }
  }
  return clusters;
}

// distinct clusters can share a canonical name (a cluster takes at most one item per entry),
// and duplicate keys corrupt keyed list rendering
function dedupeTargetKeys<T extends TargetBase>(targets: T[]): T[] {
  const seen = new Map<string, number>();
  for (const target of targets) {
    const count = seen.get(target.key) ?? 0;
    seen.set(target.key, count + 1);
    if (count > 0)
      target.key = `${target.key}:${count + 1}`;
  }
  return targets;
}

/** Groups numeric scalars across entries into candidate targets (coverage >= 2). */
export function matchScalarTargets(entries: ComparisonEntryNodes[]): ScalarTarget[] {
  const items: ClusterItem<ScalarBinding>[] = [];
  for (const entry of entries) {
    for (const scalar of entry.scalars) {
      if (!isNumericType(scalar.valueType))
        continue;
      items.push({
        entryId: entry.entryId,
        name: scalar.name,
        units: scalar.units,
        payload: {
          entryId: entry.entryId,
          path: scalar.path,
          name: scalar.name,
          friendlyPath: scalar.friendlyPath,
          units: scalar.units,
          value: scalar.value,
        },
      });
    }
  }
  return dedupeTargetKeys(clusterByName(items)
    .filter((cluster) => cluster.entryIds.size >= 2)
    .map((cluster) => ({
      kind: 'scalar' as const,
      key: `scalar:${normalizeName(cluster.canonicalName)}`,
      displayName: cluster.canonicalName,
      confidence: cluster.confidence,
      unitsWarning: cluster.unitsWarning,
      bindings: cluster.items.map((item) => item.payload),
      coverage: cluster.entryIds.size,
      total: entries.length,
    })));
}

/**
 * Groups numeric columns across entries into candidate targets.
 * Only tables with a user-defined index participate; the index and split columns
 * themselves are not candidates. Both maps: entryId -> (tablePath -> column name).
 */
export function matchColumnTargets(
  entries: ComparisonEntryNodes[],
  indexColumns: Map<string, Map<string, string>>,
  splitColumns?: Map<string, Map<string, string>>,
): ColumnTarget[] {
  const items: ClusterItem<ColumnBinding>[] = [];
  for (const entry of entries) {
    for (const table of entry.tables) {
      const indexColumnName = indexColumns.get(entry.entryId)?.get(table.path);
      if (!indexColumnName)
        continue;
      const splitColumnName = splitColumns?.get(entry.entryId)?.get(table.path);
      for (const column of table.columns) {
        if (column.name === indexColumnName || column.name === splitColumnName || !isNumericType(column.type))
          continue;
        items.push({
          entryId: entry.entryId,
          name: column.name,
          units: column.units,
          secondaryName: table.name,
          payload: {
            entryId: entry.entryId,
            tablePath: table.path,
            tableName: table.name,
            tableFriendlyPath: table.friendlyPath,
            columnName: column.name,
            units: column.units,
            indexColumnName,
            splitColumnName,
          },
        });
      }
    }
  }
  return dedupeTargetKeys(clusterByName(items)
    .filter((cluster) => cluster.entryIds.size >= 2)
    .map((cluster) => ({
      kind: 'column' as const,
      key: `column:${normalizeName(cluster.canonicalName)}:${normalizeName(cluster.canonicalSecondary ?? '')}`,
      displayName: cluster.canonicalName,
      confidence: cluster.confidence,
      unitsWarning: cluster.unitsWarning,
      bindings: cluster.items.map((item) => item.payload),
      coverage: cluster.entryIds.size,
      total: entries.length,
    })));
}


export type ExclusionReason =
  | 'no similar data'
  | 'index not set';

export interface EntryStatus {
  entryId: string;
  matched: boolean;
  reason?: ExclusionReason;
}

/** Per-entry participation status for a selected target. */
export function getEntryStatuses(
  entries: ComparisonEntryNodes[],
  target: ComparisonTarget | null,
  indexColumns: Map<string, Map<string, string>>,
): EntryStatus[] {
  return entries.map((entry) => {
    if (!target)
      return {entryId: entry.entryId, matched: false, reason: 'no similar data'};
    const binding = (target.bindings as {entryId: string}[]).find((b) => b.entryId === entry.entryId);
    if (binding)
      return {entryId: entry.entryId, matched: true};
    if (target.kind === 'column' &&
      entry.tables.length > 0 &&
      entry.tables.every((table) => !indexColumns.get(entry.entryId)?.get(table.path)))
      return {entryId: entry.entryId, matched: false, reason: 'index not set'};
    return {entryId: entry.entryId, matched: false, reason: 'no similar data'};
  });
}
