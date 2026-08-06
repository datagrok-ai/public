// Pure matching engine — groups scalars and columns across runs into comparison targets
// by name similarity and units compatibility. No DG/platform dependencies, unit-testable.
// User documentation lives in help/compute/run-comparison.md

import {
  ComparisonEntryNodes, MatchConfidence, ScalarBinding, ColumnBinding, ColumnCandidate,
  CandidateOverrides, ScalarTarget, ColumnTarget, TargetBase, isNumericType, candidateId,
} from './types';

export type UnitsCompatibility = 'match' | 'warn' | 'mismatch';

export const FUZZY_NAME_THRESHOLD = 0.7;

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

export interface TableMatchKey {
  indexColumnName: string;
  splitColumnName?: string;
}

// tables are comparable when their index columns name-match and their split columns
// are either both absent or name-match: a split table charts per-category series,
// so pairing it with an unsplit one would be misleading
export function tablesCompatible(a: TableMatchKey, b: TableMatchKey): boolean {
  if (!nameMatchConfidence(a.indexColumnName, b.indexColumnName))
    return false;
  if (!a.splitColumnName && !b.splitColumnName)
    return true;
  if (!a.splitColumnName || !b.splitColumnName)
    return false;
  return nameMatchConfidence(a.splitColumnName, b.splitColumnName) != null;
}

interface ClusterItem<P> {
  entryId: string;
  name: string;
  units?: string;
  // used only to break ties between equally-confident candidates (e.g. table name for columns)
  secondaryName?: string;
  tableKey?: TableMatchKey;
  raw?: boolean;
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
      const seedKey = cluster.items[0].tableKey;
      if (seedKey && item.tableKey && !tablesCompatible(seedKey, item.tableKey))
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
 *
 * Each target carries the full list of compatible candidates: greedy-clustered items are
 * enabled (auto), other compatible items are attached disabled, and raw (standalone) items
 * are enabled in every cluster they fit. Overrides (by target key and candidate id) flip
 * individual candidates; derived fields (bindings, coverage, confidence, unitsWarning)
 * reflect the enabled subset only.
 */
export function matchColumnTargets(
  entries: ComparisonEntryNodes[],
  indexColumns: Map<string, Map<string, string>>,
  splitColumns?: Map<string, Map<string, string>>,
  overrides?: CandidateOverrides,
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
          tableKey: {indexColumnName, splitColumnName},
          raw: entry.sourceKind === 'raw',
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

  const withCandidates = clusterByName(items).map((cluster) => {
    const members = new Set(cluster.items);
    const seed = cluster.items[0];
    const candidates: ColumnCandidate[] = [];
    for (const item of items) {
      const confidence = nameMatchConfidence(cluster.canonicalName, item.name);
      if (!confidence)
        continue;
      const units = unitsCompatibility(seed.units, item.units);
      if (units === 'mismatch')
        continue;
      if (seed.tableKey && item.tableKey && !tablesCompatible(seed.tableKey, item.tableKey))
        continue;
      const auto = members.has(item);
      // raw items join every cluster whose name they share (up to normalization) — but a
      // fuzzy guess is not confident enough to become the run's pick by default
      candidates.push({
        binding: item.payload,
        confidence,
        unitsWarn: units === 'warn',
        auto,
        enabled: item.raw ? confidence !== 'fuzzy' : auto,
      });
    }
    return {cluster, candidates};
  });

  // survival counts default enablement, so user toggles can never remove a target
  const targets = dedupeTargetKeys(withCandidates
    .filter(({candidates}) =>
      new Set(candidates.filter((c) => c.enabled).map((c) => c.binding.entryId)).size >= 2)
    .map(({cluster, candidates}) => ({
      kind: 'column' as const,
      key: `column:${normalizeName(cluster.canonicalName)}:${normalizeName(cluster.canonicalSecondary ?? '')}`,
      displayName: cluster.canonicalName,
      confidence: cluster.confidence,
      unitsWarning: cluster.unitsWarning,
      candidates,
      bindings: [] as ColumnBinding[],
      coverage: 0,
      total: entries.length,
    })));

  // overrides are keyed by the final (deduped) target key, so they resolve only here.
  // radio semantics: at most one enabled candidate per run — an explicit user pick
  // beats a default one, ties resolve to the first in candidate order
  for (const target of targets) {
    const targetOverrides = overrides?.[target.key];
    const chosen = new Map<string, number>();
    target.candidates.forEach((candidate, index) => {
      const override = targetOverrides?.[candidateId(candidate.binding)];
      if (!(override ?? candidate.enabled))
        return;
      const entryId = candidate.binding.entryId;
      const current = chosen.get(entryId);
      if (current == null) {
        chosen.set(entryId, index);
        return;
      }
      const currentExplicit =
        targetOverrides?.[candidateId(target.candidates[current].binding)] === true;
      if (override === true && !currentExplicit)
        chosen.set(entryId, index);
    });
    target.candidates = target.candidates.map((candidate, index) => ({
      ...candidate,
      enabled: chosen.get(candidate.binding.entryId) === index,
    }));
    const enabled = target.candidates.filter((candidate) => candidate.enabled);
    target.bindings = enabled.map((candidate) => candidate.binding);
    target.coverage = new Set(target.bindings.map((b) => b.entryId)).size;
    target.confidence = enabled.reduce<MatchConfidence>(
      (acc, candidate) => weakerConfidence(acc, candidate.confidence), 'exact');
    target.unitsWarning = enabled.some((candidate) => candidate.unitsWarn);
  }
  return targets;
}
