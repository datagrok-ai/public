// Pure matching engine — groups scalars and columns across runs into comparison targets
// by name similarity and units compatibility. No DG/platform dependencies, unit-testable.
// User documentation lives in help/compute/run-comparison.md

import {
  ComparisonEntryNodes, MatchConfidence, ScalarBinding, ColumnBinding, ColumnCandidate,
  CandidateOverrides, NameMapping, UnitMapping, ScalarTarget, ColumnTarget, TargetBase, isNumericType,
  candidateId, TIME_UNIT_ALIASES,
} from './types';

export type UnitsCompatibility = 'match' | 'warn' | 'convert' | 'mismatch';

// multiplier from one unit into another, or null when no mapping covers the pair
export type UnitConverter = (from?: string, to?: string) => number | null;

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

/** Normalized name -> group id for user-defined mappings; pairs merge transitively. */
export function buildAliasGroups(mappings: NameMapping[]): Map<string, number> {
  const groups = new Map<string, number>();
  let nextId = 0;
  for (const {from, to} of mappings) {
    const a = normalizeName(from);
    const b = normalizeName(to);
    if (!a || !b)
      continue;
    const ga = groups.get(a);
    const gb = groups.get(b);
    if (ga == null && gb == null) {
      groups.set(a, nextId);
      groups.set(b, nextId);
      nextId++;
    } else if (ga == null) {
      groups.set(a, gb!);
    } else if (gb == null) {
      groups.set(b, ga);
    } else if (ga !== gb) {
      for (const [name, id] of groups) {
        if (id === gb)
          groups.set(name, ga);
      }
    }
  }
  return groups;
}

export function nameMatchConfidence(
  a: string, b: string, aliases?: Map<string, number>,
): MatchConfidence | null {
  if (a === b)
    return 'exact';
  if (normalizeName(a) === normalizeName(b))
    return 'normalized';
  if (aliases) {
    const group = aliases.get(normalizeName(a));
    if (group != null && group === aliases.get(normalizeName(b)))
      return 'normalized';
  }
  if (nameSimilarity(a, b) >= FUZZY_NAME_THRESHOLD)
    return 'fuzzy';
  return null;
}

const CONFIDENCE_RANK: Record<MatchConfidence, number> = {exact: 0, normalized: 1, fuzzy: 2};

export function weakerConfidence(a: MatchConfidence, b: MatchConfidence): MatchConfidence {
  return CONFIDENCE_RANK[a] >= CONFIDENCE_RANK[b] ? a : b;
}

const SUPERSCRIPTS: Record<string, string> = {
  '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9', '⁻': '-',
};

/**
 * Case-sensitive canonical units (MPa and mPa differ by 10^9): only spelling variants that
 * never change the meaning fold together — spaces, a trailing period, micro signs, superscript
 * digits, liters as l or L, and time unit aliases (sec, hours, ...). Never fuzzy.
 */
export function normalizeUnits(units?: string): string {
  const compact = (units ?? '').replace(/\s+/g, '').replace(/\.$/, '')
    .replace(/[µμ]/g, 'u')
    .replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]/g, (c) => SUPERSCRIPTS[c])
    .replace(/(^|[^A-Za-z])([munpdck]?)l(?![A-Za-z])/g, '$1$2L');
  // whole-word aliases may differ in case (Hours, SEC); one- and two-letter symbols may not
  const alias = TIME_UNIT_ALIASES.get(compact) ??
    (compact.length > 2 ? TIME_UNIT_ALIASES.get(compact.toLowerCase()) : undefined);
  return alias ?? compact;
}

/** Both units set and different after spelling normalization, with a finite positive factor. */
export function isValidUnitMapping({from, factor, to}: UnitMapping): boolean {
  const [a, b] = [normalizeUnits(from), normalizeUnits(to)];
  return !!a && !!b && a !== b && Number.isFinite(factor) && factor > 0;
}

/**
 * Converter over user unit mappings: each pair works in both directions, pairs never chain,
 * invalid mappings are ignored.
 */
export function buildUnitConverter(mappings: UnitMapping[]): UnitConverter {
  const factors = new Map<string, number>();
  for (const {from, factor, to} of mappings.filter(isValidUnitMapping)) {
    const [a, b] = [normalizeUnits(from), normalizeUnits(to)];
    factors.set(`${a}\u0000${b}`, factor);
    factors.set(`${b}\u0000${a}`, 1 / factor);
  }
  return (from, to) => factors.get(`${normalizeUnits(from)}\u0000${normalizeUnits(to)}`) ?? null;
}

export function unitsCompatibility(a?: string, b?: string, convert?: UnitConverter): UnitsCompatibility {
  const ua = normalizeUnits(a);
  const ub = normalizeUnits(b);
  if (ua === ub)
    return 'match';
  if (!ua || !ub)
    return 'warn';
  return convert?.(ua, ub) != null ? 'convert' : 'mismatch';
}

// multiplier from `units` into `targetUnits`, omitted when no conversion is involved
const conversionFactor = (units?: string, targetUnits?: string, convert?: UnitConverter) =>
  unitsCompatibility(units, targetUnits, convert) === 'convert' ? convert!(units, targetUnits)! : undefined;

interface ClusterItem<P> {
  entryId: string;
  name: string;
  units?: string;
  // used only to break ties between equally-confident candidates (e.g. table name for columns)
  secondaryName?: string;
  indexColumnName?: string;
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

function seedMatch<P>(
  cluster: Cluster<P>, item: ClusterItem<P>, aliases?: Map<string, number>, convert?: UnitConverter,
): {confidence: MatchConfidence, units: UnitsCompatibility} | null {
  const seed = cluster.items[0];
  const confidence = nameMatchConfidence(cluster.canonicalName, item.name, aliases);
  if (!confidence)
    return null;
  const units = unitsCompatibility(seed.units, item.units, convert);
  if (units === 'mismatch')
    return null;
  if (seed.indexColumnName && item.indexColumnName &&
    !nameMatchConfidence(seed.indexColumnName, item.indexColumnName, aliases))
    return null;
  return {confidence, units};
}

function bestCluster<P>(
  clusters: Cluster<P>[], item: ClusterItem<P>, allowFuzzy: boolean, aliases?: Map<string, number>,
  convert?: UnitConverter,
): {cluster: Cluster<P>, confidence: MatchConfidence} | null {
  let best: {cluster: Cluster<P>, confidence: MatchConfidence, score: number} | null = null;
  for (const cluster of clusters) {
    if (cluster.entryIds.has(item.entryId))
      continue;
    const confidence = seedMatch(cluster, item, aliases, convert)?.confidence;
    if (!confidence || (confidence === 'fuzzy' && !allowFuzzy))
      continue;
    const score = nameSimilarity(cluster.canonicalName, item.name) +
      (cluster.canonicalSecondary && item.secondaryName ?
        nameSimilarity(cluster.canonicalSecondary, item.secondaryName) * 0.5 : 0);
    if (!best || CONFIDENCE_RANK[confidence] < CONFIDENCE_RANK[best.confidence] ||
      (confidence === best.confidence && score > best.score))
      best = {cluster, confidence, score};
  }
  return best;
}

function addToCluster<P>(cluster: Cluster<P>, item: ClusterItem<P>, confidence: MatchConfidence) {
  cluster.items.push(item);
  cluster.entryIds.add(item.entryId);
  cluster.confidence = weakerConfidence(cluster.confidence, confidence);
  if (unitsCompatibility(cluster.items[0].units, item.units) === 'warn')
    cluster.unitsWarning = true;
}

function clusterByName<P>(
  items: ClusterItem<P>[], aliases?: Map<string, number>, convert?: UnitConverter,
): Cluster<P>[] {
  const clusters: Cluster<P>[] = [];
  for (const item of items) {
    const best = bestCluster(clusters, item, false, aliases, convert);
    if (best) {
      addToCluster(best.cluster, item, best.confidence);
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
  // fuzzy matches run second, so a fuzzy name never takes a run's slot from an exact one;
  // lone items join earlier clusters or established ones, keeping the first-seen name canonical
  for (const lone of clusters.filter((cluster) => cluster.items.length === 1)) {
    const position = clusters.indexOf(lone);
    const candidates = clusters.filter((cluster, i) =>
      i < position || (i > position && cluster.items.length > 1));
    const best = bestCluster(candidates, lone.items[0], true, aliases, convert);
    if (best) {
      addToCluster(best.cluster, lone.items[0], best.confidence);
      clusters.splice(clusters.indexOf(lone), 1);
    }
  }
  return clusters;
}

// distinct clusters can share a canonical name (a cluster takes at most one item per entry),
// and duplicate keys corrupt keyed list rendering
function dedupeTargetKeys<T extends TargetBase>(targets: T[]): T[] {
  const seen = new Map<string, number>();
  return targets.map((target) => {
    const count = seen.get(target.key) ?? 0;
    seen.set(target.key, count + 1);
    return count > 0 ? {...target, key: `${target.key}:${count + 1}`} : target;
  });
}

/** Groups numeric scalars across entries into candidate targets (coverage >= 2). */
export function matchScalarTargets(
  entries: ComparisonEntryNodes[], aliases?: Map<string, number>, convert?: UnitConverter,
): ScalarTarget[] {
  const items: ClusterItem<ScalarBinding>[] = [];
  for (const entry of entries) {
    for (const scalar of entry.scalars) {
      if (!isNumericType(scalar.valueType))
        continue;
      items.push({
        entryId: entry.entryId,
        name: scalar.name,
        units: scalar.units,
        // pairs repeated workflow steps by position (Day 2 with Day 2)
        secondaryName: scalar.friendlyPath,
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
  return dedupeTargetKeys(clusterByName(items, aliases, convert)
    .filter((cluster) => cluster.entryIds.size >= 2)
    .map((cluster) => ({
      kind: 'scalar' as const,
      key: `scalar:${normalizeName(cluster.canonicalName)}`,
      displayName: cluster.canonicalName,
      units: cluster.items[0].units,
      confidence: cluster.confidence,
      unitsWarning: cluster.unitsWarning,
      bindings: cluster.items.map((item) => {
        const factor = conversionFactor(item.units, cluster.items[0].units, convert);
        return factor == null ? item.payload : {...item.payload, factor};
      }),
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
  aliases?: Map<string, number>,
  convert?: UnitConverter,
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
          secondaryName: table.friendlyPath ?? table.name,
          indexColumnName,
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

  const withCandidates = clusterByName(items, aliases, convert).map((cluster) => {
    const members = new Set(cluster.items);
    const seed = cluster.items[0];
    const candidates: ColumnCandidate[] = [];
    for (const item of items) {
      const match = seedMatch(cluster, item, aliases, convert);
      if (!match)
        continue;
      const {confidence, units} = match;
      const factor = conversionFactor(item.units, seed.units, convert);
      const auto = members.has(item);
      // raw items join every cluster whose name they share (up to normalization) — but a
      // fuzzy guess is not confident enough to become the run's pick by default
      candidates.push({
        binding: {
          ...item.payload,
          ...factor != null ? {factor} : {},
        },
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
      units: cluster.items[0].units,
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
