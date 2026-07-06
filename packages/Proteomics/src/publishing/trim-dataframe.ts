import * as DG from 'datagrok-api/dg';

import {SEMTYPE} from '../utils/proteomics-types';
import {findColumn, findProteomicsColumns} from '../utils/column-detection';
import {META_COLUMNS, PUBLISHED_TAGS, PublishedMetadata, setPublishedTags} from './publish-state';

/**
 * Deep-clone + allowlist primitives that produce frozen DataFrame snapshots
 * ready for `DG.Project.save`. Pitfall 1 (stale-snapshot leak) and Pitfall 3
 * (tag-stripping serializer) are mitigated HERE, not in `publishAnalysis`.
 *
 * The orchestrator (Plan 04) never touches the source DataFrame directly — it
 * hands the source to {@link trimForPublish} and operates on the returned
 * clone. This is the only file in the package where `df.clone(...)` lives for
 * publishing.
 *
 * Spike 15-00 confirmed that `proteomics.*` tags and `Proteomics-*` semTypes
 * survive `DG.Project.save → find → open` on `release/1.27.3`, but does NOT
 * confirm they survive `df.clone(...)`. The defensive re-set is no-cost either
 * way — re-setting an already-present tag is idempotent.
 */

interface MetadataColumnSpec {
  name: string;
  value: string | number | Date | boolean | null;
  /** 'string' | 'datetime' | 'float' | 'int' — drives `columns.addNew*` */
  type: 'string' | 'datetime' | 'float' | 'int';
  /** When true, write '' (string) or 0 (numeric) instead of skipping null. */
  emptyForNull?: boolean;
}

function addMetadataColumn(df: DG.DataFrame, spec: MetadataColumnSpec): void {
  if (df.columns.contains(spec.name)) df.columns.remove(spec.name);
  let col: DG.Column;
  switch (spec.type) {
    case 'datetime':
      col = df.columns.addNewDateTime(spec.name);
      break;
    case 'float':
      col = df.columns.addNewFloat(spec.name);
      break;
    case 'int':
      col = df.columns.addNewInt(spec.name);
      break;
    case 'string':
    default:
      col = df.columns.addNewString(spec.name);
      break;
  }
  let constant: any;
  if (spec.value == null) {
    if (!spec.emptyForNull) return;
    constant = spec.type === 'string' ? '' : 0;
  } else if (spec.value instanceof Date) {
    constant = spec.value;
  } else if (typeof spec.value === 'boolean') {
    constant = spec.value ? 'true' : 'false';
  } else if (spec.type === 'string') {
    constant = String(spec.value);
  } else {
    constant = spec.value;
  }
  col.init(() => constant);
  col.setTag('.hidden', 'true');
}

/**
 * Returns a deep clone of `source` with only the protein-level allowlist
 * (Protein ID, Gene, log2FC, p-value, adj.p-value, significant, [direction])
 * + 13 belt-and-braces `_meta_*` columns + every required
 * `proteomics.published*` tag re-set explicitly. Source DataFrame is NOT
 * mutated. Throws on missing required columns (defensive — the Plan 07 menu
 * handler should have already gated via `requireDifferentialExpression`).
 */
export function trimForPublish(source: DG.DataFrame, meta: PublishedMetadata): DG.DataFrame {
  const cols = findProteomicsColumns(source);
  const proteinId = cols.proteinId;
  const geneName = cols.geneName;
  const log2fc = cols.log2fc;
  // 'adj.p-value' and 'p-value' both substring-match 'p-value'; resolve adj first
  // and then exclude its name when looking for the raw p-value column so the
  // allowlist doesn't end up with duplicates.
  const adjPValue = findColumn(source, '', ['adj.p-value', 'adj.p', 'fdr', 'q-value']);
  const pValue = source.columns.toList().find((c) => {
    const n = c.name.toLowerCase();
    return (n === 'p-value' || n === 'pvalue') && c.name !== (adjPValue?.name ?? '');
  }) ?? null;
  const significant = findColumn(source, '', ['significant', 'sig']);
  const direction = findColumn(source, '', ['direction', 'regulation', 'up_down']);

  const required: Array<[string, DG.Column | null]> = [
    ['Protein ID', proteinId],
    ['log2FC', log2fc],
    ['p-value', pValue],
    ['adj.p-value', adjPValue],
    ['significant', significant],
  ];
  for (const [label, col] of required) {
    if (col == null)
      throw new Error(`Cannot publish: missing required column [${label}]. Run Differential Expression to completion before publishing.`);
  }

  const allowlist: string[] = Array.from(new Set([
    proteinId!.name,
    geneName?.name,
    log2fc!.name,
    pValue!.name,
    adjPValue!.name,
    significant!.name,
    direction?.name,
  ].filter((n): n is string => typeof n === 'string' && n.length > 0)));

  const frozen = source.clone(null, allowlist);

  const carryForwardTags = [
    'proteomics.source',
    'proteomics.de_method',
    'proteomics.groups',
    'proteomics.de_complete',
  ];
  for (const k of carryForwardTags) {
    const v = source.getTag(k);
    if (v != null && v !== '') frozen.setTag(k, v);
  }

  setPublishedTags(frozen, meta);

  const dateSlice = meta.publishedAt instanceof Date
    ? meta.publishedAt.toISOString().slice(0, 10)
    : new Date(String(meta.publishedAt)).toISOString().slice(0, 10);
  frozen.name = `${source.name || 'analysis'}_published_${dateSlice}`;

  const semTypeAssignments: Array<[DG.Column | null, string]> = [
    [proteinId, SEMTYPE.PROTEIN_ID],
    [geneName, SEMTYPE.GENE_SYMBOL],
    [log2fc, SEMTYPE.LOG2FC],
    [pValue, SEMTYPE.P_VALUE],
  ];
  for (const [src, sem] of semTypeAssignments) {
    if (src == null) continue;
    const c = frozen.col(src.name);
    if (c != null) c.semType = sem;
  }

  const includesEnrichmentStr = meta.includesEnrichment ? 'true' : 'false';
  // Numeric metadata (FC / p threshold, version) stored as STRINGS to avoid
  // Float32 precision loss in `addNewFloat` (which is a single-precision
  // column type). The threshold values are exact doubles in the source meta;
  // string storage preserves them exactly, and `getPublishedMetadata` parses
  // them back via parseFloat / parseInt.
  const specs: MetadataColumnSpec[] = [
    {name: META_COLUMNS.PUBLISHED, value: 'true', type: 'string'},
    {name: META_COLUMNS.PUBLISHED_AT, value: meta.publishedAt, type: 'datetime'},
    {name: META_COLUMNS.PUBLISHED_BY, value: meta.publishedBy, type: 'string'},
    {name: META_COLUMNS.PUBLISHED_BY_EMAIL, value: meta.publishedByEmail ?? '', type: 'string', emptyForNull: true},
    {name: META_COLUMNS.PUBLISHED_TARGET, value: meta.target, type: 'string'},
    {name: META_COLUMNS.PUBLISHED_DE_METHOD, value: meta.deMethod, type: 'string'},
    {name: META_COLUMNS.PUBLISHED_FC_THRESHOLD, value: String(meta.fcThreshold), type: 'string'},
    {name: META_COLUMNS.PUBLISHED_P_THRESHOLD, value: String(meta.pThreshold), type: 'string'},
    {name: META_COLUMNS.PUBLISHED_VERSION, value: String(meta.version), type: 'string'},
    {name: META_COLUMNS.PUBLISHED_ID, value: meta.publishId, type: 'string'},
    {name: META_COLUMNS.PUBLISHED_INCLUDES_ENRICHMENT, value: includesEnrichmentStr, type: 'string'},
    {name: META_COLUMNS.SUPERSEDES, value: meta.supersedes ?? '', type: 'string', emptyForNull: true},
    {name: META_COLUMNS.SUPERSEDED_BY, value: meta.supersededBy ?? '', type: 'string', emptyForNull: true},
  ];
  for (const spec of specs)
    addMetadataColumn(frozen, spec);

  return frozen;
}

/**
 * Returns a deep clone of an enrichment DataFrame containing only the term /
 * source / p / adj.p / intersection allowlist (+ optional direction column
 * from the Phase 13 split). Called by the Plan 04 orchestrator opportunistically
 * (D-05): only when the source protein DataFrame carries
 * `proteomics.enrichment === 'true'` AND an enrichment DataFrame is present
 * in `grok.shell.tables`. If neither precondition holds, the orchestrator
 * skips this call entirely and the published Project ships protein-only.
 *
 * Throws on missing required columns. Direction column is soft: omitted from
 * the allowlist when absent rather than failing.
 */
export function trimEnrichmentForPublish(enrichSource: DG.DataFrame, meta: PublishedMetadata): DG.DataFrame {
  const term = findColumn(enrichSource, '', ['term name', 'term', 'pathway']);
  const sourceCol = findColumn(enrichSource, '', ['source']);
  // The g:Profiler producer (buildEnrichmentDf, enrichment.ts) emits ONE
  // FDR-corrected significance column named 'FDR' — there is no separate raw
  // p-value. Accept any single significance column (FDR / adj.p-value /
  // p-value) rather than requiring a p-value + adj.p-value pair that real
  // enrichment tables never carry.
  const sigCol = findColumn(enrichSource, '', ['fdr', 'adj.p-value', 'adj.p', 'q-value', 'p-value', 'pvalue']);
  const intersection = findColumn(enrichSource, '', ['intersection', 'genes']);
  const directionCol = findColumn(enrichSource, '', ['direction', 'regulation']);
  // Soft (omitted when absent): the published dot plot needs these to render
  // (x = Gene Ratio, size = Gene Count). Older enrichment frames may lack them.
  const geneRatioCol = findColumn(enrichSource, '', ['gene ratio']);
  const geneCountCol = findColumn(enrichSource, '', ['gene count']);

  const required: Array<[string, DG.Column | null]> = [
    ['Term Name', term],
    ['Source', sourceCol],
    ['significance (FDR / adj.p-value)', sigCol],
    ['Intersection', intersection],
  ];
  for (const [label, col] of required) {
    if (col == null)
      throw new Error(`Cannot publish enrichment: missing required column [${label}]. Re-run enrichment before publishing.`);
  }

  const allowlist: string[] = [
    term!.name,
    sourceCol!.name,
    sigCol!.name,
    intersection!.name,
    directionCol?.name,
    geneRatioCol?.name,
    geneCountCol?.name,
  ].filter((n): n is string => typeof n === 'string' && n.length > 0);

  const enrichClone = enrichSource.clone(null, allowlist);
  enrichClone.setTag('proteomics.enrichment', 'true');

  const dateSlice = meta.publishedAt instanceof Date
    ? meta.publishedAt.toISOString().slice(0, 10)
    : new Date(String(meta.publishedAt)).toISOString().slice(0, 10);
  enrichClone.name = `enrichment_published_${dateSlice}`;

  addMetadataColumn(enrichClone, {
    name: META_COLUMNS.PUBLISHED_ID,
    value: meta.publishId,
    type: 'string',
  });
  addMetadataColumn(enrichClone, {
    name: META_COLUMNS.PUBLISHED_INCLUDES_ENRICHMENT,
    value: 'true',
    type: 'string',
  });

  // Belt-and-braces tag so a reviewer reopening the enrichment DataFrame
  // alone can recover its provenance.
  enrichClone.setTag(PUBLISHED_TAGS.PUBLISHED, 'true');
  enrichClone.setTag(PUBLISHED_TAGS.PUBLISHED_ID, meta.publishId);
  enrichClone.setTag(PUBLISHED_TAGS.PUBLISHED_INCLUDES_ENRICHMENT, 'true');

  return enrichClone;
}
