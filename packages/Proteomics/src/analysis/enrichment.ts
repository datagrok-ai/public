import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {parseAccession} from '../panels/uniprot-panel';
import {findProteomicsColumns} from '../utils/column-detection';
import {SEMTYPE, DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD} from '../utils/proteomics-types';
import {openEnrichmentVisualization} from '../viewers/enrichment-viewers';
import {requireDifferentialExpression} from './differential-expression';
import {getOrganism, setOrganism} from './experiment-setup';
import {ORGANISM_LIST, detectOrganismCode, organismDisplayForCode} from '../utils/organisms';

// ORGANISM_LIST now lives in utils/organisms (single source); re-exported here so
// existing importers (and the enrichment tests) keep resolving it from this module.
export {ORGANISM_LIST};

// --- Constants ---

const GPROFILER_BASE = 'https://biit.cs.ut.ee/gprofiler';

// --- Types ---

export interface ConvertResult {
  incoming: string;
  converted: string;
  name: string;
  description: string;
  namespaces: string;
  n_incoming: number;
  n_converted: number;
}

export interface GostResult {
  native: string;
  name: string;
  source: string;
  p_value: number;
  significant: boolean;
  term_size: number;
  query_size: number;
  intersection_size: number;
  effective_domain_size: number;
  precision: number;
  recall: number;
  intersections: string[][];
}

// --- API Client ---

async function fetchWithTimeout(
  url: string, options: RequestInit, timeoutMs: number = 30000,
): Promise<Response> {
  // External URLs go through grok.dapi.fetchProxy to avoid browser CORS.
  // Datagrok's proxy does not currently honor AbortSignal, so we race the
  // request against a timeout promise to enforce the bound.
  const proxied = grok.dapi.fetchProxy(url, options);
  const timeout = new Promise<never>((_, reject) =>
    setTimeout(() => reject(new Error('g:Profiler API request timed out. The service may be temporarily unavailable.')), timeoutMs),
  );
  const resp = await Promise.race([proxied, timeout]);
  if (!resp.ok)
    throw new Error(`g:Profiler API returned status ${resp.status}`);
  return resp;
}

export async function gConvert(
  accessions: string[],
  organism: string = 'hsapiens',
): Promise<ConvertResult[]> {
  const resp = await fetchWithTimeout(`${GPROFILER_BASE}/api/convert/convert/`, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({
      organism: organism,
      query: accessions,
      target: 'ENSG',
      numeric_ns: '',
      output: 'json',
    }),
  });
  const data = await resp.json();
  return data.result ?? [];
}

export async function gGOSt(
  queryGenes: string[],
  backgroundGenes: string[],
  organism: string,
  sources: string[],
  threshold: number = DEFAULT_P_THRESHOLD,
): Promise<GostResult[]> {
  const resp = await fetchWithTimeout(`${GPROFILER_BASE}/api/gost/profile/`, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({
      organism: organism,
      query: queryGenes,
      sources: sources,
      user_threshold: threshold,
      significance_threshold_method: 'fdr',
      domain_scope: 'custom',
      background: backgroundGenes,
      all_results: true,
      ordered: false,
      no_evidences: false,
      combined: false,
      measure_underrepresentation: false,
      no_iea: false,
      numeric_ns: '',
      output: 'json',
    }),
  });
  const data = await resp.json();
  // g:GOSt nests results: data.result is an array of per-query result objects,
  // each containing a `.result` array of term results.
  if (!data.result)
    return [];
  // Handle array of query results (standard structure)
  if (Array.isArray(data.result) && data.result.length > 0 && data.result[0]?.result)
    return data.result[0].result as GostResult[];
  // Fallback: if data.result is directly the results array
  if (Array.isArray(data.result))
    return data.result as GostResult[];
  return [];
}

// --- Smart Pathway Filter (R5/D-13/D-14) ---

// LOCKED CLIENT CONTRACT — port verbatim from
// ~/Downloads/ck/CKomics_tool2.py:4685-4736 (apply_smart_pathway_filtering).
// Do not re-derive heuristics; any change requires a CONTEXT update.

export const GENERIC_PARENT_TERMS = [
  'localization', 'cellular component organization', 'transport',
  'cellular process', 'biological process', 'metabolic process',
] as const;

export const SPECIFIC_CHILD_TERMS = ['actin', 'vesicle', 'endocytosis', 'cytoskeleton'] as const;

export interface SmartFilterStats {
  total: number;
  kept: number;
  droppedParents: number;
  cappedAtN: number;
}

/** Verbatim port of CK-omics apply_smart_pathway_filtering (lines 4685-4736).
 *
 * Algorithm:
 *  1. Partition into GO:BP rows and `other_data` (everything else).
 *  2. Sort GO:BP by p_value ASCENDING.
 *  3. Walk GO:BP in sorted order: drop a row when its name contains any
 *     GENERIC_PARENT term AND the already-kept list already contains a
 *     SPECIFIC_CHILD-term row. Otherwise keep it. Stop at `maxPerSource`.
 *  4. Non-GO:BP rows are sorted by p_value ASC and `.head(maxPerSource)` —
 *     a COMBINED cap across all non-GO:BP sources, not per-source
 *     (Assumption A4; CK-omics line 4731 explicitly uses combined .head). */
export function applySmartPathwayFilter(
  results: GostResult[], maxPerSource = 15,
): {kept: GostResult[]; stats: SmartFilterStats} {
  if (results.length === 0) {
    return {kept: [], stats: {total: 0, kept: 0, droppedParents: 0, cappedAtN: maxPerSource}};
  }

  const goBp = results.filter((r) => r.source === 'GO:BP');
  const other = results.filter((r) => r.source !== 'GO:BP');

  let droppedParents = 0;
  const filteredGoBp: GostResult[] = [];
  if (goBp.length > 0) {
    const sorted = [...goBp].sort((a, b) => a.p_value - b.p_value);
    for (const r of sorted) {
      const name = (r.name ?? '').toLowerCase();
      const isGeneric = GENERIC_PARENT_TERMS.some((g) => name.includes(g));
      if (isGeneric) {
        const hasSpecificChild = filteredGoBp.some((existing) =>
          SPECIFIC_CHILD_TERMS.some((s) => (existing.name ?? '').toLowerCase().includes(s)));
        if (hasSpecificChild) { droppedParents++; continue; }
      }
      filteredGoBp.push(r);
      if (filteredGoBp.length >= maxPerSource) break;
    }
  }

  const otherSorted = [...other].sort((a, b) => a.p_value - b.p_value);
  const otherCapped = otherSorted.slice(0, maxPerSource);

  const kept = [...filteredGoBp, ...otherCapped].sort((a, b) => a.p_value - b.p_value);
  return {
    kept,
    stats: {
      total: results.length,
      kept: kept.length,
      droppedParents,
      cappedAtN: maxPerSource,
    },
  };
}

// --- Result Builder ---

export function buildEnrichmentDf(
  results: GostResult[],
  queryGenes: string[],
  pThreshold: number = DEFAULT_P_THRESHOLD,
  direction?: 'Up' | 'Down',
): DG.DataFrame {
  const n = results.length;
  const df = DG.DataFrame.create(n);
  df.name = 'Enrichment Results';

  // Pre-derive the comma-joined member-gene string per row so the column init
  // doesn't repeat the intersections-array walk per cell.
  const memberGeneStrs = results.map((r) => {
    if (!r.intersections) return '';
    const members: string[] = [];
    const lim = Math.min(queryGenes.length, r.intersections.length);
    for (let j = 0; j < lim; j++) {
      if (r.intersections[j] && r.intersections[j].length > 0)
        members.push(queryGenes[j]);
    }
    return members.join(', ');
  });

  const sourceCol = df.columns.addNewString('Source');
  const termIdCol = df.columns.addNewString('Term ID');
  const termNameCol = df.columns.addNewString('Term Name');
  // g:GOSt with significance_threshold_method='fdr' returns the FDR-corrected
  // p-value (and no separate raw p-value) — we expose it under the FDR column
  // to keep semantics honest. Downstream consumers should read FDR.
  const fdrCol = df.columns.addNewFloat('FDR');
  const geneCountCol = df.columns.addNewInt('Gene Count');
  const geneRatioCol = df.columns.addNewFloat('Gene Ratio');
  const intersectionCol = df.columns.addNewString('Intersection');

  sourceCol.init((i) => results[i].source);
  termIdCol.init((i) => results[i].native);
  termNameCol.init((i) => results[i].name);
  fdrCol.init((i) => results[i].p_value);
  geneCountCol.init((i) => results[i].intersection_size);
  geneRatioCol.init((i) => results[i].precision);
  intersectionCol.init((i) => memberGeneStrs[i]);

  const fdrRaw = fdrCol.getRawData() as Float32Array | Float64Array;
  const sigCol = df.columns.addNewBool('Significant');
  sigCol.init((i) => fdrRaw[i] !== DG.FLOAT_NULL && fdrRaw[i] <= pThreshold);

  // Optional directional label (R2 up/down split). Categorical, bulk-init
  // (memory feedback_dg_column_bulk_init — never per-row set).
  if (direction) {
    const directionCol = df.columns.addNewString('Direction');
    directionCol.init(() => direction);
  }

  // Color coding on FDR column
  fdrCol.setTag('color-coding-type', 'Linear');
  fdrCol.setTag('color-coding-linear', '{"0":"#2ecc71","0.05":"#f39c12","1":"#e74c3c"}');

  df.setTag('proteomics.enrichment', 'true');
  return df;
}

/**
 * Splits the detected genes into up/down significant query sets and the shared
 * all-detected background. Ported from CK-omics run_gprofiler_analysis: strict
 * `fc > 0` / `fc < 0`; background is every detected gene (used by BOTH
 * directional g:Profiler calls). Null fc/p rows count toward background only.
 */
export function splitGenesByDirection(
  geneForRow: Map<number, string>,
  fcRaw: Float32Array | Float64Array,
  pRaw: Float32Array | Float64Array,
  fcThreshold: number,
  pThreshold: number,
): {upGenes: string[]; downGenes: string[]; background: string[]} {
  const up: Set<string> = new Set();
  const down: Set<string> = new Set();
  const background: Set<string> = new Set();
  for (const [row, gene] of geneForRow) {
    background.add(gene);
    const fc = fcRaw[row];
    const adjP = pRaw[row];
    if (fc !== DG.FLOAT_NULL && adjP !== DG.FLOAT_NULL &&
        adjP <= pThreshold && Math.abs(fc) >= fcThreshold) {
      if (fc > 0) up.add(gene);
      else if (fc < 0) down.add(gene);
    }
  }
  return {upGenes: [...up], downGenes: [...down], background: [...background]};
}

// --- Significant Protein Counter ---

export function countSignificantProteins(
  df: DG.DataFrame,
  fcThreshold: number,
  pThreshold: number,
  log2fcCol: DG.Column,
  pValCol: DG.Column,
): {significant: number; total: number} {
  // Bulk-read via getRawData so we don't pay a 3× JS/native bridge hop per row.
  const n = df.rowCount;
  const fcRaw = log2fcCol.getRawData() as Float32Array | Float64Array;
  const pRaw = pValCol.getRawData() as Float32Array | Float64Array;
  let significant = 0;
  for (let i = 0; i < n; i++) {
    const fc = fcRaw[i];
    const p = pRaw[i];
    if (fc === DG.FLOAT_NULL || p === DG.FLOAT_NULL) continue;
    if (Math.abs(fc) >= fcThreshold && p <= pThreshold)
      significant++;
  }
  return {significant, total: n};
}

// --- Pipeline Orchestrator ---

export async function runEnrichmentPipeline(
  df: DG.DataFrame,
  fcThreshold: number,
  pThreshold: number,
  organismCode: string,
  sources: string[],
  smartFilterEnabled: boolean = true,
  maxTermsPerSource: number = 15,
): Promise<{enrichmentDf: DG.DataFrame; mapped: number; total: number; unmapped: number}> {
  const cols = findProteomicsColumns(df);
  if (!cols.proteinId)
    throw new Error('No protein ID column found');
  if (!cols.log2fc || !cols.pValue)
    throw new Error('No log2FC or adjusted p-value columns found. Run Differential Expression first.');

  // Step 1: Get gene symbols (existing or mapped)
  let geneForRow: Map<number, string>;
  let mapped = 0;
  let total = 0;
  let unmapped = 0;

  if (cols.geneName) {
    // Use existing gene name column. Truthiness check covers null/empty;
    // skipping the redundant isNone() halves the bridge hops.
    geneForRow = new Map();
    const geneNameCol = cols.geneName;
    for (let i = 0; i < df.rowCount; i++) {
      const val = geneNameCol.get(i) as string | null;
      if (val) geneForRow.set(i, val);
    }
    mapped = geneForRow.size;
    total = df.rowCount;
    unmapped = total - mapped;
  } else {
    // Extract clean accessions and map via g:Convert
    const accToRows = new Map<string, number[]>();
    for (let i = 0; i < df.rowCount; i++) {
      const raw = cols.proteinId!.get(i) as string;
      const acc = parseAccession(raw);
      if (acc) {
        if (!accToRows.has(acc))
          accToRows.set(acc, []);
        accToRows.get(acc)!.push(i);
      }
    }

    const uniqueAccessions = [...accToRows.keys()];
    const convertResults = await gConvert(uniqueAccessions, organismCode);

    // Build accession -> gene symbol map (first match per accession)
    const accToGene = new Map<string, string>();
    for (const r of convertResults) {
      if (r.incoming && r.name && r.name !== 'N/A' && !accToGene.has(r.incoming))
        accToGene.set(r.incoming, r.name);
    }

    // Map rows to gene symbols
    geneForRow = new Map();
    for (const [acc, rows] of accToRows) {
      const gene = accToGene.get(acc);
      if (gene) {
        for (const row of rows)
          geneForRow.set(row, gene);
      }
    }

    // Add gene symbol column to main table. Materialize the row->gene array
    // first, then bulk-init the column to skip per-row set() bridge hops.
    const geneArr: string[] = new Array(df.rowCount).fill('');
    for (const [row, gene] of geneForRow)
      geneArr[row] = gene;

    const geneCol = df.columns.addNewString('Gene Symbol (mapped)');
    geneCol.semType = SEMTYPE.GENE_SYMBOL;
    geneCol.init((i) => geneArr[i]);

    mapped = accToGene.size;
    total = uniqueAccessions.length;
    unmapped = total - mapped;
  }

  // Step 2: Split significant genes by fc sign; background = all detected,
  // shared by both directional queries (CK-omics run_gprofiler_analysis).
  // Bulk-read the numeric columns once via getRawData (no per-row bridge hops).
  const fcRaw = cols.log2fc.getRawData() as Float32Array | Float64Array;
  const pRaw = cols.pValue.getRawData() as Float32Array | Float64Array;
  const {upGenes, downGenes, background: bgArray} =
    splitGenesByDirection(geneForRow, fcRaw, pRaw, fcThreshold, pThreshold);

  if (upGenes.length === 0 && downGenes.length === 0)
    throw new Error('No significant proteins found with the given thresholds');

  // Step 3: One g:GOSt call per non-empty direction, identical background.
  // A direction with zero genes is skipped (no empty-query request).
  // Pitfall 8: when smart filter is enabled, run it BEFORE buildEnrichmentDf so
  // the per-row Intersection strings are built from the kept-set indices and
  // never get reindexed after the fact.
  const directionDfs: DG.DataFrame[] = [];
  let upStats: SmartFilterStats | null = null;
  let downStats: SmartFilterStats | null = null;
  if (upGenes.length > 0) {
    let upResults = await gGOSt(upGenes, bgArray, organismCode, sources, pThreshold);
    if (smartFilterEnabled) {
      const filtered = applySmartPathwayFilter(upResults, maxTermsPerSource);
      upResults = filtered.kept;
      upStats = filtered.stats;
    }
    directionDfs.push(buildEnrichmentDf(upResults, upGenes, pThreshold, 'Up'));
  }
  if (downGenes.length > 0) {
    let downResults = await gGOSt(downGenes, bgArray, organismCode, sources, pThreshold);
    if (smartFilterEnabled) {
      const filtered = applySmartPathwayFilter(downResults, maxTermsPerSource);
      downResults = filtered.kept;
      downStats = filtered.stats;
    }
    directionDfs.push(buildEnrichmentDf(downResults, downGenes, pThreshold, 'Down'));
  }

  // Step 4: Merge directional results into one DataFrame. Phase-9 column shape
  // is intact (Direction is just an extra categorical). append() returns a NEW
  // frame and may drop frame/column metadata, so re-apply the enrichment tag,
  // name, and FDR color-coding the viewers/cross-link key on.
  let enrichmentDf = directionDfs[0];
  for (let i = 1; i < directionDfs.length; i++)
    enrichmentDf = enrichmentDf.append(directionDfs[i]);
  // Name the result after its source table (e.g. "MyStudy — Enrichment") so it's
  // distinguishable when several analyses are open — mirrors the "<source> — QC"
  // convention. Falls back to the generic name if the source is unnamed.
  enrichmentDf.name = df.name ? `${df.name} — Enrichment` : 'Enrichment Results';
  enrichmentDf.setTag('proteomics.enrichment', 'true');
  const fdrCol = enrichmentDf.col('FDR');
  if (fdrCol) {
    fdrCol.setTag('color-coding-type', 'Linear');
    fdrCol.setTag('color-coding-linear', '{"0":"#2ecc71","0.05":"#f39c12","1":"#e74c3c"}');
  }

  // D-14 — set the smart-filter banner tags only when the filter actually
  // dropped rows (kept < total). Banner is suppressed when the filter was off
  // OR when kept == total. Tag names mirror UI-SPEC §"Enrichment table:
  // smart pathway filter banner".
  if (smartFilterEnabled) {
    const totalInput = (upStats?.total ?? 0) + (downStats?.total ?? 0);
    const totalKept = (upStats?.kept ?? 0) + (downStats?.kept ?? 0);
    const totalDropped = (upStats?.droppedParents ?? 0) + (downStats?.droppedParents ?? 0);
    if (totalKept < totalInput) {
      enrichmentDf.setTag('proteomics.enrichment_smart_filtered', 'true');
      enrichmentDf.setTag('proteomics.enrichment_smart_filtered_kept', String(totalKept));
      enrichmentDf.setTag('proteomics.enrichment_smart_filtered_total', String(totalInput));
      enrichmentDf.setTag('proteomics.enrichment_smart_filtered_dropped_parents', String(totalDropped));
      enrichmentDf.setTag('proteomics.enrichment_smart_filtered_cap', String(maxTermsPerSource));
    }
  }

  return {enrichmentDf, mapped, total, unmapped};
}

// --- Dialog ---

export function showEnrichmentDialog(df: DG.DataFrame): void {
  if (!requireDifferentialExpression(df, 'Run Differential Expression first')) return;

  const cols = findProteomicsColumns(df);
  if (!cols.log2fc || !cols.pValue) {
    grok.shell.warning('No log2FC or p-value columns found');
    return;
  }

  const fcInput = ui.input.float('|log2FC| threshold', {value: DEFAULT_FC_THRESHOLD});
  fcInput.setTooltip('Minimum absolute fold change for significance');

  const pInput = ui.input.float('Adj. p-value threshold', {value: DEFAULT_P_THRESHOLD});
  pInput.setTooltip('Maximum adjusted p-value for significance');

  // Seed from the chosen/persisted organism, else auto-detect from the data's
  // organism column (covers the Candidates path, which skips Annotate), else human.
  // Sticky + shared with the subcellular-location fetch via the proteomics.organism tag.
  const seededDisplay = organismDisplayForCode(getOrganism(df) ?? detectOrganismCode(df));
  const organismInput = ui.input.choice('Organism', {
    value: (seededDisplay ?? 'Homo sapiens (Human)') as string,
    items: ORGANISM_LIST.map((o) => o.display) as unknown as string[],
    nullable: false,
  });

  // Enrichment source checkboxes
  const goBpInput = ui.input.bool('GO: Biological Process', {value: true});
  const goMfInput = ui.input.bool('GO: Molecular Function', {value: true});
  const goCcInput = ui.input.bool('GO: Cellular Component', {value: true});
  const keggInput = ui.input.bool('KEGG Pathways', {value: true});
  const reactomeInput = ui.input.bool('Reactome Pathways', {value: true});
  const wpInput = ui.input.bool('WikiPathways', {value: true});

  // D-14 — smart pathway filter checkbox, default checked.
  const smartFilterInput = ui.input.bool('Apply smart pathway filter', {value: true});
  smartFilterInput.setTooltip(
    'Drops generic GO parent terms when more specific child terms are present; ' +
    'caps each source at top-N by FDR. Recommended for cleaner results; ' +
    'uncheck for raw g:Profiler output.');

  // M3 — the per-source top-N cap the smart filter applies (was hardcoded 15).
  // Only meaningful when the smart filter is on; greyed out otherwise.
  const maxTermsInput = ui.input.int('Max terms per source', {value: 15, min: 1});
  maxTermsInput.setTooltip(
    'Upper bound on terms kept per enrichment source (top-N by FDR). ' +
    'Applies only when the smart pathway filter is on.');
  const syncMaxTermsEnabled = () => maxTermsInput.enabled = smartFilterInput.value ?? true;
  smartFilterInput.onChanged.subscribe(syncMaxTermsEnabled);
  syncMaxTermsEnabled();

  // Live count of significant proteins
  const countDiv = ui.divText('');
  const updateCount = () => {
    const fc = fcInput.value ?? DEFAULT_FC_THRESHOLD;
    const p = pInput.value ?? DEFAULT_P_THRESHOLD;
    const result = countSignificantProteins(df, fc, p, cols.log2fc!, cols.pValue!);
    const totalFmt = result.total.toLocaleString();
    countDiv.textContent =
      `${result.significant} of ${totalFmt} proteins are significant with these thresholds`;
  };
  updateCount();
  fcInput.onChanged.subscribe(updateCount);
  pInput.onChanged.subscribe(updateCount);

  ui.dialog('Enrichment Analysis')
    .add(countDiv)
    .add(fcInput)
    .add(pInput)
    .add(organismInput)
    .add(goBpInput)
    .add(goMfInput)
    .add(goCcInput)
    .add(keggInput)
    .add(reactomeInput)
    .add(wpInput)
    .add(smartFilterInput)
    .add(maxTermsInput)
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('Running enrichment analysis...');
      try {
        const fc = fcInput.value ?? DEFAULT_FC_THRESHOLD;
        const p = pInput.value ?? DEFAULT_P_THRESHOLD;

        // Resolve organism code from display name
        const selectedDisplay = organismInput.value;
        const organism = ORGANISM_LIST.find((o) => o.display === selectedDisplay);
        const organismCode = organism?.code ?? 'hsapiens';
        // Persist so the subcellular-location fetch (and next enrichment run)
        // narrow to the same species.
        setOrganism(df, organismCode);

        // Build sources array from checkbox states
        const selectedSources: string[] = [];
        if (goBpInput.value) selectedSources.push('GO:BP');
        if (goMfInput.value) selectedSources.push('GO:MF');
        if (goCcInput.value) selectedSources.push('GO:CC');
        if (keggInput.value) selectedSources.push('KEGG');
        if (reactomeInput.value) selectedSources.push('REAC');
        if (wpInput.value) selectedSources.push('WP');

        if (selectedSources.length === 0) {
          grok.shell.warning('Please select at least one enrichment source');
          return;
        }

        const result = await runEnrichmentPipeline(
          df, fc, p, organismCode, selectedSources, smartFilterInput.value ?? true,
          maxTermsInput.value ?? 15);

        // Show mapping stats
        const pct = result.total > 0 ? (result.mapped / result.total * 100).toFixed(1) : '0.0';
        grok.shell.info(
          `${result.mapped}/${result.total} proteins mapped (${pct}%). ${result.unmapped} unmapped.`,
        );

        // openEnrichmentVisualization creates (or reuses) the enrichment results
        // view itself — adding one here too spawned a duplicate "Enrichment
        // Results (2)" tab every run.
        openEnrichmentVisualization(result.enrichmentDf, df);
      } catch (e: any) {
        grok.shell.error(`Enrichment analysis failed: ${e.message}`);
      } finally {
        pi.close();
      }
    })
    .show();
}
