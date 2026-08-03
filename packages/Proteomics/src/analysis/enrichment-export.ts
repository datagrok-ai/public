import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {findProteomicsColumns} from '../utils/column-detection';
import {DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD} from '../utils/proteomics-types';
import {splitGenesByDirection, countSignificantProteins} from './enrichment';

/**
 * Exports the exact inputs an enrichment run would use — the up-regulated gene
 * list, the down-regulated gene list, and the background (every detected gene) —
 * so the analyst can reproduce or hand off the g:Profiler query outside Datagrok
 * (CK-omics CKomics_tool2.py:4531-4640 wrote the same three lists).
 *
 * Deliberately network-free: it reads the resolved gene-name column (falling
 * back to protein IDs) rather than re-running g:Convert, mirroring the pipeline's
 * non-mapping path via the SAME `splitGenesByDirection` so the exported lists
 * match what Enrichment Analysis would query at identical thresholds.
 */

/** Tidy long-format CSV: `gene,list` with list ∈ {up, down, background}. Each
 * list is filter-and-copy friendly for g:Profiler; a gene appears under
 * `background` and (if significant) also under `up`/`down` — faithful to
 * g:Profiler semantics where the background is the full universe. */
export function buildEnrichmentInputsCsv(
  df: DG.DataFrame, fcThreshold: number, pThreshold: number,
): {csv: string; counts: {up: number; down: number; background: number}} {
  const cols = findProteomicsColumns(df);
  if (!cols.log2fc || !cols.pValue)
    throw new Error('No log2FC or adjusted p-value columns found. Run Differential Expression first.');

  // Prefer resolved gene symbols; fall back to protein IDs so the export still
  // works before gene mapping (honest — the header says what was used).
  const idCol = cols.geneName ?? cols.proteinId;
  if (!idCol)
    throw new Error('No gene-name or protein-ID column found to export.');

  const geneForRow = new Map<number, string>();
  for (let i = 0; i < df.rowCount; i++) {
    const v = idCol.get(i) as string | null;
    if (v) geneForRow.set(i, String(v));
  }

  const fcRaw = cols.log2fc.getRawData() as Float32Array | Float64Array;
  const pRaw = cols.pValue.getRawData() as Float32Array | Float64Array;
  const {upGenes, downGenes, background} =
    splitGenesByDirection(geneForRow, fcRaw, pRaw, fcThreshold, pThreshold);

  const esc = (s: string): string =>
    /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  const rows: string[] = ['gene,list'];
  for (const g of upGenes) rows.push(`${esc(g)},up`);
  for (const g of downGenes) rows.push(`${esc(g)},down`);
  for (const g of background) rows.push(`${esc(g)},background`);

  return {
    csv: rows.join('\n'),
    counts: {up: upGenes.length, down: downGenes.length, background: background.length},
  };
}

/** Small dialog to pick significance thresholds (seeded from the package
 * defaults, matching the Enrichment dialog) and download the input lists as one
 * CSV. Immediate download on OK via `DG.Utils.download`. */
export function showEnrichmentInputExportDialog(df: DG.DataFrame): void {
  const cols = findProteomicsColumns(df);
  if (!cols.log2fc || !cols.pValue) {
    grok.shell.warning('No log2FC or p-value columns found. Run Differential Expression first.');
    return;
  }

  const fcInput = ui.input.float('|log2FC| threshold', {value: DEFAULT_FC_THRESHOLD});
  fcInput.setTooltip('Minimum absolute fold change for a gene to count as up/down.');
  const pInput = ui.input.float('Adj. p-value threshold', {value: DEFAULT_P_THRESHOLD});
  pInput.setTooltip('Maximum adjusted p-value for a gene to count as up/down.');

  const countDiv = ui.divText('');
  const updateCount = () => {
    const r = countSignificantProteins(df, fcInput.value ?? DEFAULT_FC_THRESHOLD,
      pInput.value ?? DEFAULT_P_THRESHOLD, cols.log2fc!, cols.pValue!);
    countDiv.textContent = `${r.significant} of ${r.total.toLocaleString()} proteins significant`;
  };
  updateCount();
  fcInput.onChanged.subscribe(updateCount);
  pInput.onChanged.subscribe(updateCount);

  ui.dialog('Export Enrichment Inputs')
    .add(countDiv)
    .add(fcInput)
    .add(pInput)
    .onOK(() => {
      const fc = fcInput.value ?? DEFAULT_FC_THRESHOLD;
      const p = pInput.value ?? DEFAULT_P_THRESHOLD;
      const {csv, counts} = buildEnrichmentInputsCsv(df, fc, p);
      const base = df.name ? `${df.name} — enrichment inputs` : 'enrichment inputs';
      DG.Utils.download(`${base}.csv`, csv, 'text/csv');
      grok.shell.info(
        `Exported ${counts.up} up, ${counts.down} down, ${counts.background} background genes.`);
    })
    .show();
}
