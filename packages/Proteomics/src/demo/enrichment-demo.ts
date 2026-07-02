import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {setGroups} from '../analysis/experiment-setup';
import {runEnrichmentPipeline} from '../analysis/enrichment';
import {openEnrichmentVisualization} from '../viewers/enrichment-viewers';
import {createVolcanoPlot} from '../viewers/volcano';
import {focusProtein} from '../panels/protein-focus';
import {SEMTYPE} from '../utils/proteomics-types';

/**
 * Single-organism pathway-enrichment demo. Loads a small engineered HUMAN
 * differential-expression result (Treatment vs Control) whose regulated proteins
 * form two coherent, textbook-enriched sets — cell-cycle / mitosis UP and
 * oxidative phosphorylation DOWN — over a diverse non-significant background, so
 * g:Profiler returns crisp, recognizable terms. It opens the volcano, runs
 * enrichment (GO / KEGG / Reactome / WikiPathways), and docks the enrichment
 * charts cross-linked to the volcano: selecting a term highlights its proteins.
 *
 * This is the counterpart to the main `Proteomics Demo`, which uses the
 * three-species HYE benchmark — a mix that is deliberately unsuitable for
 * enrichment (g:Profiler assumes a single organism). The dataset ships at
 * `files/demo/enrichment-demo.csv` (see `files/demo/README.md`).
 *
 * Enrichment calls g:Profiler over the network; if the service is unreachable
 * the demo still opens the DE table + volcano and explains how to retry.
 */
export async function runEnrichmentDemo(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Loading enrichment demo…');
  try {
    pi.update(10, 'Loading human DE result…');
    const text = await _package.files.readAsText('demo/enrichment-demo.csv');
    const df = DG.DataFrame.fromCsv(text);
    df.name = 'Enrichment demo (Treatment vs Control)';

    // Semantic-type the columns so the pipeline and panels find them by type,
    // not by name — mirrors what the parsers do on a real import.
    const setSem = (name: string, sem: string): void => {
      const c = df.col(name);
      if (c) c.semType = sem;
    };
    setSem('Protein ID', SEMTYPE.PROTEIN_ID);
    setSem('Gene', SEMTYPE.GENE_SYMBOL);
    setSem('log2FC', SEMTYPE.LOG2FC);
    setSem('p-value', SEMTYPE.P_VALUE);
    setSem('adj.p-value', SEMTYPE.P_VALUE);

    // A pre-computed DE result (like a Spectronaut Candidates import): mark DE
    // complete and record group NAMES so the volcano legend reads naturally.
    // There are no per-sample columns, so the group column lists stay empty.
    df.setTag('proteomics.source', 'generic');
    df.setTag('proteomics.de_method', 'precomputed');
    df.setTag('proteomics.de_complete', 'true');
    setGroups(df, {group1: {name: 'Treatment', columns: []}, group2: {name: 'Control', columns: []}});

    const tv = grok.shell.addTableView(df);
    const volcano = createVolcanoPlot(df);
    tv.dockManager.dock(volcano, DG.DOCK_TYPE.RIGHT, null, 'Volcano', 0.5);

    // Headline protein (row 0 = a cell-cycle hit) selected for its UniProt panel.
    focusProtein(df, 0);

    pi.update(45, 'Running g:Profiler enrichment (GO / KEGG / Reactome)…');
    try {
      const {enrichmentDf, mapped, total} = await runEnrichmentPipeline(
        df, 1.0, 0.05, 'hsapiens',
        ['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'WP'], true);
      pi.update(90, 'Opening enrichment charts…');
      openEnrichmentVisualization(enrichmentDf, df);
      grok.shell.info(
        `Enrichment demo: ${enrichmentDf.rowCount} over-represented terms across ${mapped}/${total} genes. ` +
        'Up-regulated proteins are enriched for cell-cycle / mitosis; down-regulated for oxidative ' +
        'phosphorylation. Select a term row to highlight its proteins on the volcano.');
    } catch (e: any) {
      grok.shell.warning(
        `Enrichment demo: g:Profiler enrichment could not run (${e?.message ?? e}). ` +
        'The DE table and volcano are open — retry via Proteomics | Enrichment Analysis ' +
        'once the service is reachable.');
    }

    pi.update(100, 'Done');
  } finally {
    pi.close();
  }
}
