import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {parseSpectronautText} from '../parsers/spectronaut-parser';
import {getGroups} from '../analysis/experiment-setup';
import {imputeKnn} from '../analysis/imputation';
import {runDifferentialExpression} from '../analysis/differential-expression';
import {createVolcanoPlot} from '../viewers/volcano';
import {createExpressionHeatmap} from '../viewers/heatmap';
import {openQcDashboard} from '../viewers/qc-dashboard';
import {focusProtein} from '../panels/protein-focus';

/**
 * The package's single "richest endpoint" demo. It takes the HYE benchmark (a
 * Spectronaut DIA report where Human, Yeast, and E. coli proteomes are mixed at
 * known, different ratios between two conditions — HYE mix A vs B, four
 * replicates each) all the way through the analysis pipeline and lands on the
 * fully-analyzed protein table with the Volcano and Expression Heatmap docked,
 * a QC dashboard on a second tab, and the most significant protein pre-selected
 * so its UniProt entry is already showing in the context panel.
 *
 * HYE is chosen deliberately: two whole species shift abundance between the
 * conditions, so hundreds of proteins are genuinely differential and the
 * volcano shows a strong, self-evident signal — the point of a "richest
 * endpoint" demo. (A spike-in benchmark like CPTAC is real but deliberately
 * subtle: at n=3 a plain t-test can't clear FDR, leaving an all-grey volcano.)
 *
 * Enrichment is intentionally NOT part of this demo: g:Profiler enrichment
 * assumes a single organism, and HYE is a three-species mix, so GO/pathway
 * enrichment would only map the human subset and mislead. Enrichment belongs in
 * a single-organism demo.
 *
 * The Spectronaut parser auto-annotates the two conditions from R.Condition and
 * tags the data preNormalized, so the demo skips straight to imputation + DE and
 * does NOT re-normalize — re-normalizing an already-normalized export would
 * distort the known species-ratio signal. Differential expression runs the
 * client-side Welch t-test, so the demo needs no R environment. The only network
 * touches are best-effort and cached: the parser's Ensembl gene-label lookup and
 * the UniProt context panel, both of which degrade gracefully offline. The
 * dataset ships at `files/demo/spectronaut-hye-demo.tsv` (see
 * `files/demo/README.md`).
 */
export async function runProteomicsDemo(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('Loading proteomics demo…');
  try {
    pi.update(5, 'Importing HYE benchmark (Spectronaut DIA)…');
    const text = await _package.files.readAsText('demo/spectronaut-hye-demo.tsv');
    const df = await parseSpectronautText(text);
    df.name = 'HYE benchmark';
    const analysisTv = grok.shell.addTableView(df);

    // Groups (HYE mix A / B) are auto-populated by the parser from R.Condition.
    const groups = getGroups(df);
    if (!groups) {
      grok.shell.warning('Demo: expected two auto-annotated conditions — ' +
        'opened the table without running differential expression.');
      return;
    }
    const allCols = [...groups.group1.columns, ...groups.group2.columns];

    pi.update(40, 'Imputing missing values (KNN)…');
    imputeKnn(df, allCols);
    df.setTag('proteomics.imputed', 'true');

    pi.update(60, 'Differential expression (t-test)…');
    runDifferentialExpression(df, groups.group1.columns, groups.group2.columns,
      groups.group1.name, groups.group2.name);
    df.setTag('proteomics.de_method', 't-test');

    pi.update(78, 'Building volcano and heatmap…');
    const volcano = createVolcanoPlot(df);
    analysisTv.dockManager.dock(volcano, DG.DOCK_TYPE.RIGHT, null, 'Volcano', 0.5);

    const heatmap = await createExpressionHeatmap(df, {title: 'Heatmap: Top 50 DE Proteins'});
    analysisTv.dockManager.dock(heatmap, DG.DOCK_TYPE.DOWN, null, 'Heatmap', 0.45);

    // QC dashboard on a SECOND tab. openQcDashboard docks into grok.shell.tv, so
    // open a second table view over the same DataFrame (it becomes current), build
    // QC there, then land back on the analysis endpoint. Sharing one DataFrame
    // across two views is fine — QC adds its own columns; the volcano keeps its own.
    pi.update(90, 'Building QC dashboard tab…');
    const qcTv = grok.shell.addTableView(df);
    qcTv.name = 'HYE benchmark — QC';
    openQcDashboard(df);
    grok.shell.v = analysisTv;

    // Pre-select the most significant protein so the UniProt context panel is
    // already showing — a real accession from a spiked species, not the flat
    // human background. Shared with the import handlers (which focus row 0).
    focusProtein(df, topHitRow(df));

    pi.update(100, 'Done');
    grok.shell.info('Proteomics demo: HYE benchmark analyzed end to end ' +
      '(import → impute → differential expression). Yeast and E. coli proteins ' +
      'shift abundance between the two conditions, so the volcano separates the ' +
      'genuinely differential proteins from the constant human background. ' +
      'QC charts are on the second tab; the top hit is selected for its UniProt panel.');
  } finally {
    pi.close();
  }
}

/**
 * Row index of the most significant protein (smallest adj.p among the
 * significant rows) to pre-select for the UniProt panel. Falls back to row 0
 * when there is no significant protein or no adj.p column.
 */
function topHitRow(df: DG.DataFrame): number {
  const adj = df.col('adj.p-value');
  const sig = df.col('significant');
  if (!adj) return 0;

  let best = -1;
  let bestP = Infinity;
  for (let i = 0; i < df.rowCount; i++) {
    if (sig && !sig.get(i)) continue;
    if (adj.isNone(i)) continue;
    const p = adj.get(i) as number;
    if (p < bestP) { bestP = p; best = i; }
  }
  return best < 0 ? 0 : best;
}
