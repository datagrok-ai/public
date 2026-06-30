/**
 * Dynamic half of the **Proteomics** menu — the analysis groups that need a per-table
 * grey-out, built onto a table view's ribbon.
 *
 * The menu is split across two mechanisms on purpose:
 *
 *   - In the **main Datagrok menu (left bar)**: only `Proteomics | Import`, via decorator
 *     `top-menu` strings in `package.ts`. This is the always-available entry point,
 *     reachable with no table open. The platform owns it — robust, never duplicates.
 *   - In a **table view's ribbon**: the full menu (Import + Annotate / Analyze / Visualize
 *     / Share), built here onto `view.ribbonMenu`, and ONLY when a Proteomics analysis
 *     table is active. The analysis items need `isEnabled` grey-out (a Spectronaut
 *     Candidates file has no per-sample intensities, so the sample-level steps don't apply),
 *     which only works on a programmatically built menu.
 *
 * The two hosts are independent: a decorator `top-menu` populates `grok.shell.topMenu` (the
 * left-bar menu), NOT `view.ribbonMenu`. When a view's ribbon has its own 'Proteomics'
 * group, the ribbon shows THAT one — so Import has to be (re)added here for the ribbon copy.
 * They don't collide because they live in different menus.
 */
import * as DG from 'datagrok-api/dg';

/** Handlers the ribbon menu invokes — including Import, which the ribbon builds itself
 * (the decorator copies only reach the left-bar main menu). Passed in from `package.ts`
 * so this module doesn't import `PackageFunctions` (circular). */
export interface ProteomicsMenuHandlers {
  importSpectronautCandidates: () => void;
  importSpectronaut: () => void;
  importMaxQuant: () => void;
  importFragPipe: () => void;
  importGenericMatrix: () => void;
  annotateExperiment: () => void;
  normalize: () => void;
  impute: () => void;
  differentialExpression: () => void;
  enrichmentAnalysis: () => void;
  computeSpcStatus: () => void;
  showVolcanoPlot: () => void;
  showHeatmap: () => void;
  showPcaPlot: () => void;
  showGroupMeanCorrelation: () => void;
  showQcDashboard: () => void;
  showSpcDashboard: () => void;
  enrichmentCharts: () => void;
  showAllVisualizations: () => void;
  shareAnalysisForReview: () => void;
}

/** Our subgroups/items under the ribbon 'Proteomics' group. Cleared before a (re)build so
 * repeated calls replace rather than duplicate. 'Import' is ours here too — the decorator
 * Import lives only in the left-bar main menu, not the ribbon. */
const DYNAMIC_SECTIONS = ['Import', 'Annotate Experiment...', 'Analyze', 'Visualize', 'Share'];

/**
 * Is this an analysis table the full Proteomics menu should attach to? Every parser stamps
 * `proteomics.source`, so its presence is the reliable signal. Plain CSVs and the
 * sample-level PCA table (which has no source tag) get only the decorator Import.
 */
export function isProteomicsTable(df: DG.DataFrame | null | undefined): boolean {
  return !!df && !!df.getTag('proteomics.source');
}

/**
 * Disable reason for sample-level menu items on `df`, or `null` when they're fine.
 * Re-evaluated by the platform every time the menu opens, so it tracks the live table.
 * Only a Spectronaut Candidates table is blocked; the per-handler runtime guards stay in
 * place as the backstop for every other not-ready case.
 */
export function sampleLevelDisabledReason(df: DG.DataFrame | null | undefined): string | null {
  if (df && df.getTag('proteomics.source') === 'spectronaut-candidates') {
    return 'Not applicable to a Spectronaut Candidates file — it already carries computed ' +
      'differential expression and has no per-sample intensities. Use Volcano Plot, ' +
      'Enrichment Analysis or the UniProt panel, or import the matching Spectronaut Report ' +
      'for sample-level analyses.';
  }
  return null;
}

/**
 * Build the dynamic analysis groups onto a table view's ribbon, attaching to the Proteomics
 * group the decorator Import created. Returns `true` if it built (the view is a Proteomics
 * table), `false` otherwise — the caller uses that to decide whether the view is "done".
 *
 * Idempotent: clears its own sections first (never the decorator 'Import') so a rebuild
 * replaces rather than duplicates.
 */
export function buildProteomicsRibbonMenu(view: DG.TableView, h: ProteomicsMenuHandlers): boolean {
  const df = view.dataFrame;
  if (!isProteomicsTable(df))
    return false;

  const root = view.ribbonMenu.group('Proteomics');
  for (const s of DYNAMIC_SECTIONS) {
    try { root.remove(s); } catch { /* not present — fine */ }
  }

  const sampleOnly = {isEnabled: () => sampleLevelDisabledReason(view.dataFrame)};

  // Import — the ribbon's own copy (the decorator Import only reaches the left-bar menu).
  const imp = root.group('Import');
  imp.item('Spectronaut Candidates...', h.importSpectronautCandidates);
  imp.item('Spectronaut Report...', h.importSpectronaut);
  imp.item('MaxQuant...', h.importMaxQuant);
  imp.item('FragPipe...', h.importFragPipe);
  imp.item('Generic Matrix...', h.importGenericMatrix);

  root.item('Annotate Experiment...', h.annotateExperiment, null, sampleOnly);

  const analyze = root.group('Analyze');
  analyze.item('Normalize...', h.normalize, null, sampleOnly);
  analyze.item('Impute Missing Values...', h.impute, null, sampleOnly);
  analyze.item('Differential Expression...', h.differentialExpression, null, sampleOnly);
  analyze.item('Enrichment Analysis...', h.enrichmentAnalysis);
  analyze.item('Compute SPC Status', h.computeSpcStatus, null, sampleOnly);

  const viz = root.group('Visualize');
  viz.item('Volcano Plot...', h.showVolcanoPlot);
  viz.item('Heatmap...', h.showHeatmap, null, sampleOnly);
  viz.item('PCA...', h.showPcaPlot, null, sampleOnly);
  viz.item('Group-Mean Correlation...', h.showGroupMeanCorrelation, null, sampleOnly);
  viz.item('QC Dashboard...', h.showQcDashboard, null, sampleOnly);
  viz.item('SPC Dashboard...', h.showSpcDashboard);
  viz.item('Enrichment Charts...', h.enrichmentCharts);
  viz.item('Show All Visualizations...', h.showAllVisualizations, null, sampleOnly);

  const share = root.group('Share');
  share.item('Share Analysis for Review...', h.shareAnalysisForReview);
  return true;
}
