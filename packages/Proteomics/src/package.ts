/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {debounceTime} from 'rxjs/operators';

import {parseMaxQuantText} from './parsers/maxquant-parser';
import {parseSpectronautText, parseSpectronautStream} from './parsers/spectronaut-parser';
import {parseSpectronautCandidatesText, COMPARISON_COLUMNS}
  from './parsers/spectronaut-candidates-parser';
import {parseFragPipeText} from './parsers/fragpipe-parser';
import {showGenericImportDialog} from './parsers/generic-parser';
import {showAnnotationDialog, getGroups, getOrganism, setOrganism} from './analysis/experiment-setup';
import {detectOrganismCode} from './utils/organisms';
import {computeSpcMetrics, evaluateNelsonRulesAllMetrics, setSpcStatus, getRunMeta,
  defaultRulesEnabledAllMetrics, BaselineSnapshot} from './analysis/spc';
import {upsertRun, loadRuns, loadBaseline, SpcRunRow} from './analysis/spc-storage';
import {showNormalizationDialog, quantileNormalize, vsnNormalize} from './analysis/normalization';
import {showLog2ScaleDialog} from './analysis/log2-scale';
import {showEnrichmentInputExportDialog} from './analysis/enrichment-export';
import {showImputationDialog, imputeKnn, imputeZero, imputeMean, imputeMedian} from './analysis/imputation';
import {showDEDialog, requireDifferentialExpression} from './analysis/differential-expression';
import {createVolcanoPlot, recomputeVolcano, readVolcanoState, MetricKind, applyTopNLabels,
  getVolcanoTopN, showVolcanoBusy, updateVolcanoBusy, hideVolcanoBusy,
  getVolcanoAxisMax, setVolcanoAxisMax,
  LOCATION_COL} from './viewers/volcano';
import {STORE as SUBCELL_STORE} from './analysis/subcellular-location';
import {createExpressionHeatmap} from './viewers/heatmap';
import {createPcaPlot} from './viewers/pca-plot';
import {openQcDashboard} from './viewers/qc-dashboard';
import {openSpcDashboard} from './viewers/spc-dashboard';
import {createGroupMeanCorrelation} from './viewers/group-mean-correlation';
import {uniprotPanel} from './panels/uniprot-panel';
import {focusProtein} from './panels/protein-focus';
import {publishedAnalysisPanel} from './panels/published-analysis-panel';
import {showShareForReviewDialog} from './publishing/share-dialog';
import {recoverPublishedProject} from './publishing/post-open-recovery';
import {isPublished} from './publishing/publish-state';
import {showEnrichmentDialog} from './analysis/enrichment';
import {openEnrichmentVisualization} from './viewers/enrichment-viewers';
import {findColumn} from './utils/column-detection';
import {SEMTYPE, DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD} from './utils/proteomics-types';
import {buildProteomicsRibbonMenu} from './menu';
import {runProteomicsDemo} from './demo/proteomics-demo';
import {runEnrichmentDemo} from './demo/enrichment-demo';

export const _package = new DG.Package();
export * from './package.g';

/** Temporary polyfill */

function getDecoratorFunc() {
  return function(args: any) {
    return function(
      target: any,
      propertyKey: string,
      descriptor: PropertyDescriptor,
    ) { };
  };
}

if (!grok.decorators)
  (grok as any).decorators = {};

const decorators = [
  'func', 'init', 'param', 'panel', 'editor', 'demo', 'app',
  'appTreeBrowser', 'fileHandler', 'fileExporter', 'model', 'viewer', 'filter', 'cellRenderer', 'autostart',
  'dashboard', 'folderViewer', 'semTypeDetector', 'packageSettingsEditor', 'functionAnalysis', 'converter',
  'fileViewer', 'treeBrowser', 'polyfill',
];

decorators.forEach((decorator) => {
  if (!(grok.decorators as any)[decorator])
    (grok.decorators as any)[decorator] = getDecoratorFunc();
});

/** End temporary polyfill */


/** Guard for menu handlers that operate on per-sample intensity data. Returns
 * false (with a clear warning) when the open DataFrame is a pre-computed DE
 * output like a Spectronaut Candidates report — those have no per-sample
 * intensity columns, no group annotations, and no raw quantification to feed
 * into Annotate / Normalize / Impute / DE / Heatmap / PCA / QC. */
function requireSampleLevelData(df: DG.DataFrame, action: string): boolean {
  const source = df.getTag('proteomics.source');
  if (source === 'spectronaut-candidates') {
    grok.shell.warning(
      `${action} needs per-sample intensities, but this table is a Spectronaut Candidates ` +
      `report (one row per protein per comparison — pre-computed DE only). ` +
      `Import the matching Spectronaut Report (Proteomics | Import | Spectronaut Report) ` +
      `to enable per-sample analyses, or stick to Volcano / Enrichment / UniProt panel on this one.`,
    );
    return false;
  }
  return true;
}


/** Precursor/fragment-level signature columns. A Spectronaut report carrying
 * any of these is a long-format precursor export that must go through the
 * streaming aggregator (the V8 string ceiling makes file.text() OOM on the
 * 2.6 GB real report). A PG-level long report has none of these and keeps the
 * proven file.text() path. */
const PRECURSOR_SIGNATURE_COLUMNS = ['EG.ModifiedPeptide', 'FG.Charge', 'PEP.StrippedSequence'];

/** D-01 header sniff: streams only until the first newline (with a ~1 MB sanity
 * guard), then cancels the reader to release the stream — the rest of the file
 * is never read here. Returns true iff the header carries a precursor signature
 * column, in which case `importSpectronaut` routes to `parseSpectronautStream`. */
export async function sniffIsPrecursor(file: File): Promise<boolean> {
  const reader = file.stream()
    .pipeThrough(new TextDecoderStream('utf-8'))
    .getReader();
  try {
    let buffer = '';
    const maxBytes = 1024 * 1024; // 1 MB pre-newline sanity guard
    while (buffer.indexOf('\n') < 0 && buffer.length < maxBytes) {
      const {value, done} = await reader.read();
      if (done) break;
      buffer += value;
    }
    const nl = buffer.indexOf('\n');
    let header = nl >= 0 ? buffer.slice(0, nl) : buffer;
    if (header.endsWith('\r')) header = header.slice(0, -1);
    return PRECURSOR_SIGNATURE_COLUMNS.some((c) => header.includes(c));
  } finally {
    await reader.cancel();
  }
}


/** Module-level subscription store for the unified Filters viewer's
 * search-match → df.selection wiring. Disposed on re-entry so a second
 * multi-contrast import does not leak the prior handler. */
let activeFilterSubscriptions: rxjs.Subscription[] = [];

/** Sets `proteomics.organism` from the data's organism column at import when it
 * resolves to a single supported species, so enrichment and the subcellular-
 * location fetch narrow to the right organism without waiting for a dialog.
 * No-op when already set, absent, or ambiguous (multi-species) — the user then
 * picks in the Annotate/Enrichment dialog. */
function autoDetectOrganism(df: DG.DataFrame): void {
  if (getOrganism(df)) return;
  const code = detectOrganismCode(df);
  if (code) setOrganism(df, code);
}

/**
 * R4/D-07 + G4 + D-05: when a Spectronaut Candidates file carries more than one
 * distinct Comparison, dock a native Datagrok Filters viewer carrying typed
 * per-column filters: a categorical filter on Comparison plus free-text search
 * boxes on Display Name and Source ID (D-05 unified protein search). A single
 * distinct comparison docks nothing. Shell-only orchestration — the parser
 * stays pure. Returns whether a filter was docked.
 *
 * G4 root-cause fix (14-RESEARCH.md §"Pattern 1" / §"Pitfall 1"):
 * Phase 13 observed that the `columnNames` allowlist silently extends itself
 * with the boolean `Flags` column via the combined-boolean filter. We migrate
 * to the typed `filters` array and set `showBoolCombinedFilter: false` to
 * exclude Flags explicitly. A runtime verification gate falls back to an
 * explicit setOptions() override if the platform still leaks Flags into the
 * look config.
 */
export function dockComparisonFilterIfMultiContrast(
  tv: DG.TableView, df: DG.DataFrame,
): boolean {
  const cmpCol = COMPARISON_COLUMNS.reduce<DG.Column | null>(
    (found, name) => found ?? df.col(name), null);
  if (!cmpCol) return false;
  const distinct = new Set<string>();
  for (let i = 0; i < df.rowCount; i++) {
    const v = cmpCol.get(i);
    if (typeof v === 'string' && v.length > 0) distinct.add(v);
  }
  if (distinct.size <= 1) return false;

  // Display Name + Source ID come from Plan 14-01's gene-label resolver. On
  // DataFrames predating Plan 01 the columns may be absent — fall back to the
  // Protein ID column so the search-by-gene affordance still works.
  const displayNameCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name']);
  const sourceIdCol = findColumn(df, SEMTYPE.SOURCE_ID, ['source id']);
  const proteinIdCol = findColumn(df, SEMTYPE.PROTEIN_ID,
    ['primary protein id', 'protein id', 'uniprot', 'accession']);

  // Live UAT (2026-06-04) found the typed `filters: filterSpecs` shape doesn't
  // round-trip through the platform serializer — getOptions().look ends up
  // with ONLY `showBoolCombinedFilter`, no `filters` or `columnNames`, so the
  // viewer docks empty. The legacy `columnNames` array shape is the stable
  // contract: it materializes the requested columns as categorical filters
  // (free-text/categorical type is inferred per column). Display Name +
  // Source ID become string-categorical filters that include a search box
  // each — they replace the typed free-text filters with no UX regression.
  const columnNames: string[] = [cmpCol.name];
  if (displayNameCol) columnNames.push(displayNameCol.name);
  if (sourceIdCol) columnNames.push(sourceIdCol.name);
  if (!displayNameCol && proteinIdCol) columnNames.push(proteinIdCol.name);

  const filters = DG.Viewer.filters(df, {
    columnNames,
    showBoolCombinedFilter: false, // G4 Flags exclusion (14-RESEARCH §"Pitfall 1")
    showHeader: true,
    showSearchBox: true,
  } as any);

  tv.dockManager.dock(filters, DG.DOCK_TYPE.RIGHT, null, 'Filters', 0.3);

  // D-05 search-match wiring: free-text search boxes on Display Name / Source ID
  // mutate df.filter (hide non-matching rows). UI-SPEC requires highlight-not-
  // hide so the NS cloud stays visible. Defensive capture-restore path per
  // 14-RESEARCH §"Open Q 1 (RESOLVED)":
  //   1. Save df.filter into savedBitSet before subscribing.
  //   2. On df.onFilterChanged (debounced): capture matched indices from the
  //      current (search-mutated) df.filter.
  //   3. Restore df.filter via copyFrom(savedBitSet) so non-matched rows
  //      stay visible.
  //   4. Write df.selection.set(idx) on the matched indices.
  //   5. applyTopNLabels(df, sp, N, matched) so the top-N labels coexist with
  //      the search-matched rows (Pitfall 7). Labels are decoupled from
  //      selection — step 4's selection is the search highlight only.
  // The `restoring` flag guards re-entry from the copyFrom step.
  // Assumption A5: df.selection.set does NOT re-trigger df.filter.onChanged.
  for (const s of activeFilterSubscriptions) s.unsubscribe();
  activeFilterSubscriptions = [];

  const findVolcano = (): DG.ScatterPlotViewer | null => {
    for (const v of tv.viewers) {
      if (v.type === DG.VIEWER.SCATTER_PLOT &&
        (v as DG.ScatterPlotViewer).props.yColumnName === 'negLog10P')
        return v as DG.ScatterPlotViewer;
    }
    return null;
  };

  const savedBitSet = DG.BitSet.create(df.rowCount);
  savedBitSet.copyFrom(df.filter);
  let restoring = false;

  activeFilterSubscriptions.push(df.onFilterChanged.pipe(debounceTime(100)).subscribe(() => {
    if (restoring) return; // guard against the restore call re-triggering this handler

    // Capture matched indices from the current (search-mutated) df.filter.
    const matched: number[] = [];
    for (let i = 0; i < df.rowCount; i++) {
      if (df.filter.get(i)) matched.push(i);
    }

    // If nothing was filtered out (filter equals saved state), the change came
    // from elsewhere — nothing to highlight, just refresh the saved state.
    if (matched.length === df.rowCount) {
      savedBitSet.copyFrom(df.filter);
      return;
    }

    // Restore df.filter to its pre-search state.
    restoring = true;
    try {
      df.filter.copyFrom(savedBitSet);
    } finally {
      restoring = false;
    }

    // Write df.selection from the match (Assumption A5: no onFilterChanged re-trigger).
    df.selection.setAll(false, false);
    for (const idx of matched) df.selection.set(idx, true, false);
    df.selection.fireChanged();

    // Pitfall 7: the top-N labels must coexist with the search match. Labels
    // are decoupled from selection now, so label the top-N PLUS the matched
    // rows (the matched rows also stay highlighted via the selection above).
    const sp = findVolcano();
    if (sp) applyTopNLabels(df, sp, getVolcanoTopN(df), matched);
  }));

  return true;
}

export class PackageFunctions {
  @grok.decorators.init({tags: ['init']})
  static async initProteomics(): Promise<void> {
  }

  /**
   * Adds the dynamic Proteomics analysis groups (Annotate / Analyze / Visualize / Share)
   * to a table view's ribbon — see `src/menu.ts`. Import is NOT here; it is a decorator
   * `top-menu` (above) so the platform keeps it in the main left-bar menu and on every
   * ribbon. These groups need `isEnabled` grey-out, which only works on a programmatic
   * menu, so we build them onto the SAME ribbon 'Proteomics' group (get-or-create) the
   * decorator Import created — but only for tables that actually hold Proteomics data.
   *
   * Built once per view (tracked in `done`): `isEnabled` re-evaluates on every menu open,
   * so there is no need to rebuild on view activation — which is what caused the earlier
   * duplicate-menu bug. A view whose table isn't yet a Proteomics table (builder returns
   * false) is left untracked so a later activation retries.
   */
  @grok.decorators.func({tags: ['autostart'], meta: {autostartImmediate: 'true'}})
  static async buildProteomicsTopMenu(): Promise<void> {
    const handlers = {
      importSpectronautCandidates: () => PackageFunctions.importSpectronautCandidates(),
      importSpectronaut: () => PackageFunctions.importSpectronaut(),
      importMaxQuant: () => PackageFunctions.importMaxQuant(),
      importFragPipe: () => PackageFunctions.importFragPipe(),
      importGenericMatrix: () => PackageFunctions.importGenericMatrix(),
      annotateExperiment: () => PackageFunctions.annotateExperiment(),
      setLog2Scale: () => PackageFunctions.setLog2Scale(),
      normalize: () => PackageFunctions.normalizeProteomics(),
      impute: () => PackageFunctions.imputeMissingValues(),
      differentialExpression: () => PackageFunctions.differentialExpression(),
      enrichmentAnalysis: () => PackageFunctions.enrichmentAnalysis(),
      exportEnrichmentInputs: () => PackageFunctions.exportEnrichmentInputs(),
      computeSpcStatus: () => PackageFunctions.computeSpcStatus(),
      showVolcanoPlot: () => PackageFunctions.showVolcanoPlot(),
      showHeatmap: () => PackageFunctions.showHeatmap(),
      showPcaPlot: () => PackageFunctions.showPcaPlot(),
      showGroupMeanCorrelation: () => PackageFunctions.showGroupMeanCorrelation(),
      showQcDashboard: () => PackageFunctions.showQcDashboard(),
      showSpcDashboard: () => PackageFunctions.showSpcDashboard(),
      enrichmentCharts: () => PackageFunctions.enrichmentCharts(),
      showAllVisualizations: () => PackageFunctions.showAllVisualizations(),
      shareAnalysisForReview: () => PackageFunctions.shareAnalysisForReview(),
    };

    const done = new WeakSet<object>();
    const ensure = (view: any): void => {
      if (view?.type !== DG.VIEW_TYPE.TABLE_VIEW || done.has(view)) return;
      // Defer so the platform has built the decorator Import into the ribbon's Proteomics
      // group first; our get-or-create then attaches to it instead of racing a duplicate.
      setTimeout(() => {
        if (done.has(view)) return;
        try {
          if (buildProteomicsRibbonMenu(view as DG.TableView, handlers))
            done.add(view);
        } catch { /* ribbon not ready — a later activation retries */ }
      }, 50);
    };

    try { for (const v of grok.shell.views) ensure(v); } catch { /* no views yet */ }
    grok.events.onViewAdded.subscribe(ensure);
    grok.events.onCurrentViewChanged.subscribe(() => ensure(grok.shell.v));
  }

  // Import lives on decorator `top-menu` (not the dynamic ribbon builder) so it shows in
  // the main Datagrok menu (left bar) as an always-available entry point. Declaration
  // order sets menu order: Candidates first, then Report, then the rest.
  @grok.decorators.func({'top-menu': 'Proteomics | Import | Spectronaut Candidates...'})
  static async importSpectronautCandidates(): Promise<void> {
    DG.Utils.openFile({
      accept: '.tsv,.txt,.csv',
      open: async (file: File) => {
        try {
          const text = await file.text();
          const df = await parseSpectronautCandidatesText(text);
          df.name = file.name.replace(/\.[^.]+$/, '');
          autoDetectOrganism(df);
          const tv = grok.shell.addTableView(df);
          grok.shell.info(`Imported ${df.rowCount} candidates from Spectronaut`);
          dockComparisonFilterIfMultiContrast(tv, df);
          focusProtein(df);
        } catch (e: any) {
          grok.shell.error(`Failed to import Spectronaut Candidates file: ${e?.message ?? e}`);
        }
      },
    });
  }

  @grok.decorators.func({'top-menu': 'Proteomics | Import | Spectronaut Report...'})
  static async importSpectronaut(): Promise<void> {
    DG.Utils.openFile({
      accept: '.tsv,.txt,.csv',
      open: async (file: File) => {
        try {
          const df = (await sniffIsPrecursor(file)) ?
            await parseSpectronautStream(file) : await parseSpectronautText(await file.text());
          df.name = file.name.replace(/\.[^.]+$/, '');
          autoDetectOrganism(df);
          grok.shell.addTableView(df);
          grok.shell.info(`Imported ${df.rowCount} protein groups from Spectronaut`);
          focusProtein(df);
        } catch (e: any) {
          grok.shell.error(`Failed to import Spectronaut file: ${e?.message ?? e}`);
          grok.shell.warning(
            `If this is a very large precursor-level Spectronaut report that fails to ` +
            `stream in-browser, pre-aggregate it to the importable protein-group shape ` +
            `with the bundled duckdb fallback: tools/spectronaut-aggregate.sh <input.tsv>, ` +
            `then re-import the resulting file. See files/demo/README.md for usage.`,
          );
        }
      },
    });
  }

  @grok.decorators.func({'top-menu': 'Proteomics | Import | MaxQuant...'})
  static async importMaxQuant(): Promise<void> {
    DG.Utils.openFile({
      accept: '.txt,.tsv',
      open: async (file: File) => {
        try {
          const text = await file.text();
          const df = await parseMaxQuantText(text);
          df.name = file.name.replace(/\.[^.]+$/, '');
          autoDetectOrganism(df);
          grok.shell.addTableView(df);
          grok.shell.info(`Imported ${df.rowCount} protein groups`);
          focusProtein(df);
        } catch (e: any) {
          grok.shell.error(`Failed to import MaxQuant file: ${e?.message ?? e}`);
        }
      },
    });
  }

  @grok.decorators.func({'top-menu': 'Proteomics | Import | FragPipe...'})
  static async importFragPipe(): Promise<void> {
    DG.Utils.openFile({
      accept: '.tsv,.txt',
      open: async (file: File) => {
        try {
          const text = await file.text();
          const df = await parseFragPipeText(text);
          df.name = file.name.replace(/\.[^.]+$/, '');
          autoDetectOrganism(df);
          grok.shell.addTableView(df);
          grok.shell.info(`Imported ${df.rowCount} protein groups from FragPipe`);
          focusProtein(df);
        } catch (e: any) {
          grok.shell.error(`Failed to import FragPipe file: ${e?.message ?? e}`);
        }
      },
    });
  }

  @grok.decorators.func({'top-menu': 'Proteomics | Import | Generic Matrix...'})
  static async importGenericMatrix(): Promise<void> {
    showGenericImportDialog();
  }

  @grok.decorators.func()
  static async annotateExperiment(): Promise<void> {
    const df = grok.shell.tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Annotate Experiment')) return;
    showAnnotationDialog(df);
  }

  @grok.decorators.func()
  static async setLog2Scale(): Promise<void> {
    const df = grok.shell.tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Set Log2 Scale')) return;
    showLog2ScaleDialog(df);
  }

  @grok.decorators.func()
  static async normalizeProteomics(): Promise<void> {
    const df = grok.shell.tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Normalize')) return;
    showNormalizationDialog(df);
  }

  @grok.decorators.func()
  static async imputeMissingValues(): Promise<void> {
    const df = grok.shell.tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Impute Missing Values')) return;
    showImputationDialog(df);
  }

  @grok.decorators.func()
  static async differentialExpression(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Differential Expression')) return;
    showDEDialog(df, () => {
      // Auto-open volcano plot after DE completes. Title is synthesized by
      // createVolcanoPlot from proteomics.groups per the G1 contract.
      const sp = createVolcanoPlot(df);
      tv.addViewer(sp);
    });
  }

  @grok.decorators.func()
  static async computeSpcStatus(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) {
      grok.shell.warning('Open an analyzed file first, then run Compute SPC Status.');
      return;
    }
    if (df.getTag('proteomics.source') === 'spectronaut-candidates') {
      grok.shell.warning('SPC requires per-sample intensities. Re-import this analysis ' +
        'from the Spectronaut PG report (not the Candidates report) to compute SPC.');
      return;
    }
    const runMeta = getRunMeta(df);
    if (!runMeta || !runMeta.instrument_id || !runMeta.acquisition_datetime) {
      grok.shell.warning('Open Annotate Experiment to set instrument + acquisition datetime.');
      return;
    }
    const groups = getGroups(df);
    if (!groups) {
      grok.shell.warning('Annotate experimental groups first (Proteomics → Annotate Experiment) ' +
        '— SPC needs Group 1 to compute control-replicate correlation.');
      return;
    }
    if (groups.group1.columns.length < 2) {
      grok.shell.info("Group 1 has fewer than 2 samples — control-replicate correlation can't " +
        'be computed. The other three metrics will still be recorded.');
    }

    const pi = DG.TaskBarProgressIndicator.create('Computing SPC status...');
    try {
      const metrics = computeSpcMetrics(df, groups, runMeta);
      const baseline = await loadBaseline(runMeta.instrument_id);
      const priorRuns = await loadRuns(runMeta.instrument_id);
      const priorRow = priorRuns.find((r) =>
        r.acquisition_datetime === runMeta.acquisition_datetime);
      const priorStatus = priorRow?.status ?? null;

      let ruleResult: {status: 'pass' | 'flagged' | 'out_of_control'; rulesTripped: string[]};
      if (baseline === null) {
        ruleResult = {status: 'pass', rulesTripped: []};
        grok.shell.info(`No baseline locked for ${runMeta.instrument_id} yet — recorded as pass. ` +
          'Define a baseline once at least 4 runs have been computed.');
      } else {
        const priorSeries = {
          median_intensity: priorRuns.map((r) => r.median_intensity),
          missing_pct: priorRuns.map((r) => r.missing_pct),
          control_corr: priorRuns.map((r) => r.control_corr),
          protein_count: priorRuns.map((r) => r.protein_count),
        };
        ruleResult = evaluateNelsonRulesAllMetrics(
          metrics, priorSeries,
          baseline.metrics as BaselineSnapshot,
          baseline.rules_enabled ?? defaultRulesEnabledAllMetrics(),
        );
      }

      setSpcStatus(df, metrics, ruleResult);

      const row: Omit<SpcRunRow, 'run_id'> & {run_id?: string} = {
        instrument_id: runMeta.instrument_id,
        acquisition_datetime: runMeta.acquisition_datetime,
        run_label: df.name,
        median_intensity: metrics.median_intensity,
        missing_pct: metrics.missing_pct,
        control_corr: metrics.control_corr,
        protein_count: metrics.protein_count,
        status: ruleResult.status,
        rules_tripped: ruleResult.rulesTripped,
        source_project_id: null,
        source_df_name: df.name,
        computed_at: metrics.computed_at,
      };
      await upsertRun(row);

      if (priorStatus !== null)
        grok.shell.info(`Updated SPC for ${df.name} (previous status: ${priorStatus}).`);
      else if (ruleResult.status === 'pass')
        grok.shell.info(`SPC computed: ${df.name} — status: pass.`);
      else if (ruleResult.status === 'flagged')
        grok.shell.info(`SPC computed: ${df.name} — flagged on ${ruleResult.rulesTripped.length} ` +
          `rule(s): ${ruleResult.rulesTripped.join(', ')}.`);
      else
        grok.shell.info(`SPC computed: ${df.name} — OUT OF CONTROL on ` +
          `${ruleResult.rulesTripped.length} rule(s): ${ruleResult.rulesTripped.join(', ')}.`);
    } finally {
      pi.close();
    }
  }

  @grok.decorators.func()
  static async showVolcanoPlot(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireDifferentialExpression(df,
      'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;

    // Single menu item: open the options dialog. It reconfigures an existing
    // volcano, or — when none exists yet — CREATES it on OK (so nothing is
    // drawn until the user confirms; Cancel draws nothing). The old separate
    // 'Volcano Options...' item is folded in here.
    await PackageFunctions.volcanoOptions();
  }

  /** Volcano Plot options dialog. No longer its own menu item — showVolcanoPlot
   * delegates here. Reconfigures the existing volcano, or creates one on OK when
   * none exists. `spArg` is an already-resolved viewer (skips the lookup). */
  static async volcanoOptions(spArg?: DG.ScatterPlotViewer): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireDifferentialExpression(df,
      'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;

    // The volcano is the scatter plot whose Y is the -log10(metric) column.
    // It may NOT exist yet — in that case sp stays undefined and the OK handler
    // creates it, so the plot is only drawn once the user confirms.
    let sp: DG.ScatterPlotViewer | undefined = spArg;
    if (!sp) {
      for (const v of tv.viewers) {
        if (v.type === DG.VIEWER.SCATTER_PLOT &&
          (v as DG.ScatterPlotViewer).props.yColumnName === 'negLog10P') {
          sp = v as DG.ScatterPlotViewer;
          break;
        }
      }
    }

    // G2 dialog-state preload: snapshot current viewer state at open time and
    // seed every input from it. Pitfall 2 — getOptions is a snapshot, not a
    // live binding; OK uses the input values, never a re-read. With no existing
    // plot, seed from defaults.
    const state: {metric: MetricKind; colorDim: 'significance' | 'location'} = sp ?
      readVolcanoState(df, sp) : {metric: 'adj.p-value', colorDim: 'significance'};
    const hasPValue = df.col('p-value') != null;
    const initialMetric: MetricKind =
      hasPValue ? state.metric : 'adj.p-value';
    const metricInput = ui.input.choice('Significance metric', {
      value: initialMetric,
      items: hasPValue ? ['adj.p-value', 'p-value'] : ['adj.p-value'],
      nullable: false,
    });
    metricInput.setTooltip(hasPValue ?
      'Y axis, classification and threshold line all switch together' :
      'p-value column not present — only adj.p-value is available');
    const colorInput = ui.input.choice('Color by', {
      value: state.colorDim,
      items: ['significance', 'location'],
      nullable: false,
    });
    colorInput.setTooltip(
      'significance = Enriched in g1 / Enriched in g2 / Not significant; ' +
      'location = UniProt subcellular location');
    const labelTopNInput = ui.input.int('Label top N points', {
      value: getVolcanoTopN(df),
      min: 0,
    });
    labelTopNInput.setTooltip(
      'How many of the most significant proteins get name labels. 0 = none. ' +
      'Labels no longer touch your row selection.');

    // Optional axis-max overrides so volcanoes from different contrasts can be
    // pinned to a shared scale and compared side-by-side. Empty = auto-scale.
    const axis = getVolcanoAxisMax(df);
    const xMaxInput = ui.input.float('X-axis max (|log2FC|)', {
      value: axis.xMax ?? undefined,
      min: 0,
      nullable: true,
    });
    xMaxInput.setTooltip(
      'Pin the X axis to ±this |log2FC|. Leave empty to auto-scale. ' +
      'Set the same value on two volcanoes to compare them side-by-side.');
    const yMaxInput = ui.input.float('Y-axis max (−log10 p)', {
      value: axis.yMax ?? undefined,
      min: 0,
      nullable: true,
    });
    yMaxInput.setTooltip(
      'Pin the Y axis to 0…this −log10(p). Leave empty to auto-scale. ' +
      'Set the same value on two volcanoes to compare them side-by-side.');

    ui.dialog('Volcano Options')
      .add(metricInput)
      .add(colorInput)
      .add(labelTopNInput)
      .add(xMaxInput)
      .add(yMaxInput)
      .onOK(async () => {
        const metric = (metricInput.value ?? 'adj.p-value') as MetricKind;
        const colorDim = colorInput.value === 'location' ? 'location' : 'significance';

        // Persist the axis-max overrides before recompute — recomputeVolcano
        // ends with applyVolcanoAxisBounds, which reads these tags back.
        setVolcanoAxisMax(df, xMaxInput.value ?? null, yMaxInput.value ?? null);

        // First-time path: create the volcano NOW (on OK), so it isn't drawn
        // until the user confirms. topNLabels=0 here — applyTopNLabels below
        // sets the chosen count after the metric/color recompute.
        if (!sp) {
          sp = createVolcanoPlot(df, {topNLabels: 0});
          tv.addViewer(sp);
        }

        // Cache-aware pre-OK toast — the column-present path is sub-second
        // (short-circuit in ensureLocationColumn) so no toast is needed; the
        // column-absent + cold-cache path takes 30-90 s on a real Spectronaut
        // file and the user deserves a heads-up; the column-absent + warm-cache
        // path is a few seconds at most and deserves a shorter toast.
        if (colorDim === 'location' && !df.col(LOCATION_COL)) {
          const cache = await grok.dapi.userDataStorage.get(SUBCELL_STORE).catch(() => null);
          const cacheEntryCount = cache ?
            Object.keys(cache).filter((k) => k !== '__schema_v').length : 0;
          if (cacheEntryCount === 0) {
            grok.shell.info(
              'Fetching subcellular locations from UniProt — may take a minute or two ' +
              'on first use without cache; subsequent toggles use the cache.');
          } else {
            grok.shell.info('Resolving subcellular locations from cache.');
          }
        }

        // G3 wording: progress label is color-specific so a slow Color =
        // Location click reads as classification work, not generic "updating".
        const initialLabel = colorDim === 'location' ?
          'Classifying subcellular locations…' : 'Updating volcano…';
        const pi = DG.TaskBarProgressIndicator.create(initialLabel);

        // 13-10: also show progress on the volcano viewer itself — the
        // TaskBarProgressIndicator at the bottom of the platform shell is too
        // easy to miss while staring at a stale-looking chart. Attach when the
        // recompute can actually fetch from UniProt (colorDim === 'location').
        // On the 13-08 warm-cache short-circuit the overlay either flashes for
        // a single frame or appears very briefly; hideVolcanoBusy in finally
        // cleans up either way.
        const willFetchLocation = colorDim === 'location';
        if (willFetchLocation) showVolcanoBusy(sp!, initialLabel);

        // Map ProgressCb phases onto a human-readable pi.update label.
        const progress = (done: number, total: number, phase: string) => {
          const label = phase === 'fetch-acc' ? 'Fetching subcellular locations' :
            phase === 'fetch-gene' ? 'Resolving by gene name' :
              phase === 'init-column' ? 'Classifying subcellular locations' :
                colorDim === 'location' ? 'Classifying subcellular locations' :
                  'Updating volcano';
          const pct = total > 0 ? Math.round((done / total) * 100) : 0;
          pi.update(pct, `${label}: ${done}/${total}`);
          if (willFetchLocation) {
            const detail = total > 1 ? `${done}/${total} (${pct}%)` : '';
            updateVolcanoBusy(sp!, label, detail);
          }
        };
        try {
          await recomputeVolcano(df, sp!, metric, colorDim,
            DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD, progress);
          // Re-rank labels against the (possibly new) metric, applying the
          // chosen count. Decoupled from selection — leaves the user's rows be.
          applyTopNLabels(df, sp!, labelTopNInput.value ?? getVolcanoTopN(df));
        } catch (e: any) {
          grok.shell.error(`Volcano update failed: ${e?.message ?? e}`);
        } finally {
          pi.close();
          if (willFetchLocation && sp) hideVolcanoBusy(sp);
        }
      })
      .show();
  }

  @grok.decorators.func()
  static async showHeatmap(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Heatmap')) return;
    if (!requireDifferentialExpression(df,
      'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;

    // M4 — top-N is now user-configurable (was hardcoded 50). Small dialog so the
    // analyst can widen/narrow the heatmap; the ellipsis menu item implies it.
    const topNInput = ui.input.int('Top N proteins', {value: 50, min: 1});
    topNInput.setTooltip('How many of the most significant proteins (by adj. p-value) to show.');
    ui.dialog('Heatmap Options')
      .add(topNInput)
      .onOK(async () => {
        const topN = topNInput.value ?? 50;
        const pi = DG.TaskBarProgressIndicator.create('Creating heatmap...');
        try {
          const grid = await createExpressionHeatmap(df,
            {topN, title: `Heatmap: Top ${topN} DE Proteins`});
          tv.addViewer(grid);
        } finally {
          pi.close();
        }
      })
      .show();
  }

  @grok.decorators.func()
  static async showPcaPlot(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'PCA')) return;
    // PCA does NOT require DE -- it works after import/normalization
    const groups = getGroups(df);
    if (!groups) {
      grok.shell.warning('Annotate experimental groups first (Proteomics | Annotate Experiment)');
      return;
    }
    const allCols = [...groups.group1.columns, ...groups.group2.columns];
    const {viewer: sp, pcaDf} = createPcaPlot(df, allCols, groups, 'PCA: All Groups');
    pcaDf.name = `PCA: ${df.name}`;
    // PCA creates its own sample-level DataFrame -- open in separate table view
    const pcaTv = grok.shell.addTableView(pcaDf);
    pcaTv.addViewer(sp);
  }

  @grok.decorators.func()
  static async showGroupMeanCorrelation(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireDifferentialExpression(df,
      'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;
    const groups = getGroups(df);
    if (!groups) {
      grok.shell.warning('Annotate experimental groups first (Proteomics | Annotate Experiment)');
      return;
    }
    const sp = createGroupMeanCorrelation(df);
    tv.addViewer(sp);
  }

  @grok.decorators.func()
  static async showQcDashboard(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'QC Dashboard')) return;
    openQcDashboard(df);
  }

  @grok.decorators.func()
  static async showSpcDashboard(): Promise<void> {
    await openSpcDashboard();
  }

  @grok.decorators.func()
  static async showAllVisualizations(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!tv || !df) { grok.shell.warning('No table open'); return; }
    if (!requireSampleLevelData(df, 'Show All Visualizations')) return;
    if (!requireDifferentialExpression(df, 'Run Differential Expression first')) return;
    // Volcano -- docked in current table view (protein-level)
    const groups = getGroups(df);
    const volcano = createVolcanoPlot(df);
    tv.addViewer(volcano);
    // Heatmap -- docked in current table view (protein-level)
    const heatmapPi = DG.TaskBarProgressIndicator.create('Creating heatmap...');
    try {
      const heatmap = await createExpressionHeatmap(df, {title: 'Heatmap: Top 50 DE Proteins'});
      tv.addViewer(heatmap);
    } finally {
      heatmapPi.close();
    }
    // PCA (if groups available) -- opens in SEPARATE table view (sample-level, different row count)
    if (groups) {
      const allCols = [...groups.group1.columns, ...groups.group2.columns];
      const {viewer: pcaViewer, pcaDf} = createPcaPlot(df, allCols, groups, 'PCA: All Groups');
      pcaDf.name = `PCA: ${df.name}`;
      const pcaTv = grok.shell.addTableView(pcaDf);
      pcaTv.addViewer(pcaViewer);
    }
  }

  @grok.decorators.func()
  static async enrichmentAnalysis(): Promise<void> {
    const df = grok.shell.tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    showEnrichmentDialog(df);
  }

  @grok.decorators.func()
  static async exportEnrichmentInputs(): Promise<void> {
    const df = grok.shell.tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    if (!requireDifferentialExpression(df,
      'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;
    showEnrichmentInputExportDialog(df);
  }

  @grok.decorators.func()
  static async enrichmentCharts(): Promise<void> {
    const tv = grok.shell.tv;
    const df = tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }
    if (df.getTag('proteomics.enrichment') !== 'true') {
      grok.shell.warning('Run Enrichment Analysis first (Proteomics | Enrichment Analysis)');
      return;
    }
    // Find protein DataFrame(s) by proteomics.de_complete tag. If multiple are open,
    // pick the first but warn — the user may have meant a different one.
    const candidates = grok.shell.tables.filter((t: DG.DataFrame) => t.getTag('proteomics.de_complete') === 'true');
    if (candidates.length === 0) {
      grok.shell.warning('No protein table with DE results found. Enrichment charts will open without volcano linking.');
      openEnrichmentVisualization(df, df);
      return;
    }
    if (candidates.length > 1)
      grok.shell.warning(`Multiple protein tables with DE results found; linking enrichment to "${candidates[0].name}"`);
    openEnrichmentVisualization(df, candidates[0]);
  }

  @grok.decorators.func({
    name: 'Proteomics Demo',
    description: 'End-to-end proteomics differential-expression analysis on the HYE benchmark: ' +
      'import → impute → differential expression → volcano + heatmap, with a QC dashboard tab ' +
      'and the top hit pre-selected for its UniProt panel',
    meta: {demoPath: 'Proteomics | Differential Expression', isDemoDashboard: 'true'},
  })
  static async proteomicsDemo(): Promise<void> {
    await runProteomicsDemo();
  }

  @grok.decorators.func({
    name: 'Proteomics Enrichment Demo',
    description: 'Pathway enrichment (g:Profiler GO / KEGG / Reactome / WikiPathways) on a human ' +
      'differential-expression result — cell-cycle up, oxidative-phosphorylation down — with the ' +
      'enrichment charts cross-linked to the volcano',
    meta: {demoPath: 'Proteomics | Enrichment Analysis', isDemoDashboard: 'true'},
  })
  static async proteomicsEnrichmentDemo(): Promise<void> {
    await runEnrichmentDemo();
  }

  @grok.decorators.panel({
    name: 'Proteomics | UniProt',
    description: 'UniProt protein details',
    meta: {role: 'widgets'},
  })
  static uniprotPanelWidget(
    @grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string,
  ): DG.Widget {
    return uniprotPanel(proteinId);
  }

  @grok.decorators.func()
  static async shareAnalysisForReview(): Promise<void> {
    const tv = grok.shell.tv;
    let df = tv?.dataFrame;
    if (!df) { grok.shell.warning('No table open'); return; }

    // Share may be invoked from the Enrichment Results tab, whose df is the
    // enrichment frame (no de_complete) rather than the analyzed protein table.
    // Resolve the protein DataFrame: use the current df if it's DE-complete,
    // otherwise find an analyzed table among open tables (mirrors the
    // enrichment-charts handler) so the user doesn't have to switch tabs first.
    if (df.getTag('proteomics.de_complete') !== 'true') {
      const candidates = grok.shell.tables.filter(
        (t: DG.DataFrame) => t.getTag('proteomics.de_complete') === 'true');
      if (candidates.length === 0) {
        requireDifferentialExpression(df,
          'Run Differential Expression first (Proteomics | Analyze | Differential Expression)');
        return;
      }
      if (candidates.length > 1)
        grok.shell.warning(`Multiple analyzed tables open; sharing "${candidates[0].name}".`);
      df = candidates[0];
    }
    await showShareForReviewDialog(df);
  }

  @grok.decorators.panel({
    name: 'Proteomics | Shared Analysis',
    description: 'Audit context for a shared analysis snapshot',
    meta: {role: 'widgets'},
  })
  static publishedAnalysisPanelWidget(
    @grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string,
  ): DG.Widget {
    return publishedAnalysisPanel(proteinId);
  }

  @grok.decorators.func({tags: ['autostart'], meta: {autostartImmediate: 'true'}})
  static async recoverPublishedProjectsOnStartup(): Promise<void> {
    const tryEvent: any =
      (grok.events as any).onProjectOpened ??
      (grok.events as any).onCurrentProjectChanged ??
      null;
    const handleOpen = async (): Promise<void> => {
      await new Promise((r) => setTimeout(r, 200));
      const tables = (grok.shell as any).tables as DG.DataFrame[] | undefined;
      if (!Array.isArray(tables)) return;
      for (const df of tables) {
        if (isPublished(df)) {
          try { await recoverPublishedProject(df); } catch (e: any) {
            grok.shell.warning(`Could not auto-recover shared analysis on open: ${e?.message ?? e}`);
          }
        }
      }
    };
    if (tryEvent && typeof tryEvent.subscribe === 'function') {
      tryEvent.subscribe(() => { void handleOpen(); });
    } else {
      const onViewAdded: any = (grok.events as any).onViewAdded;
      if (onViewAdded && typeof onViewAdded.subscribe === 'function') {
        onViewAdded.subscribe((view: DG.View) => {
          if ((view as any) instanceof DG.TableView) {
            const df = (view as any).dataFrame;
            if (df && isPublished(df)) {
              recoverPublishedProject(df).catch((e: any) => {
                grok.shell.warning(`Could not auto-recover shared analysis on open: ${e?.message ?? e}`);
              });
            }
          }
        });
      }
    }
  }
}
