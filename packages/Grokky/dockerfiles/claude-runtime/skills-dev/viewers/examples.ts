/**
 * Compile-time verification of every fenced example in SKILL.md.
 *
 * Every `datagrok-exec` block referenced from SKILL.md appears here as a
 * named function with the same code. If any example uses a non-existent
 * method or wrong arg shape, this file fails to type-check.
 *
 * Each example is wrapped in a function with the same global signature the
 * runtime injects into a `datagrok-exec` block: `(grok, ui, DG, view, t)`.
 * No `grokky` parameter — the viewers skill is now pure raw js-api.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

type ExampleFn = (
  grok: typeof import('datagrok-api/grok'),
  ui: typeof import('datagrok-api/ui'),
  DG: typeof import('datagrok-api/dg'),
  view: DG.TableView,
  t: DG.DataFrame,
) => Promise<unknown> | unknown;

export const examples: Array<{name: string; fn: ExampleFn}> = [
  // -----------------------------------------------------------------------
  // Inspecting a viewer's property schema
  // -----------------------------------------------------------------------

  {name: 'inspect-scatter-properties', fn: async (_grok, ui, DG, view, _t) => {
    const v = view.addViewer(DG.VIEWER.SCATTER_PLOT);
    const names = v.getProperties().map((p) => p.name);
    return ui.divText(names.join(', '));
  }},

  // -----------------------------------------------------------------------
  // Adding viewers — primary forms
  // -----------------------------------------------------------------------

  {name: 'add-scatter-mw-logp-color-activity', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.SCATTER_PLOT, {
      xColumnName: 'MW',
      yColumnName: 'LogP',
      colorColumnName: 'activity',
      showRegressionLine: true,
    });
  }},

  {name: 'add-histogram-split', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.HISTOGRAM, {
      valueColumnName: 'activity',
      splitColumnName: 'compound_class',
    });
  }},

  {name: 'add-line-chart-multi-y', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.LINE_CHART, {
      xColumnName: 'date',
      yColumnNames: ['revenue', 'expenses', 'profit'],
    });
  }},

  {name: 'add-bar-chart', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.BAR_CHART, {
      splitColumnName: 'category',
      valueColumnName: 'count',
      valueAggrType: 'sum',
    });
  }},

  {name: 'add-box-plot', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.BOX_PLOT, {
      valueColumnName: 'activity',
      categoryColumnName: 'compound_class',
    });
  }},

  {name: 'add-correlation-plot', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.CORR_PLOT, {
      columnNames: ['MW', 'cLogP', 'pIC50', 'PSA'],
    });
  }},

  // -----------------------------------------------------------------------
  // Wrong-vs-right view scoping
  // -----------------------------------------------------------------------

  {name: 'wrong-shell-tv', fn: async (grok, _ui, DG, _view, _t) => {
    // Wrong: hits the globally active TableView, which may not be the one
    // the user was looking at when the prompt fired.
    const tv = grok.shell.tv;
    tv.addViewer(DG.VIEWER.SCATTER_PLOT);
  }},

  {name: 'right-view-in-scope', fn: async (_grok, _ui, DG, view, _t) => {
    view.addViewer(DG.VIEWER.SCATTER_PLOT, {xColumnName: 'mw', yColumnName: 'logp'});
  }},

  // -----------------------------------------------------------------------
  // Configuring an existing viewer
  // -----------------------------------------------------------------------

  {name: 'configure-recolor', fn: async (_grok, _ui, DG, view, _t) => {
    const sp = Array.from(view.viewers).slice(1)
      .find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
    if (sp)
      sp.setOptions({colorColumnName: 'logP'});
  }},

  {name: 'configure-logarithmic', fn: async (_grok, _ui, DG, view, _t) => {
    const h = Array.from(view.viewers).slice(1)
      .find((v) => v.type === DG.VIEWER.HISTOGRAM);
    if (h)
      h.setOptions({yAxisType: 'logarithmic'});
  }},

  {name: 'create-then-configure', fn: async (_grok, _ui, DG, view, _t) => {
    const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT, {
      xColumnName: 'height',
      yColumnName: 'weight',
      showRegressionLine: true,
    });
    sp.setOptions({xMin: 150, xMax: 200, colorColumnName: 'age'});
  }},

  // -----------------------------------------------------------------------
  // Finding viewers
  // -----------------------------------------------------------------------

  {name: 'find-count-non-grid', fn: async (_grok, ui, _DG, view, _t) => {
    return ui.tableFromMap({count: Array.from(view.viewers).slice(1).length});
  }},

  {name: 'find-first-scatter', fn: async (_grok, _ui, DG, view, _t) => {
    const sp = Array.from(view.viewers).slice(1)
      .find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
    return sp ?? null;
  }},

  {name: 'find-all-histograms', fn: async (_grok, _ui, DG, view, _t) => {
    return Array.from(view.viewers).slice(1)
      .filter((v) => v.type === DG.VIEWER.HISTOGRAM);
  }},

  // -----------------------------------------------------------------------
  // Closing viewers
  // -----------------------------------------------------------------------

  {name: 'close-by-handle', fn: async (_grok, _ui, DG, view, _t) => {
    const sp = Array.from(view.viewers).slice(1)
      .find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
    if (sp) {
      try { sp.close(); } catch (e) { console.warn(e); }
    }
  }},

  {name: 'close-by-type', fn: async (_grok, _ui, DG, view, _t) => {
    Array.from(view.viewers).slice(1)
      .filter((v) => v.type === DG.VIEWER.SCATTER_PLOT)
      .forEach((v) => { try { v.close(); } catch { /* swallow */ } });
  }},

  {name: 'close-all-non-grid', fn: async (_grok, _ui, _DG, view, _t) => {
    Array.from(view.viewers).slice(1)
      .forEach((v) => { try { v.close(); } catch { /* swallow */ } });
  }},

  // -----------------------------------------------------------------------
  // Close-and-replace idioms
  // -----------------------------------------------------------------------

  {name: 'replace-scatter-with-histogram', fn: async (_grok, _ui, DG, view, _t) => {
    Array.from(view.viewers).slice(1)
      .filter((v) => v.type === DG.VIEWER.SCATTER_PLOT)
      .forEach((v) => { try { v.close(); } catch { /* swallow */ } });
    view.addViewer(DG.VIEWER.HISTOGRAM, {valueColumnName: 'activity'});
  }},

  // -----------------------------------------------------------------------
  // Plugin viewer (Chem:chemSpaceTopMenu) — function, not viewer type
  // -----------------------------------------------------------------------

  {name: 'chem-space-plugin', fn: async (grok, _ui, _DG, _view, t) => {
    const sp = await grok.functions.call('Chem:chemSpaceTopMenu', {
      table: t,
      molecules: t.col('smiles'),
      methodName: 'UMAP',
      similarityMetric: 'Tanimoto',
      plotEmbeddings: true,
    }) as DG.Viewer | undefined;
    if (sp)
      sp.setOptions({colorColumnName: 'activity'});
  }},

  // -----------------------------------------------------------------------
  // Standalone (detached) viewer
  // -----------------------------------------------------------------------

  {name: 'standalone-fromType', fn: async (_grok, _ui, DG, _view, t) => {
    const sp = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, t, {
      xColumnName: 'mw', yColumnName: 'logp',
    });
    return sp.root;
  }},

  // -----------------------------------------------------------------------
  // Dock layout — explicit position override
  // -----------------------------------------------------------------------

  {name: 'dock-right-with-title', fn: async (_grok, _ui, DG, view, _t) => {
    const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT, {xColumnName: 'mw', yColumnName: 'logp'});
    view.dockManager.dock(sp, 'right', null, 'Structure–activity', 0.4);
  }},

  // Keep the imports referenced even when no example uses them at module scope.
  {name: '_noop-imports', fn: () => {
    void grok; void ui; void DG;
  }},
];
