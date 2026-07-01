/**
 * Phase 16 SPC Dashboard.
 *
 * Sibling-TableView host for the runs.csv-derived DataFrame (per
 * ARCHITECTURE.md §"Cross-DataFrame Selection") — the dashboard NEVER docks
 * onto the protein DataFrame's view. Composes Plan 02's math + Plan 03's
 * AppData I/O into a live operator surface: instrument picker, 4 I-charts
 * + 4 MR-charts with UCL/CL/LCL formula lines, baseline status sidebar with
 * a per-metric Nelson rule toggle grid, Define / Rebuild baseline modal,
 * and drill-down to the source Project on flagged-point click.
 *
 * Pareto chart (SPC-08, P2) is Plan 16-07's deliverable — `aggregateParetoCounts`
 * lives there. The rules-tripped tag is set unconditionally so Plan 07 can
 * land as a hotfix without unwinding this plan.
 */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  WHEELER_D4,
  evaluateNelsonRulesAllMetrics,
  defaultRulesEnabledAllMetrics,
  classifyStatus,
  BaselineSnapshot, SpcMetrics, SpcRulesEnabled,
} from '../analysis/spc';
import {
  LockedBaseline, SpcRunRow,
  loadRuns, loadBaseline, saveBaseline,
  iterativeOutlierRemoval, computeBaselineStats,
  runsToDataFrame,
} from '../analysis/spc-storage';

type MetricKey = 'median_intensity' | 'missing_pct' | 'control_corr' | 'protein_count';

const METRIC_KEYS: MetricKey[] = [
  'median_intensity', 'missing_pct', 'control_corr', 'protein_count',
];

const METRIC_LABELS: Record<MetricKey, string> = {
  median_intensity: 'Median log2 intensity',
  missing_pct: 'Missing %',
  control_corr: 'Control replicate correlation',
  protein_count: 'Protein count',
};

const FORMULA_LINE_COLORS = {
  cl: '#888888',
  bound: '#666666',
};

const LAST_INSTRUMENT_STORAGE_KEY = 'spc.last_picked_instrument';

// Five canonical false-alarm-rate disclosures per UI-SPEC §"False-alarm-rate".
const FALSE_ALARM_DISCLOSURE: Record<string, string> = {
  'rule_1_only': 'Expected false-alarm rate at this rule set: ~0.27%. (Rule 1 only — a point beyond 3σ.)',
  'rules_1_5':   'Expected false-alarm rate at this rule set: ~0.4%. (Rules 1 + 5 — a single point beyond 3σ, or 2 of 3 beyond 2σ.)',
  'rules_1_5_2': 'Enabling rules 1+5+2 raises the expected false-alarm rate to ~1.1%.',
  'rules_all_8': 'Expected false-alarm rate at this rule set: ~2.65%. Enabling all 8 rules will flag roughly 1 in 38 normal-operation runs.',
  'rules_other': 'More rules enabled → higher false-alarm rate. The Nelson 1+5 default is calibrated for ~0.4%.',
};

function disclosureKey(enabled: SpcRulesEnabled): string {
  const sample = enabled.median_intensity;
  const onIds = Object.entries(sample).filter(([_, on]) => on).map(([id]) => id).sort();
  if (onIds.length === 1 && onIds[0] === 'nelson_1') return 'rule_1_only';
  if (onIds.length === 2 && onIds[0] === 'nelson_1' && onIds[1] === 'nelson_5') return 'rules_1_5';
  if (onIds.length === 3 && onIds.includes('nelson_1') && onIds.includes('nelson_5') && onIds.includes('nelson_2'))
    return 'rules_1_5_2';
  if (onIds.length === 8) return 'rules_all_8';
  return 'rules_other';
}

// =====================================================================================
// Pure helpers — exported because they're directly testable.
// =====================================================================================

export interface ControlLines {
  ucl: number;
  cl: number;
  lcl: number;
}

/** UCL = mean + 3sd, CL = mean, LCL = mean - 3sd. Used for the I-chart formula lines
 *  AND tested by SPC:formula_lines. Returns NaN-bearing lines when the baseline
 *  metric is unavailable; the caller skips formula-line addition in that case. */
export function computeControlLines(
  baseline: {metrics: BaselineSnapshot} | null, metric: MetricKey,
): ControlLines {
  if (!baseline) return {ucl: NaN, cl: NaN, lcl: NaN};
  const m = baseline.metrics[metric];
  if (!m || !Number.isFinite(m.mean) || !Number.isFinite(m.sd))
    return {ucl: NaN, cl: NaN, lcl: NaN};
  return {ucl: m.mean + 3 * m.sd, cl: m.mean, lcl: m.mean - 3 * m.sd};
}

export interface DrillDownInput {
  run_id: string;
  run_label: string;
  source_project_id: string | null;
  instrument_id?: string;
}

export type DrillDownResult =
  | {kind: 'opened'; id: string}
  | {kind: 'toast'; message: string};

/** Pure drill-down resolver. Caller injects a `projectsFind` that returns either
 *  a project (with `.open()`) or null. The handler classifies into:
 *   - 'opened' on success (project found + opened)
 *   - 'toast' for both the null-id case AND the find()-returns-null case
 *  Both toast messages are biologist-readable per UI-SPEC §"Drill-down". */
export async function resolveDrillDown(
  row: DrillDownInput,
  deps: {projectsFind: (id: string) => Promise<{open: () => any} | null>},
): Promise<DrillDownResult> {
  if (row.source_project_id == null) {
    return {
      kind: 'toast',
      message: `Source for '${row.run_label}' is not currently available. ` +
        `Open the run's analyzed file and re-run Compute SPC Status.`,
    };
  }
  const project = await deps.projectsFind(row.source_project_id);
  if (!project) {
    return {
      kind: 'toast',
      message: `Source for '${row.run_label}' is not currently available. ` +
        `Open the run's analyzed file and re-run Compute SPC Status.`,
    };
  }
  await project.open();
  return {kind: 'opened', id: row.source_project_id};
}

// =====================================================================================
// createSpcChartPanel — pure factory returning {iChart, mrChart}. The MR column
// must already exist on runsDf (openSpcDashboard adds it before docking).
// =====================================================================================

export function createSpcChartPanel(
  runsDf: DG.DataFrame,
  baseline: {metrics: BaselineSnapshot} | null,
  metric: MetricKey,
): {iChart: DG.Viewer; mrChart: DG.Viewer} {
  const iChart: any = DG.Viewer.lineChart(runsDf, {
    xColumnName: 'acquisition_datetime',
    yColumnNames: [metric],
    title: METRIC_LABELS[metric],
  } as any);

  if (baseline) {
    const lines = computeControlLines(baseline, metric);
    if (Number.isFinite(lines.cl)) {
      const fl: any = iChart?.meta?.formulaLines;
      if (fl && typeof fl.addLine === 'function') {
        // Spike outcome (Plan 16-01) lands at runtime — both branches are
        // documented in the spike test; using `style:'dashed'` is the
        // expected PASS path. FAIL path is the constant-column overlay
        // fallback described in 16-RESEARCH Pattern 4.
        try {
          fl.addLine({
            formula: `\${${metric}} = ${lines.ucl}`,
            color: FORMULA_LINE_COLORS.bound, width: 1, style: 'dashed',
            visible: true, title: 'UCL',
          });
          fl.addLine({
            formula: `\${${metric}} = ${lines.cl}`,
            color: FORMULA_LINE_COLORS.cl, width: 1, style: 'solid',
            visible: true, title: 'CL',
          });
          fl.addLine({
            formula: `\${${metric}} = ${lines.lcl}`,
            color: FORMULA_LINE_COLORS.bound, width: 1, style: 'dashed',
            visible: true, title: 'LCL',
          });
        } catch (_e) {
          // Spike FAIL path — formula-line API didn't honor the style. The
          // dashboard still renders; the constant-overlay fallback is
          // deferred to a follow-up if the spike returns FAIL in practice.
        }
      }
    }
  }

  const mrColName = `MR_${metric}`;
  const mrChart: any = DG.Viewer.lineChart(runsDf, {
    xColumnName: 'acquisition_datetime',
    yColumnNames: [mrColName],
    title: `MR ${METRIC_LABELS[metric]}`,
  } as any);

  if (baseline && runsDf.col(mrColName)) {
    const mrCol = runsDf.col(mrColName)!;
    let sum = 0;
    let n = 0;
    for (let i = 0; i < runsDf.rowCount; i++) {
      if (mrCol.isNone(i)) continue;
      const v = Number(mrCol.get(i));
      if (Number.isFinite(v)) { sum += v; n++; }
    }
    if (n > 0) {
      const mrMean = sum / n;
      const mrUcl = WHEELER_D4 * mrMean;
      const flMr: any = (mrChart as any)?.meta?.formulaLines;
      if (flMr && typeof flMr.addLine === 'function') {
        try {
          flMr.addLine({
            formula: `\${${mrColName}} = ${mrUcl}`,
            color: FORMULA_LINE_COLORS.bound, width: 1, style: 'dashed',
            visible: true, title: 'MR UCL',
          });
          flMr.addLine({
            formula: `\${${mrColName}} = ${mrMean}`,
            color: FORMULA_LINE_COLORS.cl, width: 1, style: 'solid',
            visible: true, title: 'MR CL',
          });
        } catch (_e) { /* spike FAIL fallback — same as I-chart */ }
      }
    }
  }

  return {iChart, mrChart};
}

// =====================================================================================
// Pareto aggregation (SPC-08, P2) — Plan 16-07.
// =====================================================================================

export interface ParetoEntry {
  rule: string;
  metric: string;
  count: number;
}

/** Aggregates rules-tripped trip counts across a slice of runs and returns the
 *  entries sorted DESCENDING by count. Accepts an array shape (test-friendly:
 *  rules_tripped as a string[]) — the dashboard caller adapts from CSV strings
 *  via JSON.parse-with-try/catch before passing in. */
export function aggregateParetoCounts(
  runs: Array<{rules_tripped: string[]}>,
): ParetoEntry[] {
  const counts = new Map<string, number>();
  for (const r of runs) {
    const tripped = Array.isArray(r.rules_tripped) ? r.rules_tripped : [];
    for (const t of tripped) {
      const prev = counts.get(t) ?? 0;
      counts.set(t, prev + 1);
    }
  }
  const entries: ParetoEntry[] = [];
  for (const [key, count] of counts.entries()) {
    const at = key.indexOf('@');
    if (at < 0) entries.push({rule: key, metric: '', count});
    else entries.push({rule: key.slice(0, at), metric: key.slice(at + 1), count});
  }
  entries.sort((a, b) => b.count - a.count);
  return entries;
}

/** Returns a barChart viewer ordered DESC by trip count, or null when there
 *  are fewer than 4 runs (caller renders the empty-state body in that case). */
export function createParetoChart(
  runsDf: DG.DataFrame, instrumentId: string,
): DG.Viewer | null {
  if (runsDf.rowCount < 4) return null;
  const rulesCol = runsDf.col('rules_tripped');
  if (!rulesCol) return null;
  const totalRuns = runsDf.rowCount;
  const inputs: Array<{rules_tripped: string[]}> = [];
  for (let i = 0; i < totalRuns; i++) {
    const raw = String(rulesCol.get(i) ?? '');
    let parsed: string[] = [];
    if (raw.length > 0) {
      try { parsed = JSON.parse(raw) as string[]; } catch { parsed = []; }
    }
    inputs.push({rules_tripped: parsed});
  }
  const entries = aggregateParetoCounts(inputs);
  if (entries.length === 0) return null;

  const ruleMetricCol = DG.Column.fromStrings('rule_metric',
    entries.map((e) => e.metric ? `${e.rule}@${e.metric}` : e.rule));
  const tripCountArr = new Int32Array(entries.map((e) => e.count));
  const totalRunsArr = new Int32Array(entries.map(() => totalRuns));
  const pctArr = new Float32Array(entries.map((e) => (e.count / totalRuns) * 100));
  const aggDf = DG.DataFrame.fromColumns([
    ruleMetricCol,
    DG.Column.fromInt32Array('trip_count', tripCountArr),
    DG.Column.fromInt32Array('total_runs', totalRunsArr),
    DG.Column.fromFloat32Array('pct', pctArr),
  ]);

  const chart: any = DG.Viewer.barChart(aggDf, {
    splitColumnName: 'rule_metric',
    valueColumnName: 'trip_count',
    title: `Which rules trip most often — ${instrumentId}`,
  } as any);
  // Axis labels per UI-SPEC §"Pareto chart labels".
  try {
    chart.props.xLabel = '(rule, metric)';
    chart.props.yLabel = 'Times tripped';
  } catch (_e) { /* prop names may shift across platform versions — non-fatal */ }
  return chart;
}

// =====================================================================================
// MR column derivation (4 sibling columns on the runs df, computed once at open).
// =====================================================================================

function addMrColumns(runsDf: DG.DataFrame): void {
  for (const metric of METRIC_KEYS) {
    const srcCol = runsDf.col(metric);
    if (!srcCol) continue;
    const name = `MR_${metric}`;
    if (runsDf.columns.contains(name)) runsDf.columns.remove(name);
    const mr = runsDf.columns.addNewFloat(name);
    const rc = runsDf.rowCount;
    let prev: number | null = null;
    for (let i = 0; i < rc; i++) {
      if (i === 0) { mr.set(i, null); prev = null; continue; }
      if (srcCol.isNone(i) || srcCol.isNone(i - 1)) { mr.set(i, null); continue; }
      const v = Number(srcCol.get(i));
      const p = Number(srcCol.get(i - 1));
      if (!Number.isFinite(v) || !Number.isFinite(p)) { mr.set(i, null); continue; }
      mr.set(i, Math.abs(v - p));
      prev = v; void prev;
    }
  }
}

// =====================================================================================
// Sidebar — instrument readout + baseline status + Nelson rule toggle grid
// =====================================================================================

interface SidebarRefs {
  el: HTMLElement;
  setDisclosure: (text: string) => void;
}

function buildSidebar(
  instrumentId: string,
  baseline: LockedBaseline | null,
  onRebuildClick: () => void,
  onRuleToggle: (metric: MetricKey, rule: keyof SpcRulesEnabled['median_intensity'], value: boolean) => Promise<void>,
): SidebarRefs {
  const heading = ui.h2(`Instrument: ${instrumentId}`);

  let statusBlock: HTMLElement;
  if (baseline) {
    statusBlock = ui.divV([
      ui.divText(`Baseline locked: ${baseline.locked_at}`),
      ui.divText(`${baseline.included_run_ids.length} runs included · ` +
        `${baseline.excluded_run_ids.length} runs excluded (3σ outlier).`),
      ui.button('Rebuild baseline...', onRebuildClick),
    ]);
  } else {
    statusBlock = ui.divV([
      ui.divText('No baseline locked yet. The dashboard turns on once a baseline exists.'),
    ]);
  }

  const ruleIds: Array<keyof SpcRulesEnabled['median_intensity']> = [
    'nelson_1', 'nelson_2', 'nelson_3', 'nelson_4',
    'nelson_5', 'nelson_6', 'nelson_7', 'nelson_8',
  ];

  const enabled = baseline?.rules_enabled ?? defaultRulesEnabledAllMetrics();
  const disclosure = ui.divText(FALSE_ALARM_DISCLOSURE[disclosureKey(enabled)]);
  const setDisclosure = (text: string) => { disclosure.textContent = text; };

  const gridRows: HTMLElement[] = [];
  gridRows.push(ui.divH([
    ui.divText(''),
    ...ruleIds.map((r) => ui.divText(r.replace('nelson_', 'R'))),
  ]));
  for (const metric of METRIC_KEYS) {
    const row: HTMLElement[] = [ui.divText(METRIC_LABELS[metric])];
    for (const rid of ruleIds) {
      const checkbox = ui.input.bool('', {
        value: enabled[metric][rid],
        tooltipText: ruleTooltip(rid),
      });
      checkbox.onChanged.subscribe(async (v: any) => {
        await onRuleToggle(metric, rid, !!v);
      });
      row.push(checkbox.root);
    }
    gridRows.push(ui.divH(row));
  }

  const el = ui.divV([
    heading,
    statusBlock,
    ui.h2('Nelson rules'),
    ui.divV(gridRows),
    disclosure,
  ]);
  return {el, setDisclosure};
}

function ruleTooltip(rule: keyof SpcRulesEnabled['median_intensity']): string {
  switch (rule) {
  case 'nelson_1': return 'Nelson rule 1: a single point beyond 3σ from the centerline.';
  case 'nelson_2': return 'Nelson rule 2: 9 consecutive points on the same side of the centerline.';
  case 'nelson_3': return 'Nelson rule 3: 6 consecutive points all trending up, or all trending down.';
  case 'nelson_4': return 'Nelson rule 4: 14 consecutive points alternating above and below the centerline.';
  case 'nelson_5': return 'Nelson rule 5: 2 out of 3 consecutive points beyond 2σ on the same side.';
  case 'nelson_6': return 'Nelson rule 6: 4 out of 5 consecutive points beyond 1σ on the same side.';
  case 'nelson_7': return 'Nelson rule 7: 15 consecutive points within 1σ of the centerline (unusually stable).';
  case 'nelson_8': return 'Nelson rule 8: 8 consecutive points more than 1σ from the centerline on either side.';
  }
}

// =====================================================================================
// Drill-down subscription — wires the resolveDrillDown helper to a live TableView.
// =====================================================================================

function wireDrillDown(tv: DG.TableView, df: DG.DataFrame): void {
  tv.dataFrame.onCurrentRowChanged.subscribe(async () => {
    const i = df.currentRowIdx;
    if (i < 0) return;
    const status = String(df.col('status')?.get(i) ?? '');
    if (status !== 'flagged' && status !== 'out_of_control') return;
    const sourceRaw = df.col('source_project_id')?.get(i);
    const source = (sourceRaw == null || sourceRaw === '') ? null : String(sourceRaw);
    const result = await resolveDrillDown(
      {
        run_id: String(df.col('run_id')?.get(i) ?? ''),
        run_label: String(df.col('run_label')?.get(i) ?? ''),
        source_project_id: source,
      },
      {projectsFind: async (id: string) => {
        const p: any = await grok.dapi.projects.find(id);
        return p ? {open: () => p.open()} : null;
      }},
    );
    if (result.kind === 'toast') grok.shell.warning(result.message);
  });
}

// =====================================================================================
// Status recomputation — invoked on every rule toggle to update point colors live.
// =====================================================================================

function recomputeStatusesInPlace(
  runsDf: DG.DataFrame, baseline: LockedBaseline,
): void {
  const rowCount = runsDf.rowCount;
  if (rowCount === 0) return;
  const cols = {
    median_intensity: runsDf.col('median_intensity'),
    missing_pct: runsDf.col('missing_pct'),
    control_corr: runsDf.col('control_corr'),
    protein_count: runsDf.col('protein_count'),
  };
  const statusCol = runsDf.col('status');
  const rulesCol = runsDf.col('rules_tripped');
  if (!statusCol || !rulesCol) return;
  for (let i = 0; i < rowCount; i++) {
    const metrics: SpcMetrics = {
      median_intensity: Number(cols.median_intensity?.get(i)),
      missing_pct: Number(cols.missing_pct?.get(i)),
      control_corr: Number(cols.control_corr?.get(i)),
      protein_count: Number(cols.protein_count?.get(i)),
      sample_count: 0, computed_at: '',
    };
    const priorSeries = {
      median_intensity: [] as number[],
      missing_pct: [] as number[],
      control_corr: [] as number[],
      protein_count: [] as number[],
    };
    for (let j = 0; j < i; j++) {
      priorSeries.median_intensity.push(Number(cols.median_intensity?.get(j)));
      priorSeries.missing_pct.push(Number(cols.missing_pct?.get(j)));
      priorSeries.control_corr.push(Number(cols.control_corr?.get(j)));
      priorSeries.protein_count.push(Number(cols.protein_count?.get(j)));
    }
    const r = evaluateNelsonRulesAllMetrics(
      metrics, priorSeries, baseline.metrics, baseline.rules_enabled,
    );
    statusCol.set(i, r.status);
    rulesCol.set(i, JSON.stringify(r.rulesTripped));
  }
  void classifyStatus;
}

// =====================================================================================
// Define / Rebuild baseline modal
// =====================================================================================

export function showDefineBaselineDialog(
  instrumentId: string,
  runs: SpcRunRow[],
  onLocked: (baseline: LockedBaseline) => void | Promise<void>,
  existing?: LockedBaseline,
): void {
  const isRebuild = !!existing;
  const title = isRebuild
    ? `Rebuild baseline for ${instrumentId}`
    : `Lock baseline for ${instrumentId}`;

  const ordered = runs.slice().sort((a, b) =>
    a.acquisition_datetime.localeCompare(b.acquisition_datetime));

  // Default check selection per CONTEXT.md D-04.
  const checkedSet = new Set<string>();
  if (isRebuild) {
    for (const id of existing!.included_run_ids) checkedSet.add(id);
  } else {
    const n = Math.min(7, ordered.length);
    for (let i = 0; i < n; i++) checkedSet.add(ordered[i].run_id);
  }

  const runCheckboxes: Array<{run: SpcRunRow; input: DG.InputBase<boolean>}> = [];
  const listEl = ui.divV([]);
  listEl.style.maxHeight = '300px';
  listEl.style.overflowY = 'auto';
  for (const r of ordered) {
    const cb = ui.input.bool('', {value: checkedSet.has(r.run_id)});
    runCheckboxes.push({run: r, input: cb});
    listEl.appendChild(ui.divH([
      cb.root,
      ui.divText(r.acquisition_datetime.slice(0, 10)),
      ui.divText(r.run_label),
    ]));
  }

  const iterateToggle = ui.input.bool('Iterate 3σ outlier removal (recommended)', {
    value: isRebuild ? (existing!.iteration_trace.length > 0) : true,
    tooltipText: 'When on, runs more than 3 standard deviations from the baseline mean ' +
      'are dropped from the baseline. Capped at 2 passes.',
  });

  const initialEnabled: SpcRulesEnabled = existing?.rules_enabled ?? defaultRulesEnabledAllMetrics();
  const ruleIds: Array<keyof SpcRulesEnabled['median_intensity']> = [
    'nelson_1', 'nelson_2', 'nelson_3', 'nelson_4',
    'nelson_5', 'nelson_6', 'nelson_7', 'nelson_8',
  ];
  const ruleInputs: Record<MetricKey, Record<string, DG.InputBase<boolean>>> = {
    median_intensity: {} as any, missing_pct: {} as any,
    control_corr: {} as any, protein_count: {} as any,
  };
  const gridRows: HTMLElement[] = [
    ui.divH([ui.divText(''), ...ruleIds.map((r) => ui.divText(r.replace('nelson_', 'R')))]),
  ];
  const disclosureEl = ui.divText(FALSE_ALARM_DISCLOSURE[disclosureKey(initialEnabled)]);
  for (const metric of METRIC_KEYS) {
    const row: HTMLElement[] = [ui.divText(METRIC_LABELS[metric])];
    for (const rid of ruleIds) {
      const cb = ui.input.bool('', {
        value: initialEnabled[metric][rid],
        tooltipText: ruleTooltip(rid),
      });
      cb.onChanged.subscribe(() => {
        const snapshot = currentEnabled();
        disclosureEl.textContent = FALSE_ALARM_DISCLOSURE[disclosureKey(snapshot)];
      });
      ruleInputs[metric][rid] = cb;
      row.push(cb.root);
    }
    gridRows.push(ui.divH(row));
  }
  function currentEnabled(): SpcRulesEnabled {
    const out: any = {};
    for (const m of METRIC_KEYS) {
      out[m] = {} as any;
      for (const rid of ruleIds) out[m][rid] = !!ruleInputs[m][rid].value;
    }
    return out as SpcRulesEnabled;
  }

  const bodyEls: HTMLElement[] = [];
  if (isRebuild) {
    bodyEls.push(ui.divText(
      `This replaces the locked baseline from ${existing!.locked_at}. ` +
      'The previous baseline is not kept.',
    ));
  } else {
    bodyEls.push(ui.divText('Pick the runs that represent normal operation:'));
  }
  bodyEls.push(listEl);
  bodyEls.push(iterateToggle.root);
  bodyEls.push(ui.h2('Nelson rules'));
  bodyEls.push(ui.divV(gridRows));
  bodyEls.push(disclosureEl);

  const dialog = ui.dialog(title);
  for (const el of bodyEls) dialog.add(el);
  dialog.onOK(async () => {
    const checked = runCheckboxes.filter((c) => c.input.value === true).map((c) => c.run);
    if (checked.length < 3) {
      grok.shell.warning('At least 3 runs must remain after outlier removal — ' +
        'too few to compute a baseline.');
      return;
    }
    const exclusionUnion = new Set<string>();
    const traceUnion: Array<{iter: number; dropped: string[]}> = [];
    if (iterateToggle.value) {
      for (const metric of METRIC_KEYS) {
        const series = checked.map((r) => ({run_id: r.run_id, value: r[metric] as number}));
        const res = iterativeOutlierRemoval(series, 2);
        for (const id of res.excluded_run_ids) exclusionUnion.add(id);
        for (const t of res.iteration_trace) traceUnion.push(t);
      }
    }
    const retained = checked.filter((r) => !exclusionUnion.has(r.run_id));
    if (retained.length < 3) {
      grok.shell.warning('At least 3 runs must remain after outlier removal — ' +
        'too few to compute a baseline.');
      return;
    }
    const metricsSnapshot = computeBaselineStats(retained);
    const baseline: LockedBaseline = {
      instrument_id: instrumentId,
      locked_at: new Date().toISOString(),
      included_run_ids: retained.map((r) => r.run_id),
      excluded_run_ids: Array.from(exclusionUnion),
      iteration_trace: traceUnion,
      metrics: metricsSnapshot,
      rules_enabled: currentEnabled(),
    };
    await onLocked(baseline);
    if (isRebuild) {
      grok.shell.info(`Baseline rebuilt for ${instrumentId}: ${retained.length} runs included, ` +
        `${exclusionUnion.size} excluded. Previous baseline replaced.`);
    } else {
      grok.shell.info(`Baseline locked for ${instrumentId}: ${retained.length} runs included, ` +
        `${exclusionUnion.size} excluded.`);
    }
  });
  dialog.show();
}

// =====================================================================================
// openSpcDashboard — orchestrator
// =====================================================================================

export async function openSpcDashboard(): Promise<void> {
  const allRuns = await loadRuns();
  if (allRuns.length === 0) {
    grok.shell.warning('No SPC runs recorded yet. Open an analyzed file, then run ' +
      'Analyze → Compute SPC Status. The dashboard turns on once at least 4 runs have been computed.');
    return;
  }
  const instruments = Array.from(new Set(allRuns.map((r) => r.instrument_id)));
  // userDataStorage is soft-deprecated; matches v1.3 codebase (subcellular-location.ts:246).
  // Migrate to grok.userSettings as cross-cutting cleanup in v1.5.
  const lastPicked = await grok.dapi.userDataStorage.get(LAST_INSTRUMENT_STORAGE_KEY);
  const last = (lastPicked && lastPicked['instrument_id']) || instruments[0];
  const baselineByInstrument = new Map<string, LockedBaseline | null>();
  for (const inst of instruments)
    baselineByInstrument.set(inst, await loadBaseline(inst));

  const items = instruments.map((i) => {
    const b = baselineByInstrument.get(i);
    return b ? `${i} (baseline locked)` : `${i} (no baseline)`;
  });
  const itemToId = new Map<string, string>();
  for (let i = 0; i < instruments.length; i++) itemToId.set(items[i], instruments[i]);
  const initialItem = items[instruments.indexOf(last)] ?? items[0];

  const picker = ui.input.choice('Instrument', {
    value: initialItem,
    items,
    nullable: false,
    tooltipText: 'Switches the trend charts and baseline to a different instrument. ' +
      'Only instruments with computed SPC runs are listed.',
  } as any);

  ui.dialog('SPC Dashboard').add(picker).onOK(async () => {
    const picked = itemToId.get(String(picker.value)) ?? instruments[0];
    // userDataStorage is soft-deprecated; matches v1.3 codebase.
    await grok.dapi.userDataStorage.put(LAST_INSTRUMENT_STORAGE_KEY,
      {instrument_id: picked});
    await openDashboardForInstrument(picked);
  }).show();
}

async function openDashboardForInstrument(instrumentId: string): Promise<void> {
  const slice = await loadRuns(instrumentId);
  let runsDf: DG.DataFrame;
  if (slice.length === 0)
    runsDf = DG.DataFrame.create(0);
  else
    runsDf = runsToDataFrame(slice);
  runsDf.name = `SPC — ${instrumentId}`;
  addMrColumns(runsDf);

  const tv = grok.shell.addTableView(runsDf);
  let baseline = await loadBaseline(instrumentId);

  const refreshCharts = (): void => {
    for (const metric of METRIC_KEYS) {
      const {iChart, mrChart} = createSpcChartPanel(runsDf, baseline, metric);
      tv.dockManager.dock(iChart, DG.DOCK_TYPE.RIGHT, null, `I-chart ${METRIC_LABELS[metric]}`, 0.5);
      tv.dockManager.dock(mrChart, DG.DOCK_TYPE.DOWN, null, `MR-chart ${METRIC_LABELS[metric]}`, 0.5);
    }
    // SPC-08 Pareto panel (Plan 16-07). Empty-state body when slice < 4.
    const pareto = createParetoChart(runsDf, instrumentId);
    if (pareto !== null) {
      tv.dockManager.dock(pareto, DG.DOCK_TYPE.DOWN, null, 'Pareto', 0.3);
    } else {
      const empty = ui.divV([
        ui.divText('Pareto chart appears once 4 or more runs have been recorded for this instrument.'),
      ]);
      tv.dockManager.dock(empty, DG.DOCK_TYPE.DOWN, null, 'Pareto', 0.3);
    }
  };

  const onRebuildClick = (): void => {
    showDefineBaselineDialog(instrumentId, slice, async (newBaseline) => {
      await saveBaseline(instrumentId, newBaseline);
      baseline = newBaseline;
      recomputeStatusesInPlace(runsDf, newBaseline);
      refreshCharts();
    }, baseline ?? undefined);
  };

  const onRuleToggle = async (
    metric: MetricKey,
    rule: keyof SpcRulesEnabled['median_intensity'],
    value: boolean,
  ) => {
    if (!baseline) return;
    baseline.rules_enabled[metric][rule] = value;
    await saveBaseline(instrumentId, baseline);
    recomputeStatusesInPlace(runsDf, baseline);
  };

  const sidebar = buildSidebar(instrumentId, baseline, onRebuildClick, onRuleToggle);
  tv.dockManager.dock(sidebar.el, DG.DOCK_TYPE.LEFT, null, 'SPC Controls', 0.25);

  if (baseline === null) {
    const banner = ui.divV([
      ui.h2(`No baseline locked for ${instrumentId}.`),
      ui.divText('Pick the runs that represent normal operation, then click Lock baseline. ' +
        'The dashboard turns on once a baseline exists.'),
      ui.button('Define baseline...', () => {
        showDefineBaselineDialog(instrumentId, slice, async (newBaseline) => {
          await saveBaseline(instrumentId, newBaseline);
          baseline = newBaseline;
          recomputeStatusesInPlace(runsDf, newBaseline);
          refreshCharts();
        });
      }),
    ]);
    tv.dockManager.dock(banner, DG.DOCK_TYPE.RIGHT, null, 'Baseline', 0.75);
  } else {
    refreshCharts();
  }

  wireDrillDown(tv, runsDf);
}
