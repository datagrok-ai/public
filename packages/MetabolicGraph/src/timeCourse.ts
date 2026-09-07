/* eslint-disable camelcase */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import type {BuilderType} from '../escher_src/src/Builder';
import type {CobraModelData, ReactionBounds, SamplingFunctionResult} from '../escher_src/src/ts/types';
import type {PrecomputedExtremePoints} from './cobra/sampler-wrapper';
import {
  sampleFluxAveragesWasm, sampleFluxAveragesPython, computeExtremePointsPython, samplesToDataFrame,
} from './utils';
import type {FluxSamplingResult} from './utils';

/** Sampling parameters reused verbatim from the main "Sample Reactions" dialog. */
export type TimeCourseSamplingParams = {
  samples: number;
  thinning: number;
  bins: number;
  addDf: boolean;
  runInPython: boolean;
  usePythonFBA: boolean;
};

type StepBounds = Map<string, ReactionBounds>;

type TimeCourseStep = {
  index: number;
  reactionData: {[id: string]: number}; // per-reaction average flux -> map colors
  distribution: SamplingFunctionResult; // per-reaction flux histograms -> reaction tooltip
  bounds: {[name: string]: ReactionBounds}; // interpolated bounds applied when this step is selected
  dataSource: string;
};

const delay = (ms: number) => new Promise<void>((r) => setTimeout(r, ms));
const lerp = (a: number, b: number, t: number) => a + (b - a) * t;

/**
 * Push a step's flux distribution into the builder's sampling distribution, mutating the
 * same object the tooltip holds (set once at first load) so hovering shows this step's histogram.
 */
export function applyStepDistribution(builder: BuilderType, dist: SamplingFunctionResult) {
  const target = builder.reaction_sampling_distribution;
  if (!target)
    return;
  target.lower_bound = dist.lower_bound;
  target.upper_bound = dist.upper_bound;
  target.data.clear();
  for (const [k, v] of dist.data)
    target.data.set(k, v);
}

// --------------------------------------------------------------------------
// Managed bottom slider (only one alive at a time; removed on any other action)
// --------------------------------------------------------------------------
let activeSlider: {remove: () => void} | null = null;

/** Remove the time-course slider/animation if one is shown. Safe to call anytime. */
export function clearTimeCourseSlider() {
  if (activeSlider) {
    activeSlider.remove();
    activeSlider = null;
  }
}

/** Whether a time-course slider/animation is currently shown. */
export function isTimeCourseActive(): boolean {
  return activeSlider != null;
}

function showTimeCourseSlider(builder: BuilderType, steps: TimeCourseStep[], originalBounds: {[id: string]: ReactionBounds}) {
  clearTimeCourseSlider();
  const container = (builder.selection?.node() as HTMLElement) ??
    (document.querySelector('.d4-escher-container') as HTMLElement | null);
  if (!container || steps.length === 0)
    return;

  const n = steps.length;
  let current = 0;
  let playing = false;
  let playTimer: any = null;
  let stopOnClick: ((e: Event) => void) | null = null;

  const label = ui.divText(`Step 1/${n}`);
  label.style.minWidth = '74px';
  label.style.fontWeight = 'bold';

  const range = document.createElement('input');
  range.type = 'range';
  range.min = '0';
  range.max = String(n - 1);
  range.step = '1';
  range.value = '0';
  range.style.flex = '1';
  range.style.minWidth = '240px';

  const setStep = (i: number, updateRange = true) => {
    current = Math.max(0, Math.min(n - 1, i));
    if (updateRange)
      range.value = String(current);
    label.textContent = `Step ${current + 1}/${n}`;
    const step = steps[current];
    builder.set_reaction_bounds(step.bounds); // apply this step's interpolated bounds to the model
    applyStepDistribution(builder, step.distribution); // so hovering a reaction shows this step's histogram
    builder.set_reaction_data(step.reactionData, step.dataSource); // update the map colors
  };
  range.addEventListener('input', () => setStep(parseInt(range.value), false));

  const stopPlay = () => {
    playing = false;
    if (playTimer) {
      clearInterval(playTimer);
      playTimer = null;
    }
    if (stopOnClick) {
      document.removeEventListener('mousedown', stopOnClick, true);
      stopOnClick = null;
    }
    playBtn.textContent = '▶ Play';
  };
  const startPlay = () => {
    if (current >= n - 1)
      setStep(0);
    playing = true;
    playBtn.textContent = '⏸ Pause';
    // clicking anywhere outside the panel stops the animation
    stopOnClick = (e: Event) => {
      if (!panel.contains(e.target as Node))
        stopPlay();
    };
    document.addEventListener('mousedown', stopOnClick, true);
    playTimer = setInterval(() => {
      if (current >= n - 1) {
        stopPlay();
        return;
      }
      setStep(current + 1);
    }, 500);
  };
  const playBtn = ui.button('▶ Play', () => (playing ? stopPlay() : startPlay()));

  // quick pass through every step, downloading a PNG of the map for each
  const exportBtn = ui.button('⤓ Export PNGs', async () => {
    stopPlay();
    exportBtn.disabled = playBtn.disabled = range.disabled = true;
    try {
      for (let i = 0; i < n; i++) {
        setStep(i);
        await delay(450); // let the map redraw
        builder.map?.save_png(`time_course_step_${i + 1}_of_${n}`);
        await delay(450); // let the download + svg restore finish
      }
    } finally {
      exportBtn.disabled = playBtn.disabled = range.disabled = false;
    }
  });
  ui.tooltip.bind(exportBtn, 'Quickly step through all steps and download a PNG of the map for each');

  const closeBtn = ui.button('✕', () => clearTimeCourseSlider());

  const panel = ui.divH([label, range, playBtn, exportBtn, closeBtn]);
  Object.assign(panel.style, {
    position: 'absolute', bottom: '12px', left: '50%', transform: 'translateX(-50%)',
    zIndex: '1000', background: 'var(--grey-1, #ffffff)', padding: '8px 12px',
    borderRadius: '8px', boxShadow: '0 2px 10px rgba(0,0,0,0.25)', display: 'flex',
    alignItems: 'center', gap: '10px', width: '640px', maxWidth: '85%',
  });
  container.appendChild(panel);

  activeSlider = {remove: () => {
    stopPlay();
    panel.remove();
    // restore the baseline model bounds and clear the per-step histograms
    if (Object.keys(originalBounds).length)
      builder.set_reaction_bounds(originalBounds);
    builder.reaction_sampling_distribution?.data.clear();
  }};

  setStep(0); // apply the first step (bounds, histograms, colors)
}

// --------------------------------------------------------------------------
// Per-step average flux summary (table + line chart)
// --------------------------------------------------------------------------
/** Reaction ids that carry averages, in model order, with any unknown ids appended. */
function summaryReactionIds(cobraModel: CobraModelData, steps: TimeCourseStep[]): string[] {
  const present = new Set<string>();
  for (const s of steps) {
    for (const id of Object.keys(s.reactionData))
      present.add(id);
  }
  const ordered = cobraModel.reactions.map((r) => r.id).filter((id) => present.has(id));
  const known = new Set(ordered);
  return [...ordered, ...[...present].filter((id) => !known.has(id))];
}

/**
 * Add the time-course summary: one row per step, a numeric `step` column (1..n) followed by one
 * column per reaction holding that step's average flux. A multi-axis line chart is docked over the
 * grid so the trends are what the user sees when the view opens.
 */
function addStepAveragesView(cobraModel: CobraModelData, steps: TimeCourseStep[], highlight: string[]) {
  if (!steps.length)
    return;
  const reactionIds = summaryReactionIds(cobraModel, steps);
  if (!reactionIds.length)
    return;

  const stepCol = DG.Column.int('step', steps.length).init((i) => i + 1);
  const fluxCols = reactionIds.map((id) =>
    DG.Column.float(id, steps.length).init((i) => steps[i].reactionData[id] ?? null));
  const df = DG.DataFrame.fromColumns([stepCol, ...fluxCols]);
  df.name = 'Time-course average flux';

  const title = `Average flux per step (${steps.length} steps)`;
  // default the chart to the reactions whose bounds were interpolated — those are the ones that move
  const inDf = new Set(reactionIds);
  const preferred = highlight.filter((id) => inDf.has(id));
  const yColumnNames = (preferred.length ? preferred : reactionIds).slice(0, 2);

  const view = grok.shell.addTableView(df);
  const chart = DG.Viewer.lineChart(df, {
    title,
    xColumnName: 'step',
    yColumnNames,
    multiAxis: true,
    yAxisTitle: 'Average flux',
    multiAxisLegendPosition: 'RightCenter',
  });
  // ratio 1 => the chart takes the whole view, hiding the grid underneath
  view.dockManager.dock(chart, DG.DOCK_TYPE.FILL, null, title, 1);
}

// --------------------------------------------------------------------------
// Bounds CSV parsing & validation
// --------------------------------------------------------------------------
function parseBoundsCsv(text: string): {bounds?: StepBounds; error?: string} {
  let df: DG.DataFrame;
  try {
    df = DG.DataFrame.fromCsv(text);
  } catch (_) {
    return {error: 'could not parse CSV'};
  }
  if (df.columns.length !== 3)
    return {error: `expected 3 columns (reaction, lower_bound, upper_bound), found ${df.columns.length}`};
  const nameCol = df.columns.byIndex(0);
  const lbCol = df.columns.byIndex(1);
  const ubCol = df.columns.byIndex(2);
  if (!nameCol.isCategorical || !lbCol.isNumerical || !ubCol.isNumerical)
    return {error: 'columns must be reaction (text), lower_bound (number), upper_bound (number)'};
  const bounds: StepBounds = new Map();
  for (let i = 0; i < df.rowCount; i++) {
    const name = nameCol.get(i);
    const lb = lbCol.get(i);
    const ub = ubCol.get(i);
    if (typeof name !== 'string' || typeof lb !== 'number' || typeof ub !== 'number')
      continue;
    bounds.set(name, {lower_bound: lb, upper_bound: ub});
  }
  if (bounds.size === 0)
    return {error: 'no valid reaction rows found'};
  return {bounds};
}

function fileRow(label: string, onFile: (text: string | null) => void): HTMLElement {
  const input = document.createElement('input');
  input.type = 'file';
  input.accept = '.csv,text/csv';
  input.addEventListener('change', async () => {
    const f = input.files?.[0];
    if (!f) {
      onFile(null);
      return;
    }
    try {
      onFile(await f.text());
    } catch (_) {
      onFile(null);
    }
  });
  const lbl = ui.divText(label);
  lbl.style.minWidth = '200px';
  const row = ui.divH([lbl, input]);
  Object.assign(row.style, {alignItems: 'center', gap: '8px', margin: '6px 0'});
  return row;
}

// --------------------------------------------------------------------------
// Manual bounds editor (two side-by-side grids kept in reaction parity)
// --------------------------------------------------------------------------
/** Amount added/removed by the up-/down-regulate icons. */
const REGULATE_STEP = 10;

/** Per-action colors for the row action icons. */
const ICON_COLORS = {
  add: '#2083d5', // blue
  delete: '#d4504e', // red
  knockout: '#7d7d7d', // grey — zeroes the reaction out
  up: '#2ea84f', // green
  down: '#e08a2e', // orange
};

/** One reaction shared by both days. Day 1 and Day N bounds are edited independently. */
type ManualEntry = {reaction: string | null; d1: ReactionBounds; dn: ReactionBounds};

const newManualEntry = (): ManualEntry =>
  ({reaction: null, d1: {lower_bound: 0, upper_bound: 0}, dn: {lower_bound: 0, upper_bound: 0}});

/** Strip the (empty) label so the input sits flush inside a labelled grid. */
function gridCell(root: HTMLElement): HTMLElement {
  root.querySelector('label')?.remove();
  root.style.width = '100%';
  return root;
}

// --------------------------------------------------------------------------
// Dialog
// --------------------------------------------------------------------------
export function openTimeCourseDialog(cobraModel: CobraModelData | null, builder: BuilderType, params: TimeCourseSamplingParams) {
  if (!cobraModel || !cobraModel.reactions?.length) {
    grok.shell.error('No model loaded — cannot run time-course sampling');
    return;
  }

  // map every identifier a bounds file may reference (id/name/bigg_id) to the canonical reaction id,
  // and remember each reaction's model bounds (used as defaults when a reaction is picked)
  const idLookup = new Map<string, string>();
  const modelBoundsById = new Map<string, ReactionBounds>();
  for (const r of cobraModel.reactions) {
    idLookup.set(r.id, r.id);
    if (r.name) idLookup.set(r.name, r.id);
    if (r.bigg_id) idLookup.set(r.bigg_id, r.id);
    modelBoundsById.set(r.id, {lower_bound: r.lower_bound ?? 0, upper_bound: r.upper_bound ?? 0});
  }
  const modelIds = new Set(idLookup.keys());
  const reactionIds = cobraModel.reactions.map((r) => r.id);
  const modelBounds = (id: string | null): ReactionBounds =>
    (id && modelBoundsById.get(id)) ? {...modelBoundsById.get(id)!} : {lower_bound: 0, upper_bound: 0};

  let err1 = '';
  let err2 = '';
  // manually entered reactions, shared by both days to keep parity
  let entries: ManualEntry[] = [newManualEntry()];

  const dialog = ui.dialog('Time-course Sampling (interpolate bounds)');
  const stepsInput = ui.input.int('Steps', {
    value: 10,
    tooltipText: 'Number of interpolation steps, including both endpoints (min 2)',
  });
  const errorDiv = ui.divText('');
  Object.assign(errorDiv.style, {
    color: 'var(--red-3, #b00020)', whiteSpace: 'pre-wrap', marginTop: '8px',
    maxWidth: '560px', display: 'none',
  });

  // --- manual grids ---------------------------------------------------------
  const grid1 = ui.div([], {classes: 'mg-bounds-grid'});
  const gridN = ui.div([], {classes: 'mg-bounds-grid'});
  for (const g of [grid1, gridN]) {
    Object.assign(g.style, {
      display: 'grid', gridTemplateColumns: 'minmax(110px, 1fr) 64px 64px auto',
      alignItems: 'center', columnGap: '8px', rowGap: '6px',
    });
  }

  /** Collect a day's manual rows into a bounds map; reports duplicate reactions. */
  const entriesToBounds = (pick: (e: ManualEntry) => ReactionBounds): {bounds: StepBounds; dupes: string[]} => {
    const bounds: StepBounds = new Map();
    const dupes: string[] = [];
    for (const e of entries) {
      if (!e.reaction) continue;
      if (bounds.has(e.reaction)) dupes.push(e.reaction);
      const b = pick(e);
      bounds.set(e.reaction, {lower_bound: b.lower_bound, upper_bound: b.upper_bound});
    }
    return {bounds, dupes};
  };

  /** Resolve the effective bounds for each day from the (single source of truth) grid. */
  const resolveBounds = () => {
    const m1 = entriesToBounds((e) => e.d1);
    const mN = entriesToBounds((e) => e.dn);
    const start = m1.bounds.size ? m1.bounds : null;
    const end = mN.bounds.size ? mN.bounds : null;
    return {start, end, dupes: [...new Set([...m1.dupes, ...mN.dupes])]};
  };

  const validate = () => {
    const errors: string[] = [];
    if (err1) errors.push(`Start file: ${err1}`);
    if (err2) errors.push(`End file: ${err2}`);

    const {start, end, dupes} = resolveBounds();
    if (dupes.length) errors.push(`Reaction listed more than once: ${dupes.join(', ')}`);
    if (!start && !err1) errors.push('Provide start (Day 1) bounds — upload a CSV or add reactions below');
    if (!end && !err2) errors.push('Provide end (Day N) bounds — upload a CSV or add reactions below');

    if (start && end) {
      const missing1 = [...start.keys()].filter((k) => !modelIds.has(k));
      const missing2 = [...end.keys()].filter((k) => !modelIds.has(k));
      if (missing1.length) errors.push(`Start: reactions not in model: ${missing1.join(', ')}`);
      if (missing2.length) errors.push(`End: reactions not in model: ${missing2.join(', ')}`);
      const only1 = [...start.keys()].filter((k) => !end.has(k));
      const only2 = [...end.keys()].filter((k) => !start.has(k));
      if (only1.length || only2.length)
        errors.push(`Start and end must list the same reactions. Mismatch: ${[...only1, ...only2].join(', ')}`);
    }
    const steps = stepsInput.value ?? 0;
    if (!steps || steps < 2) errors.push('Steps must be at least 2');

    errorDiv.textContent = errors.join('\n');
    errorDiv.style.display = errors.length ? 'block' : 'none';
    const okBtn = dialog.getButton('OK');
    if (okBtn) okBtn.disabled = errors.length > 0;
  };

  /** Build one day's grid: header, a row per shared entry, and a trailing add row. */
  const renderGrid = (grid: HTMLElement, day: 'd1' | 'dn') => {
    ui.empty(grid);
    for (const caption of ['Reaction', 'Lower', 'Upper', '']) {
      const h = ui.divText(caption);
      h.style.fontWeight = 'bold';
      grid.appendChild(h);
    }

    for (const entry of entries) {
      const recInput = ui.input.choice('', {items: reactionIds, value: entry.reaction, nullable: true});
      // reaction is shared across both days; picking one seeds both days with the model's bounds,
      // then re-renders both grids to keep parity
      recInput.onChanged.subscribe(() => {
        entry.reaction = recInput.value ?? null;
        entry.d1 = modelBounds(entry.reaction);
        entry.dn = modelBounds(entry.reaction);
        renderBoth();
      });
      const loInput = ui.input.float('', {value: entry[day].lower_bound});
      loInput.onChanged.subscribe(() => {
        entry[day].lower_bound = loInput.value ?? 0;
        validate();
      });
      const upInput = ui.input.float('', {value: entry[day].upper_bound});
      upInput.onChanged.subscribe(() => {
        entry[day].upper_bound = upInput.value ?? 0;
        validate();
      });

      const delIcon = ui.icons.delete(() => {
        entries.splice(entries.indexOf(entry), 1);
        renderBoth();
      }, 'Remove reaction (from both days)');
      delIcon.style.color = ICON_COLORS.delete;
      const koIcon = ui.iconFA('ban', () => {
        entry[day] = {lower_bound: 0, upper_bound: 0};
        renderBoth();
      }, 'Knock out (set lower and upper to 0)');
      koIcon.style.color = ICON_COLORS.knockout;
      const upIcon = ui.iconFA('arrow-up', () => {
        entry[day].lower_bound += REGULATE_STEP;
        entry[day].upper_bound += REGULATE_STEP;
        renderBoth();
      }, `Up-regulate (+${REGULATE_STEP} to lower and upper)`);
      upIcon.style.color = ICON_COLORS.up;
      const downIcon = ui.iconFA('arrow-down', () => {
        entry[day].lower_bound -= REGULATE_STEP;
        entry[day].upper_bound -= REGULATE_STEP;
        renderBoth();
      }, `Down-regulate (−${REGULATE_STEP} from lower and upper)`);
      downIcon.style.color = ICON_COLORS.down;
      const icons = ui.divH([delIcon, koIcon, upIcon, downIcon]);
      icons.style.gap = '8px';

      grid.appendChild(gridCell(recInput.root));
      grid.appendChild(gridCell(loInput.root));
      grid.appendChild(gridCell(upInput.root));
      grid.appendChild(icons);
    }

    // trailing add row: a new reaction appears in both days at once (parity)
    const addIcon = ui.icons.add(() => {
      entries.push(newManualEntry());
      renderBoth();
    }, 'Add reaction (to both days)');
    addIcon.style.color = ICON_COLORS.add;
    const addRow = ui.divH([addIcon]);
    addRow.style.gridColumn = '1 / -1';
    grid.appendChild(addRow);
  };

  const renderBoth = () => {
    renderGrid(grid1, 'd1');
    renderGrid(gridN, 'dn');
    validate();
  };

  const section = (title: string, fileEl: HTMLElement, grid: HTMLElement) => {
    const header = ui.divText(title);
    Object.assign(header.style, {fontWeight: 'bold', marginBottom: '4px'});
    const col = ui.divV([header, fileEl, grid]);
    Object.assign(col.style, {flex: '1', minWidth: '320px', gap: '6px'});
    return col;
  };

  /** Merge an uploaded day's bounds into the shared grid, upserting rows by canonical reaction id. */
  const mergeFileBounds = (bounds: StepBounds, day: 'd1' | 'dn') => {
    for (const [rawId, b] of bounds) {
      const id = idLookup.get(rawId) ?? rawId; // canonicalize name/bigg_id -> reaction id
      let entry = entries.find((e) => e.reaction === id);
      if (!entry) {
        // new reaction: seed the other day from the model so both days stay in parity
        entry = {reaction: id, d1: modelBounds(id), dn: modelBounds(id)};
        entries.push(entry);
      }
      entry[day] = {lower_bound: b.lower_bound, upper_bound: b.upper_bound};
    }
    // drop the leftover blank starter row once real rows exist
    if (entries.length > 1)
      entries = entries.filter((e) => e.reaction !== null);
    renderBoth();
  };

  const file1Row = fileRow('Start bounds (Day 1) CSV', (text) => {
    err1 = '';
    if (text != null) {
      const r = parseBoundsCsv(text);
      if (r.bounds) mergeFileBounds(r.bounds, 'd1');
      err1 = r.error ?? '';
    }
    validate();
  });
  const file2Row = fileRow('End bounds (Day N) CSV', (text) => {
    err2 = '';
    if (text != null) {
      const r = parseBoundsCsv(text);
      if (r.bounds) mergeFileBounds(r.bounds, 'dn');
      err2 = r.error ?? '';
    }
    validate();
  });
  stepsInput.onChanged.subscribe(() => validate());

  const hint = ui.divText('Upload CSVs and/or edit reaction bounds below. Uploading a CSV fills the grid for that day; you can then edit individual rows.');
  Object.assign(hint.style, {color: 'var(--grey-5)', maxWidth: '700px', marginBottom: '8px'});

  const divider = ui.div([]);
  Object.assign(divider.style, {width: '1px', alignSelf: 'stretch', background: 'var(--grey-3, #d0d0d0)'});

  const grids = ui.divH([
    section('Day 1 (start)', file1Row, grid1),
    divider,
    section('Day N (end)', file2Row, gridN),
  ]);
  grids.style.gap = '24px';
  grids.style.alignItems = 'stretch';

  dialog.add(ui.divV([hint, grids, stepsInput.root, errorDiv]));
  renderBoth();

  dialog.onOK(async () => {
    const {start, end} = resolveBounds();
    if (!start || !end || (stepsInput.value ?? 0) < 2)
      return;
    await runTimeCourseSampling(cobraModel, builder, start, end, stepsInput.value!, params);
  });

  // persist/restore the manually entered bounds and step count via dialog history.
  // history only round-trips flat primitive values, so the rows are serialized to a JSON string.
  dialog.history(
    () => ({
      steps: stepsInput.value ?? 10,
      entries: JSON.stringify(entries.map((e) => ({
        reaction: e.reaction, d1Lower: e.d1.lower_bound, d1Upper: e.d1.upper_bound,
        dnLower: e.dn.lower_bound, dnUpper: e.dn.upper_bound,
      }))),
    }),
    (x) => {
      if (!x) return;
      if (typeof x.steps === 'number') stepsInput.value = x.steps;
      let parsed: any[] | null = null;
      try {
        if (typeof x.entries === 'string') parsed = JSON.parse(x.entries);
        else if (Array.isArray(x.entries)) parsed = x.entries;
      } catch (_) { parsed = null; }
      if (parsed) {
        entries = parsed.length ? parsed.map((e: any) => ({
          reaction: e.reaction ?? null,
          d1: {lower_bound: e.d1Lower ?? 0, upper_bound: e.d1Upper ?? 0},
          dn: {lower_bound: e.dnLower ?? 0, upper_bound: e.dnUpper ?? 0},
        })) : [newManualEntry()];
      }
      renderBoth();
    },
  );

  dialog.show({center: true});
  dialog.root.style.minWidth = '720px';
  const okBtn = dialog.getButton('OK');
  if (okBtn) okBtn.disabled = true;
  // keep keystrokes from reaching escher's key manager
  dialog.root.addEventListener('keydown', (e) => e.stopPropagation());
  validate();
}

// --------------------------------------------------------------------------
// Runner: interpolate bounds, sample every step, then show the slider
// --------------------------------------------------------------------------
export async function runTimeCourseSampling(
  cobraModel: CobraModelData, builder: BuilderType,
  bounds1: StepBounds, bounds2: StepBounds, steps: number, params: TimeCourseSamplingParams,
): Promise<{succeeded: boolean, failedSteps: number[]}> {
  clearTimeCourseSlider();
  const reactionNames = [...bounds1.keys()];

  // snapshot original bounds of the affected reactions so we can restore them afterwards
  const affected = cobraModel.reactions.filter((r) =>
    reactionNames.some((nm) => nm === r.id || nm === r.name || nm === r.bigg_id));
  const original: {[id: string]: ReactionBounds} = {};
  for (const r of affected)
    original[r.id] = {lower_bound: r.lower_bound ?? 0, upper_bound: r.upper_bound ?? 0};

  const pg = DG.TaskBarProgressIndicator.create('Time-course sampling');
  const stepBoundsArr: {[name: string]: ReactionBounds}[] = new Array(steps);
  type StepData = {reactionData: {[id: string]: number}; distribution: SamplingFunctionResult};
  const sampled: (StepData | null)[] = new Array(steps).fill(null);
  const failed: number[] = [];
  try {
    for (let i = 0; i < steps; i++) {
      const t = steps === 1 ? 0 : i / (steps - 1);
      const stepBounds: {[name: string]: ReactionBounds} = {};
      for (const nm of reactionNames) {
        const b1 = bounds1.get(nm)!;
        const b2 = bounds2.get(nm)!;
        stepBounds[nm] = {
          lower_bound: lerp(b1.lower_bound, b2.lower_bound, t),
          upper_bound: lerp(b1.upper_bound, b2.upper_bound, t),
        };
      }
      stepBoundsArr[i] = stepBounds; // keep the interpolated bounds even if this step fails to sample
      pg.update(Math.round((100 * i) / steps), `Sampling step ${i + 1}/${steps}`);
      // isolate each step: an infeasible/unstable step must not abort the whole run
      try {
        builder.set_reaction_bounds(stepBounds);
        const model = builder.model_data!;
        let res: FluxSamplingResult | null;
        if (params.runInPython)
          res = await sampleFluxAveragesPython(model, params.bins, params.samples, params.thinning);
        else {
          let precomputed: PrecomputedExtremePoints | undefined;
          if (params.usePythonFBA)
            precomputed = await computeExtremePointsPython(model);
          res = await sampleFluxAveragesWasm(model, params.bins, params.samples, params.thinning, precomputed);
        }
        if (!res)
          throw new Error('infeasible reaction bounds');
        sampled[i] = {reactionData: res.reactionData, distribution: res.distribution};
        if (params.addDf) {
          const df = res.df ?? samplesToDataFrame(model, res.results!, params.samples);
          df.name = `Samples step ${i + 1} of ${steps}`;
          grok.shell.addTableView(df);
        }
      } catch (e) {
        console.error(`Time-course step ${i + 1} failed:`, e);
        failed.push(i);
      }
    }
  } finally {
    pg.close();
  }

  if (sampled.every((d) => d === null)) {
    // nothing feasible: restore the baseline model and bail
    if (Object.keys(original).length)
      builder.set_reaction_bounds(original);
    grok.shell.error('Time-course sampling failed for every step — check that the interpolated bounds are feasible');
    return {succeeded: false, failedSteps: failed};
  }

  // Keep exactly `steps` slider positions: fill a failed step with the nearest sampled one.
  const nearest = (i: number): StepData => {
    for (let d = 0; d < steps; d++) {
      if (i - d >= 0 && sampled[i - d]) return sampled[i - d]!;
      if (i + d < steps && sampled[i + d]) return sampled[i + d]!;
    }
    return {reactionData: {}, distribution: {upper_bound: 0, lower_bound: 0, data: new Map<string, number[]>()}};
  };
  const stepResults: TimeCourseStep[] = [];
  for (let i = 0; i < steps; i++) {
    const ok = sampled[i] !== null;
    const data = ok ? sampled[i]! : nearest(i);
    stepResults.push({
      index: i,
      reactionData: data.reactionData,
      distribution: data.distribution,
      bounds: stepBoundsArr[i],
      dataSource: ok ?
        `Time-course step ${i + 1}/${steps}` :
        `Time-course step ${i + 1}/${steps} (infeasible — showing nearest sampled step)`,
    });
  }
  if (failed.length)
    grok.shell.warning(`Steps ${failed.map((i) => i + 1).join(', ')} had infeasible bounds and could not be sampled; their frames reuse the nearest sampled step.`);

  addStepAveragesView(cobraModel, stepResults, affected.map((r) => r.id));

  // expose the managed scrubber; it selects step 1 (applying its bounds, histograms and colors)
  showTimeCourseSlider(builder, stepResults, original);
  return {succeeded: true, failedSteps: failed};
}
