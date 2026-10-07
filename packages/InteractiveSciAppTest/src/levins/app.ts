// Levins Metapopulation Model — Application (Coordinator + UI)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULTS, validate, solve, validateOptimize,
  LevinsParams, LevinsSolution, InputId, WorkerTask, WorkerResult,
} from './core';

import '../../css/levins.css';

const DEBOUNCE_MS = 50;
const OPTIMIZE_POINTS = 10000;

export function levinsMetapopulationApp(): void {
  // --- State ---
  let computationsBlocked = false;
  let debounceTimer: ReturnType<typeof setTimeout> | null = null;
  const subs: {unsubscribe(): void}[] = [];
  let activeWorkers: Worker[] = [];
  let lineChart!: DG.Viewer;

  // --- Initial DataFrame ---
  const initSolution = solve(DEFAULTS);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromFloat64Array('t', initSolution.t),
    DG.Column.fromFloat64Array('p', initSolution.p),
  ]);
  df.name = 'Levins Metapopulation';

  const view = grok.shell.addTableView(df);
  view.name = 'Levins Metapopulation Model';

  // --- Rho badge (created before controls so onValueChanged callbacks can reference it) ---
  const rhoBadge = ui.div([], 'd4-tag levins-rho-badge');
  ui.tooltip.bind(rhoBadge, 'Extinction-to-colonization rate ratio. \u03C1 < 1 \u2014 metapopulation persists, \u03C1 \u2265 1 \u2014 extinction.');

  // --- Controls ---

  // Initial condition
  const ctrlP0 = ui.input.float('Initial patch fraction p₀', {
    value: DEFAULTS.p0, nullable: false,
    min: 0.001, max: 1,
    tooltipText: 'Fraction of patches occupied at t=0. If p₀=0, the population cannot recover — computation is skipped.',
    onValueChanged: () => debouncedRun(),
  });

  // Parameters
  const ctrlM = ui.input.float('Colonization rate m', {
    value: DEFAULTS.m, nullable: false,
    min: 0.001, max: 100,
    tooltipText: 'How fast empty patches are colonized from occupied ones. Units: 1/time.',
    onValueChanged: () => { updateRhoBadge(); debouncedRun(); },
  });

  const ctrlE0 = ui.input.float('Extinction rate e₀', {
    value: DEFAULTS.e0, nullable: false,
    min: 0.001, max: 100,
    tooltipText: 'Base local extinction rate of a subpopulation in a patch. With rescue effect — decreases as p grows. Units: 1/time.',
    onValueChanged: () => { updateRhoBadge(); debouncedRun(); },
  });

  const ctrlRescue = ui.input.toggle('Rescue effect', {
    value: DEFAULTS.rescueEffect,
    tooltipText: 'When enabled, extinction rate depends on p: e(p) = e₀·(1−p). More occupied patches — lower local extinction.',
    onValueChanged: () => { updateRescueLabel(); runPrimary(); },
  });

  // Argument
  const ctrlTStart = ui.input.float('Start t₀', {
    value: DEFAULTS.t_start, nullable: false,
    min: 0, max: 10000,
    tooltipText: 'Simulation start time. Usually 0.',
    onValueChanged: () => { updateArgRanges(); debouncedRun(); },
  });

  const ctrlTEnd = ui.input.float('End t_end', {
    value: DEFAULTS.t_end, nullable: false,
    min: 0.1, max: 10000,
    tooltipText: 'Simulation end time. Recommended ≥ 5/e₀ so the system reaches equilibrium.',
    onValueChanged: () => { updateArgRanges(); debouncedRun(); },
  });

  const ctrlTStep = ui.input.float('Step Δt', {
    value: DEFAULTS.t_step, nullable: false,
    min: 0.001, max: 1000,
    tooltipText: 'Grid step of the numerical solution. Affects chart detail, not stability (MRT is an implicit method).',
    onValueChanged: () => debouncedRun(),
  });

  // Solver
  const ctrlTolerance = ui.input.float('Tolerance', {
    value: DEFAULTS.tolerance, nullable: false,
    min: 1e-12, max: 1e-2,
    tooltipText: 'MRT method numerical tolerance. Lower — more precise but slower. Recommended: 1e-6 … 1e-9.',
    onValueChanged: () => debouncedRun(),
  });

  // Set formats (block computations to avoid spurious runs from format-triggered events)
  computationsBlocked = true;
  ctrlP0.format = '0.000';
  ctrlM.format = '0.000';
  ctrlE0.format = '0.000';
  ctrlTStart.format = '0.0';
  ctrlTEnd.format = '0.0';
  ctrlTStep.format = '0.000';
  ctrlTolerance.format = '0.##E+0';
  computationsBlocked = false;

  // --- Input map for validators ---
  const inputMap: Record<InputId, DG.InputBase> = {
    'ctrl_p0': ctrlP0,
    'ctrl_m': ctrlM,
    'ctrl_e0': ctrlE0,
    'ctrl_rescue': ctrlRescue,
    'ctrl_t_start': ctrlTStart,
    'ctrl_t_end': ctrlTEnd,
    'ctrl_t_step': ctrlTStep,
    'ctrl_tolerance': ctrlTolerance,
  };

  function updateRhoBadge(): void {
    const m = ctrlM.value ?? DEFAULTS.m;
    const e0 = ctrlE0.value ?? DEFAULTS.e0;
    const rho = e0 / m;
    const persists = rho < 1;
    rhoBadge.textContent = `ρ = e₀/m = ${rho.toFixed(3)}`;
    rhoBadge.classList.toggle('levins-rho-badge--persists', persists);
    rhoBadge.classList.toggle('levins-rho-badge--extinct', !persists);
  }
  updateRhoBadge();

  // --- Rescue effect label reactivity ---
  function updateRescueLabel(): void {
    if (ctrlRescue.value) {
      ctrlE0.caption = 'Base extinction rate e₀';
      ctrlE0.setTooltip('Base local extinction rate. Effective rate: e(p) = e₀·(1−p)');
    } else {
      ctrlE0.caption = 'Extinction rate e₀';
      ctrlE0.setTooltip('Base local extinction rate of a subpopulation in a patch. Units: 1/time.');
    }
  }

  // --- Argument range reactivity ---
  // Note: Datagrok InputBase does not have setOptions for changing min/max at runtime.
  // Range validation is handled by the complex validator instead.
  function updateArgRanges(): void {
    // Ranges are enforced through validation (val_06, val_07)
  }

  // --- Gather current inputs ---
  function getInputs(): LevinsParams {
    return {
      p0: ctrlP0.value ?? DEFAULTS.p0,
      m: ctrlM.value ?? DEFAULTS.m,
      e0: ctrlE0.value ?? DEFAULTS.e0,
      rescueEffect: ctrlRescue.value ?? DEFAULTS.rescueEffect,
      t_start: ctrlTStart.value ?? DEFAULTS.t_start,
      t_end: ctrlTEnd.value ?? DEFAULTS.t_end,
      t_step: ctrlTStep.value ?? DEFAULTS.t_step,
      tolerance: ctrlTolerance.value ?? DEFAULTS.tolerance,
    };
  }

  // --- Validators ---
  function addValidators(): void {
    const validatorFor = (id: InputId) => {
      return () => {
        const inputs = getInputs();
        const errors = validate(inputs);
        return errors.get(id) ?? null;
      };
    };

    ctrlP0.addValidator(validatorFor('ctrl_p0'));
    ctrlM.addValidator(validatorFor('ctrl_m'));
    ctrlE0.addValidator(validatorFor('ctrl_e0'));
    ctrlTStart.addValidator(validatorFor('ctrl_t_start'));
    ctrlTEnd.addValidator(validatorFor('ctrl_t_end'));
    ctrlTStep.addValidator(validatorFor('ctrl_t_step'));
    ctrlTolerance.addValidator(validatorFor('ctrl_tolerance'));
  }
  addValidators();

  // --- Color coding ---
  function updateColorCoding(): void {
    const m = ctrlM.value ?? DEFAULTS.m;
    const e0 = ctrlE0.value ?? DEFAULTS.e0;
    const threshold = e0 / m;
    const pCol = view.dataFrame.col('p');
    if (pCol == null) return;

    const rules: Record<string, string> = {};
    rules['<' + threshold] = '#F44336';
    rules['>=' + threshold] = '#4CAF50';
    pCol.meta.colors.setConditional(rules);
  }

  // --- Grid column header tooltip (via onCellTooltip, as in EDA) ---
  function setupGridTooltip(): void {
    view.grid.onCellTooltip((cell, x, y) => {
      if (!cell.isColHeader)
        return false;

      const colName = cell.tableColumn?.name;
      if (colName === 'p') {
        const m = ctrlM.value ?? DEFAULTS.m;
        const e0 = ctrlE0.value ?? DEFAULTS.e0;
        const rho = (e0 / m).toFixed(3);
        ui.tooltip.show(ui.divV([
          ui.h2('Occupied patch fraction p(t)'),
          ui.divText(`Color: green — persistence zone (p ≥ e₀/m)`),
          ui.divText(`red — extinction threat zone (p < e₀/m)`),
          ui.divText(`Threshold: e₀/m = ${rho}`),
        ]), x, y);
        return true;
      }

      return false;
    });
  }

  // --- Primary pipeline ---
  function runPrimary(): void {
    if (computationsBlocked)
      return;

    const inputs = getInputs();
    const errors = validate(inputs);

    // Clear previous errors on all inputs
    for (const input of Object.values(inputMap))
      input.input?.classList.remove('d4-invalid');

    if (errors.size > 0) {
      errors.forEach((_msg, id) => {
        const input = inputMap[id];
        if (input)
          input.input?.classList.add('d4-invalid');
      });
      clearResults();
      return;
    }

    try {
      const result = solve(inputs);
      updateDataFrame(result);
      updateColorCoding();
    } catch (err) {
      clearResults();
      const msg = err instanceof Error ? err.message : 'Computation error';
      grok.shell.error(msg);
    }
  }

  function debouncedRun(): void {
    if (debounceTimer !== null)
      clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => runPrimary(), DEBOUNCE_MS);
  }

  // --- Update DataFrame ---
  function updateDataFrame(result: LevinsSolution): void {
    const newDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('t', result.t),
      DG.Column.fromFloat64Array('p', result.p),
    ]);
    newDf.name = 'Levins Metapopulation';
    view.dataFrame = newDf;
    lineChart.dataFrame = newDf;
  }

  function clearResults(): void {
    const emptyDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('t', new Float64Array(0)),
      DG.Column.fromFloat64Array('p', new Float64Array(0)),
    ]);
    emptyDf.name = 'Levins Metapopulation';
    view.dataFrame = emptyDf;
    lineChart.dataFrame = emptyDf;
  }

  // --- Optimization task ---
  let optimizeBtn: HTMLElement;

  function setOptimizeBtnEnabled(enabled: boolean): void {
    optimizeBtn.classList.toggle('levins-btn--disabled', !enabled);
  }

  async function runOptimization(mMin: number, mMax: number): Promise<void> {
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Internal error: invalid task parameters. Check the inputs and try again.');
      return;
    }

    setOptimizeBtnEnabled(false);

    const pi = DG.TaskBarProgressIndicator.create('Optimizing m...');
    let canceled = false;
    let completed = 0;
    let errorCount = 0;

    const results: {m_i: number; p_end: number}[] = [];
    const workerCount = Math.max(1, (navigator.hardwareConcurrency ?? 4) - 2);

    // Generate m values
    const mValues: number[] = [];
    for (let i = 0; i < OPTIMIZE_POINTS; i++)
      mValues.push(mMin + i * (mMax - mMin) / (OPTIMIZE_POINTS - 1));

    // Worker pool
    const workerUrl = _package.webRoot + 'dist/optimize-worker.js';

    try {
      await new Promise<void>((resolve, reject) => {
        const taskQueue = [...mValues];
        activeWorkers = [];

        const createWorker = (): Worker | null => {
          try {
            const worker = new Worker(workerUrl);
            activeWorkers.push(worker);
            return worker;
          } catch (_err) {
            return null;
          }
        };

        const processNext = (worker: Worker) => {
          if (canceled) {
            terminateWorkers();
            resolve();
            return;
          }

          if (taskQueue.length === 0) {
            worker.terminate();
            activeWorkers = activeWorkers.filter((w) => w !== worker);
            if (activeWorkers.length === 0)
              resolve();
            return;
          }

          const m_i = taskQueue.shift()!;
          const task: WorkerTask = {
            m_i,
            p0: inputs.p0,
            e0: inputs.e0,
            rescueEffect: inputs.rescueEffect,
            t_start: inputs.t_start,
            t_end: inputs.t_end,
            t_step: inputs.t_step,
            tolerance: inputs.tolerance,
          };
          worker.postMessage(task);
        };

        const handleResult = (worker: Worker, event: MessageEvent<WorkerResult>) => {
          const result = event.data;
          completed++;
          pi.update(Math.round(completed / OPTIMIZE_POINTS * 100), `${completed}/${OPTIMIZE_POINTS}`);

          if (result.error)
            errorCount++;
          else
            results.push({m_i: result.m_i, p_end: result.p_end});

          processNext(worker);
        };

        // Create worker pool
        for (let i = 0; i < workerCount; i++) {
          const worker = createWorker();
          if (worker == null) {
            if (i === 0) {
              reject(new Error('Failed to start parallel computations. Try again later.'));
              return;
            }
            break;
          }

          worker.onmessage = (event) => handleResult(worker, event);
          worker.onerror = () => {
            completed++;
            errorCount++;
            pi.update(Math.round(completed / OPTIMIZE_POINTS * 100), `${completed}/${OPTIMIZE_POINTS}`);
            processNext(worker);
          };

          processNext(worker);
        }
      });
    } catch (err) {
      grok.shell.error(err instanceof Error ? err.message : 'Failed to start parallel computations.');
      pi.close();
      setOptimizeBtnEnabled(true);
      return;
    }

    pi.close();
    terminateWorkers();
    setOptimizeBtnEnabled(true);

    if (canceled)
      return;

    // Handle results
    if (errorCount > 0 && errorCount < OPTIMIZE_POINTS)
      grok.shell.warning(`${errorCount} of ${OPTIMIZE_POINTS} points failed to compute. Result based on ${OPTIMIZE_POINTS - errorCount} points.`);

    if (results.length === 0) {
      grok.shell.error('Failed to compute any point. Check the parameters.');
      return;
    }

    // Find optimal
    let best = results[0];
    for (const r of results) {
      if (r.p_end > best.p_end)
        best = r;
    }

    // Batch update: block primary, write m_optimal, unblock and run once
    try {
      computationsBlocked = true;
      ctrlM.value = best.m_i;
      computationsBlocked = false;
      runPrimary();
      grok.shell.info(`Optimal m = ${best.m_i.toFixed(3)}\np(t_end) = ${best.p_end.toFixed(3)}`);
    } catch (_err) {
      computationsBlocked = false;
      grok.shell.warning(`Optimal m = ${best.m_i.toFixed(3)}, but failed to update the field automatically. Enter the value manually.`);
    }
  }

  function terminateWorkers(): void {
    for (const w of activeWorkers)
      w.terminate();
    activeWorkers = [];
  }

  // --- Optimize dialog ---
  function showOptimizeDialog(): void {
    const dlgMMin = ui.input.float('Minimum m', {
      value: 0.1, nullable: false,
      min: 0.001, max: 100,
      tooltipText: 'Lower bound of the m search range. Must be less than the maximum value.',
    });
    dlgMMin.format = '0.000';

    const dlgMMax = ui.input.float('Maximum m', {
      value: 1.0, nullable: false,
      min: 0.001, max: 100,
      tooltipText: 'Upper bound of the m search range. Must be greater than the minimum value.',
    });
    dlgMMax.format = '0.000';

    // Cross-validation of dialog inputs
    const validateDialog = (): boolean => {
      const mMin = dlgMMin.value ?? 0.1;
      const mMax = dlgMMax.value ?? 1.0;
      const e0 = ctrlE0.value ?? DEFAULTS.e0;
      const rescueEffect = ctrlRescue.value ?? DEFAULTS.rescueEffect;

      const {errors, warning} = validateOptimize({m_min: mMin, m_max: mMax}, e0, rescueEffect);

      if (warning)
        grok.shell.warning(warning);

      return errors.size === 0;
    };

    dlgMMin.addValidator(() => {
      const mMin = dlgMMin.value ?? 0.1;
      const mMax = dlgMMax.value ?? 1.0;
      if (mMin <= 0) return 'Colonization rate must be positive';
      if (mMin >= mMax) return 'Minimum value must be less than maximum';
      return null;
    });

    dlgMMax.addValidator(() => {
      const mMax = dlgMMax.value ?? 1.0;
      if (mMax <= 0) return 'Colonization rate must be positive';
      return null;
    });

    ui.dialog('Find optimal m')
      .add(dlgMMin)
      .add(dlgMMax)
      .onOK(() => {
        if (!validateDialog())
          return;
        runOptimization(dlgMMin.value!, dlgMMax.value!);
      })
      .show();
  }

  // --- Toolbar buttons ---
  optimizeBtn = ui.iconFA('search', () => showOptimizeDialog(), 'Find the m value that maximizes the occupied patch fraction at t_end');

  const resetBtn = ui.iconFA('undo', () => {
    computationsBlocked = true;
    ctrlP0.value = DEFAULTS.p0;
    ctrlM.value = DEFAULTS.m;
    ctrlE0.value = DEFAULTS.e0;
    ctrlRescue.value = DEFAULTS.rescueEffect;
    ctrlTStart.value = DEFAULTS.t_start;
    ctrlTEnd.value = DEFAULTS.t_end;
    ctrlTStep.value = DEFAULTS.t_step;
    ctrlTolerance.value = DEFAULTS.tolerance;
    computationsBlocked = false;
    updateRhoBadge();
    runPrimary();
  }, 'Reset all parameters to default values');

  view.setRibbonPanels([[optimizeBtn, resetBtn]]);

  // --- Layout: left panel with form ---
  const form = ui.form([]);

  form.append(ui.h2('Initial Condition'));
  form.append(ctrlP0.root);

  form.append(ui.h2('Parameters'));
  form.append(ctrlM.root);
  form.append(ctrlE0.root);
  form.append(ctrlRescue.root);
  form.append(rhoBadge);

  form.append(ui.h2('Argument'));
  form.append(ctrlTStart.root);
  form.append(ctrlTEnd.root);
  form.append(ctrlTStep.root);

  form.append(ui.h2('Solver'));
  form.append(ctrlTolerance.root);

  const dockMng = view.dockManager;
  dockMng.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.2);

  // --- Line chart ---
  lineChart = view.addViewer('Line chart', {
    xColumnName: 't',
    yColumnNames: ['p'],
    title: 'p(t) Dynamics',
  });

  const gridNode = dockMng.findNode(view.grid.root);
  if (gridNode != null)
    dockMng.dock(lineChart, DG.DOCK_TYPE.RIGHT, gridNode, undefined, 0.5);

  // --- Initial color coding and tooltip ---
  updateColorCoding();
  setupGridTooltip();

  // --- Cleanup on close ---
  subs.push(grok.events.onViewRemoved.subscribe((v: any) => {
    if (v === view) {
      terminateWorkers();
      for (const sub of subs)
        sub.unsubscribe();
      if (debounceTimer !== null)
        clearTimeout(debounceTimer);
    }
  }));
}

// Package reference (set from package.ts)
let _package: DG.Package;
export function setPackage(pkg: DG.Package): void {
  _package = pkg;
}
