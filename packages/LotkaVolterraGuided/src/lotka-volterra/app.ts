// Lotka-Volterra Predator-Prey Simulation — Application (Coordinator + UI)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULTS, RANGES, validate, solve,
  LotkaVolterraParams, LotkaVolterraSolution, InputId, WorkerTask, WorkerResult,
} from './core';

import '../../css/lotka-volterra.css';

const DEBOUNCE_MS = 50;
const GRID_STEPS = 11; // 0%, 10%, ..., 100% of range
const TOTAL_GRID_POINTS = GRID_STEPS ** 4; // 14641

export function lotkaVolterraApp(_package: DG.Package): void {
  // --- State ---
  let computationsBlocked = false;
  let debounceTimer: ReturnType<typeof setTimeout> | null = null;
  const subs: {unsubscribe(): void}[] = [];
  let activeWorkers: Worker[] = [];
  let lineChart!: DG.Viewer;
  let scatterPlot!: DG.Viewer;

  // --- Initial DataFrame ---
  const initSolution = solve(DEFAULTS);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromFloat64Array('t', initSolution.t),
    DG.Column.fromFloat64Array('x', initSolution.x),
    DG.Column.fromFloat64Array('y', initSolution.y),
  ]);
  df.name = 'Lotka-Volterra';

  const view = grok.shell.addTableView(df);
  view.name = 'Lotka-Volterra Predator-Prey';

  // --- Stats panel ---
  const statsPanel = ui.div([], 'lv-stats-panel');

  function updateStatsPanel(sol: LotkaVolterraSolution): void {
    statsPanel.innerHTML = '';
    const lines = [
      `<span class="lv-stats-label">Equilibrium:</span>`,
      `<span class="lv-stats-value">  x* = ${sol.xStar.toFixed(2)}, y* = ${sol.yStar.toFixed(2)}</span>`,
      `<span class="lv-stats-label">Max prey:</span> <span class="lv-stats-value">${sol.maxPrey.toFixed(2)}</span>`,
      `<span class="lv-stats-label">Max predators:</span> <span class="lv-stats-value">${sol.maxPredators.toFixed(2)}</span>`,
      `<span class="lv-stats-label">Steps:</span> <span class="lv-stats-value">${sol.stepCount}</span>`,
    ];
    statsPanel.innerHTML = lines.join('<br>');
  }
  updateStatsPanel(initSolution);

  // --- Controls ---

  // Model Coefficients
  const ctrlAlpha = ui.input.float('Prey birth rate \u03B1', {
    value: DEFAULTS.alpha, nullable: false,
    min: RANGES.alpha.min, max: RANGES.alpha.max,
    tooltipText: 'Rate at which prey reproduce. Higher \u03B1 \u2192 faster prey growth in the absence of predators.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlBeta = ui.input.float('Predation rate \u03B2', {
    value: DEFAULTS.beta, nullable: false,
    min: RANGES.beta.min, max: RANGES.beta.max,
    tooltipText: 'Rate at which predators consume prey. Higher \u03B2 \u2192 more prey eaten per encounter, reducing prey population faster.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlDelta = ui.input.float('Predator efficiency \u03B4', {
    value: DEFAULTS.delta, nullable: false,
    min: RANGES.delta.min, max: RANGES.delta.max,
    tooltipText: 'Efficiency of converting consumed prey into predator growth. Higher \u03B4 \u2192 predators grow faster from each prey consumed.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlGamma = ui.input.float('Predator death rate \u03B3', {
    value: DEFAULTS.gamma, nullable: false,
    min: RANGES.gamma.min, max: RANGES.gamma.max,
    tooltipText: 'Natural death rate of predators. Higher \u03B3 \u2192 predators die off faster without sufficient prey.',
    onValueChanged: () => debouncedRun(),
  });

  // Initial conditions
  const ctrlX0 = ui.input.float('Initial prey x\u2080', {
    value: DEFAULTS.x0, nullable: false,
    min: RANGES.x0.min, max: RANGES.x0.max,
    tooltipText: 'Starting prey population at time t=0.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlY0 = ui.input.float('Initial predators y\u2080', {
    value: DEFAULTS.y0, nullable: false,
    min: RANGES.y0.min, max: RANGES.y0.max,
    tooltipText: 'Starting predator population at time t=0.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlT = ui.input.float('Simulation time T', {
    value: DEFAULTS.T, nullable: false,
    min: RANGES.T.min, max: RANGES.T.max,
    tooltipText: 'Total simulation time. Longer T shows more oscillation cycles.',
    onValueChanged: () => debouncedRun(),
  });

  // Set formats (block computations to avoid spurious runs from format-triggered events)
  computationsBlocked = true;
  ctrlAlpha.format = '0.00';
  ctrlBeta.format = '0.000';
  ctrlDelta.format = '0.000';
  ctrlGamma.format = '0.00';
  ctrlX0.format = '0.0';
  ctrlY0.format = '0.0';
  ctrlT.format = '0.0';
  computationsBlocked = false;

  // --- Input map for validators ---
  const inputMap: Record<InputId, DG.InputBase> = {
    'ctrl_alpha': ctrlAlpha,
    'ctrl_beta': ctrlBeta,
    'ctrl_delta': ctrlDelta,
    'ctrl_gamma': ctrlGamma,
    'ctrl_x0': ctrlX0,
    'ctrl_y0': ctrlY0,
    'ctrl_T': ctrlT,
  };

  // --- Gather current inputs ---
  function getInputs(): LotkaVolterraParams {
    return {
      alpha: ctrlAlpha.value ?? DEFAULTS.alpha,
      beta: ctrlBeta.value ?? DEFAULTS.beta,
      delta: ctrlDelta.value ?? DEFAULTS.delta,
      gamma: ctrlGamma.value ?? DEFAULTS.gamma,
      x0: ctrlX0.value ?? DEFAULTS.x0,
      y0: ctrlY0.value ?? DEFAULTS.y0,
      T: ctrlT.value ?? DEFAULTS.T,
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

    ctrlAlpha.addValidator(validatorFor('ctrl_alpha'));
    ctrlBeta.addValidator(validatorFor('ctrl_beta'));
    ctrlDelta.addValidator(validatorFor('ctrl_delta'));
    ctrlGamma.addValidator(validatorFor('ctrl_gamma'));
    ctrlX0.addValidator(validatorFor('ctrl_x0'));
    ctrlY0.addValidator(validatorFor('ctrl_y0'));
    ctrlT.addValidator(validatorFor('ctrl_T'));
  }
  addValidators();

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
      updateStatsPanel(result);
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
  function updateDataFrame(result: LotkaVolterraSolution): void {
    const newDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('t', result.t),
      DG.Column.fromFloat64Array('x', result.x),
      DG.Column.fromFloat64Array('y', result.y),
    ]);
    newDf.name = 'Lotka-Volterra';
    view.dataFrame = newDf;
    lineChart.dataFrame = newDf;
    scatterPlot.dataFrame = newDf;
  }

  function clearResults(): void {
    const emptyDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('t', new Float64Array(0)),
      DG.Column.fromFloat64Array('x', new Float64Array(0)),
      DG.Column.fromFloat64Array('y', new Float64Array(0)),
    ]);
    emptyDf.name = 'Lotka-Volterra';
    view.dataFrame = emptyDf;
    lineChart.dataFrame = emptyDf;
    scatterPlot.dataFrame = emptyDf;
  }

  // --- Optimization task ---
  const optimizeBtn = ui.div([], 'lv-optimize-btn');
  const progressBar = ui.div([], 'lv-progress-bar');
  optimizeBtn.appendChild(progressBar);
  const optimizeBtnLabel = ui.div([], '');
  optimizeBtnLabel.textContent = 'Optimize Max Prey';
  optimizeBtn.appendChild(optimizeBtnLabel);
  ui.tooltip.bind(optimizeBtn, 'Run brute-force grid search over all four model coefficients to maximize peak prey population');

  function setOptimizeEnabled(enabled: boolean): void {
    optimizeBtn.classList.toggle('lv-optimize-btn--disabled', !enabled);
  }

  function setOptimizeProgress(percent: number): void {
    progressBar.style.width = `${percent}%`;
    optimizeBtnLabel.textContent = percent < 100 ? `Optimizing... ${percent}%` : 'Optimize Max Prey';
  }

  optimizeBtn.addEventListener('click', () => {
    runOptimization();
  });

  async function runOptimization(): Promise<void> {
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Invalid parameters. Fix inputs before optimizing.');
      return;
    }

    setOptimizeEnabled(false);
    setOptimizeProgress(0);

    let completed = 0;
    let errorCount = 0;
    const results: WorkerResult[] = [];
    const workerCount = Math.max(1, (navigator.hardwareConcurrency ?? 4) - 2);

    // Generate grid points: 11 steps per parameter (10% of range each)
    const alphaValues = linspace(RANGES.alpha.min, RANGES.alpha.max, GRID_STEPS);
    const betaValues = linspace(RANGES.beta.min, RANGES.beta.max, GRID_STEPS);
    const deltaValues = linspace(RANGES.delta.min, RANGES.delta.max, GRID_STEPS);
    const gammaValues = linspace(RANGES.gamma.min, RANGES.gamma.max, GRID_STEPS);

    const tasks: WorkerTask[] = [];
    for (const a of alphaValues) {
      for (const b of betaValues) {
        for (const d of deltaValues) {
          for (const g of gammaValues) {
            tasks.push({alpha: a, beta: b, delta: d, gamma: g, x0: inputs.x0, y0: inputs.y0, T: inputs.T});
          }
        }
      }
    }

    const workerUrl = _package.webRoot + 'dist/optimize-worker.js';

    try {
      await new Promise<void>((resolve, reject) => {
        const taskQueue = [...tasks];
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
          if (taskQueue.length === 0) {
            worker.terminate();
            activeWorkers = activeWorkers.filter((w) => w !== worker);
            if (activeWorkers.length === 0)
              resolve();
            return;
          }

          const task = taskQueue.shift()!;
          worker.postMessage(task);
        };

        const handleResult = (worker: Worker, event: MessageEvent<WorkerResult>) => {
          const result = event.data;
          completed++;
          setOptimizeProgress(Math.round(completed / TOTAL_GRID_POINTS * 100));

          if (result.error)
            errorCount++;
          else
            results.push(result);

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
            setOptimizeProgress(Math.round(completed / TOTAL_GRID_POINTS * 100));
            processNext(worker);
          };

          processNext(worker);
        }
      });
    } catch (err) {
      grok.shell.error(err instanceof Error ? err.message : 'Failed to start parallel computations.');
      setOptimizeProgress(0);
      setOptimizeEnabled(true);
      return;
    }

    terminateWorkers();
    setOptimizeProgress(0);
    setOptimizeEnabled(true);

    // Handle results
    if (errorCount > 0 && errorCount < TOTAL_GRID_POINTS)
      grok.shell.warning(`${errorCount} of ${TOTAL_GRID_POINTS} points failed. Result based on ${TOTAL_GRID_POINTS - errorCount} points.`);

    if (results.length === 0) {
      grok.shell.error('Failed to compute any point. Check the parameters.');
      return;
    }

    // Find optimal
    let best = results[0];
    for (const r of results) {
      if (r.maxPrey > best.maxPrey)
        best = r;
    }

    // Batch update: block primary, write optimal values, unblock and run once
    try {
      computationsBlocked = true;
      ctrlAlpha.value = best.alpha;
      ctrlBeta.value = best.beta;
      ctrlDelta.value = best.delta;
      ctrlGamma.value = best.gamma;
      computationsBlocked = false;
      runPrimary();
    } catch (_err) {
      computationsBlocked = false;
      grok.shell.warning(`Optimal values found but could not update sliders automatically.`);
    }
  }

  function terminateWorkers(): void {
    for (const w of activeWorkers)
      w.terminate();
    activeWorkers = [];
  }

  function linspace(min: number, max: number, steps: number): number[] {
    const arr: number[] = [];
    for (let i = 0; i < steps; i++)
      arr.push(min + i * (max - min) / (steps - 1));
    return arr;
  }

  // --- Toolbar buttons ---
  const resetBtn = ui.iconFA('undo', () => {
    computationsBlocked = true;
    ctrlAlpha.value = DEFAULTS.alpha;
    ctrlBeta.value = DEFAULTS.beta;
    ctrlDelta.value = DEFAULTS.delta;
    ctrlGamma.value = DEFAULTS.gamma;
    ctrlX0.value = DEFAULTS.x0;
    ctrlY0.value = DEFAULTS.y0;
    ctrlT.value = DEFAULTS.T;
    computationsBlocked = false;
    runPrimary();
  }, 'Reset all parameters to default values');

  view.setRibbonPanels([[resetBtn]]);

  // --- Layout: left panel with form ---
  const form = ui.form([]);

  form.append(ui.h2('Model Coefficients'));
  form.append(ctrlAlpha.root);
  form.append(ctrlBeta.root);
  form.append(ctrlDelta.root);
  form.append(ctrlGamma.root);

  form.append(ui.h2('Initial Conditions'));
  form.append(ctrlX0.root);
  form.append(ctrlY0.root);
  form.append(ctrlT.root);

  form.append(optimizeBtn);

  form.append(ui.h2('Equilibrium & Stats'));
  form.append(statsPanel);

  const dockMng = view.dockManager;
  dockMng.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.2);

  // --- Line chart (time series) ---
  lineChart = view.addViewer('Line chart', {
    xColumnName: 't',
    yColumnNames: ['x', 'y'],
    title: 'Population Dynamics',
  });

  // --- Phase portrait (scatter plot) ---
  scatterPlot = view.addViewer('Scatter plot', {
    xColumnName: 'x',
    yColumnName: 'y',
    title: 'Phase Portrait',
    markerDefaultSize: 2,
  });

  // Dock: line chart top-center, phase portrait bottom-center, grid right
  const gridNode = dockMng.findNode(view.grid.root);
  if (gridNode != null) {
    dockMng.dock(lineChart, DG.DOCK_TYPE.LEFT, gridNode, undefined, 0.65);
    const lineChartNode = dockMng.findNode(lineChart.root);
    if (lineChartNode != null)
      dockMng.dock(scatterPlot, DG.DOCK_TYPE.DOWN, lineChartNode, undefined, 0.5);
  }

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
