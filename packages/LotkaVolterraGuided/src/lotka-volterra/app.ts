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
    DG.Column.fromFloat64Array('Time', initSolution.t),
    DG.Column.fromFloat64Array('Prey', initSolution.x),
    DG.Column.fromFloat64Array('Predators', initSolution.y),
  ]);
  df.name = 'Lotka-Volterra';
  df.col('Prey')!.meta.colors.setLinear([DG.Color.darkRed, DG.Color.yellow, DG.Color.darkGreen]);
  df.col('Predators')!.meta.colors.setLinear([DG.Color.darkRed, DG.Color.yellow, DG.Color.darkGreen]);

  const view = grok.shell.addTableView(df);
  view.name = 'Lotka-Volterra Predator-Prey';

  // --- Stats panel ---
  const eqLabel = ui.label('');
  const maxPreyLabel = ui.label('');
  const maxPredatorsLabel = ui.label('');
  const stepsLabel = ui.label('');

  const statsPanel = ui.divV([
    ui.divH([ui.label('Equilibrium:'), eqLabel], 'lotka-volterra-app-stats-row'),
    ui.divH([ui.label('Max prey:'), maxPreyLabel], 'lotka-volterra-app-stats-row'),
    ui.divH([ui.label('Max predators:'), maxPredatorsLabel], 'lotka-volterra-app-stats-row'),
    ui.divH([ui.label('Steps:'), stepsLabel], 'lotka-volterra-app-stats-row'),
  ], 'lotka-volterra-app-stats-panel');

  function updateStatsPanel(sol: LotkaVolterraSolution): void {
    eqLabel.textContent = `x* = ${sol.xStar.toFixed(2)}, y* = ${sol.yStar.toFixed(2)}`;
    maxPreyLabel.textContent = sol.maxPrey.toFixed(2);
    maxPredatorsLabel.textContent = sol.maxPredators.toFixed(2);
    stepsLabel.textContent = `${sol.stepCount}`;
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
  function applyColorCoding(targetDf: DG.DataFrame): void {
    targetDf.col('Prey')!.meta.colors.setLinear([DG.Color.darkRed, DG.Color.yellow, DG.Color.darkGreen]);
    targetDf.col('Predators')!.meta.colors.setLinear([DG.Color.darkRed, DG.Color.yellow, DG.Color.darkGreen]);
  }

  function updateDataFrame(result: LotkaVolterraSolution): void {
    const newDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Time', result.t),
      DG.Column.fromFloat64Array('Prey', result.x),
      DG.Column.fromFloat64Array('Predators', result.y),
    ]);
    newDf.name = 'Lotka-Volterra';
    applyColorCoding(newDf);
    view.dataFrame = newDf;
    lineChart.dataFrame = newDf;
    scatterPlot.dataFrame = newDf;
  }

  function clearResults(): void {
    const emptyDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Time', new Float64Array(0)),
      DG.Column.fromFloat64Array('Prey', new Float64Array(0)),
      DG.Column.fromFloat64Array('Predators', new Float64Array(0)),
    ]);
    emptyDf.name = 'Lotka-Volterra';
    view.dataFrame = emptyDf;
    lineChart.dataFrame = emptyDf;
    scatterPlot.dataFrame = emptyDf;
  }

  // --- Optimization task ---
  const optimizeBtn = ui.bigButton('Optimize', () => runOptimization(),
    'Run brute-force grid search over all four model coefficients to maximize peak prey population');

  function setOptimizeEnabled(enabled: boolean): void {
    optimizeBtn.classList.toggle('lotka-volterra-app-optimize-btn--disabled', !enabled);
  }

  async function runOptimization(): Promise<void> {
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Invalid parameters. Fix inputs before optimizing.');
      return;
    }

    setOptimizeEnabled(false);

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
    const nWorkers = Math.min(workerCount, tasks.length);
    const chunks: WorkerTask[][] = Array.from({length: nWorkers}, () => []);
    for (let i = 0; i < tasks.length; i++)
      chunks[i % nWorkers].push(tasks[i]);

    activeWorkers = [];
    const pi = DG.TaskBarProgressIndicator.create('Optimizing Max Prey...', {cancelable: true});

    const resolvers = new Array<(value: WorkerResult[]) => void>(chunks.length);
    const batchPromises = chunks.map((batch, i) =>
      new Promise<WorkerResult[]>((resolve, reject) => {
        resolvers[i] = resolve;
        let worker: Worker;
        try {
          worker = new Worker(workerUrl);
        } catch (_err) {
          reject(new Error('Failed to start parallel computations. Try again later.'));
          return;
        }
        activeWorkers.push(worker);

        worker.onmessage = (event: MessageEvent<WorkerResult[]>) => {
          worker.terminate();
          activeWorkers = activeWorkers.filter((w) => w !== worker);
          resolve(event.data);
        };

        worker.onerror = (err) => {
          worker.terminate();
          activeWorkers = activeWorkers.filter((w) => w !== worker);
          reject(new Error(err.message ?? 'Worker error'));
        };

        worker.postMessage(batch);
      }),
    );

    const cancelSub = pi.onCanceled.subscribe(() => {
      terminateWorkers();
      for (const resolve of resolvers)
        resolve([]);
    });

    const settled = await Promise.allSettled(batchPromises);

    cancelSub.unsubscribe();
    pi.close();
    terminateWorkers();
    setOptimizeEnabled(true);

    if (pi.canceled)
      return;

    const results: WorkerResult[] = [];
    let errorCount = 0;
    for (let i = 0; i < settled.length; i++) {
      const outcome = settled[i];
      if (outcome.status === 'fulfilled') {
        for (const r of outcome.value) {
          if (r.error)
            errorCount++;
          else
            results.push(r);
        }
      } else {
        errorCount += chunks[i].length;
      }
    }

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

  const helpBtn = ui.icons.help(
    () => window.open('https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations', '_blank'),
    'Learn more about the Predator-Prey model. Opens in a separate page.',
  );
  helpBtn.classList.add('lotka-volterra-app-help-icon');

  view.setRibbonPanels([[resetBtn], [helpBtn]]);

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
  dockMng.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.1);

  // --- Line chart (time series) ---
  lineChart = view.addViewer('Line chart', {
    xColumnName: 'Time',
    yColumnNames: ['Prey', 'Predators'],
    title: 'Population Dynamics',
  });

  // --- Phase portrait (scatter plot) ---
  scatterPlot = view.addViewer('Scatter plot', {
    xColumnName: 'Prey',
    yColumnName: 'Predators',
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

  // --- Grid column header tooltips ---
  function columnTooltip(colName: string): HTMLElement {
    const col = view.dataFrame.col(colName);
    const minVal = col?.min ?? 0;
    const maxVal = col?.max ?? 0;

    const minSpan = ui.label(`min ${minVal.toFixed(2)}`);
    minSpan.classList.add('lotka-volterra-app-tooltip-min');

    const maxSpan = ui.label(`max ${maxVal.toFixed(2)}`);
    maxSpan.classList.add('lotka-volterra-app-tooltip-max');

    return ui.divV([
      ui.h2(colName),
      ui.div([minSpan, document.createTextNode(' ... '), maxSpan], 'lotka-volterra-app-tooltip-range'),
    ], 'lotka-volterra-app-col-tooltip');
  }

  view.grid.onCellTooltip((cell, x, y) => {
    if (!cell.isColHeader)
      return false;

    const colName = cell.tableColumn?.name;
    if (colName === 'Prey' || colName === 'Predators') {
      ui.tooltip.show(columnTooltip(colName), x, y);
      return true;
    }
    return false;
  });

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
