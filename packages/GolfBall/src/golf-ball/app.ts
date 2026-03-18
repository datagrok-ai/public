// Golf Ball Flight Simulator — Application (Coordinator + UI)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULTS, RANGES, validate, solve, linspace,
  GolfBallParams, GolfBallSolution, InputId,
  OptimizeParam, OPTIMIZE_PARAMS, OptimizeTask, OptimizeResult,
  SweepParamName, SWEEP_PARAMS, SensitivityTask, SensitivityResult,
} from './core';

import '../../css/golf-ball.css';

const DEBOUNCE_MS = 50;
const OPTIMIZE_STEPS = 1000;
const SENSITIVITY_STEPS = 100;

export function golfBallApp(_package: DG.Package): void {
  // --- State ---
  let computationsBlocked = false;
  let debounceTimer: ReturnType<typeof setTimeout> | null = null;
  const subs: {unsubscribe(): void}[] = [];
  let activeWorkers: Worker[] = [];
  let trajectoryPlot!: DG.Viewer;
  let timeSeriesChart!: DG.Viewer;

  // --- Initial DataFrame ---
  const initSolution = solve(DEFAULTS);
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromFloat64Array('Time', initSolution.t),
    DG.Column.fromFloat64Array('x', initSolution.x),
    DG.Column.fromFloat64Array('y', initSolution.y),
    DG.Column.fromFloat64Array('v', initSolution.v),
  ]);
  df.name = 'Golf Ball Flight';

  const view = grok.shell.addTableView(df);
  view.name = 'Golf Ball Flight Simulator';

  // --- Stats panel ---
  const maxHeightLabel = ui.label('');
  const maxHeightTimeLabel = ui.label('');
  const flightDistLabel = ui.label('');
  const flightTimeLabel = ui.label('');
  const landingSpeedLabel = ui.label('');
  const landingAngleLabel = ui.label('');

  maxHeightLabel.classList.add('golf-ball-app-stats-value');
  maxHeightTimeLabel.classList.add('golf-ball-app-stats-value');
  flightDistLabel.classList.add('golf-ball-app-stats-value');
  flightTimeLabel.classList.add('golf-ball-app-stats-value');
  landingSpeedLabel.classList.add('golf-ball-app-stats-value');
  landingAngleLabel.classList.add('golf-ball-app-stats-value');

  const statsPanel = ui.divV([
    ui.label('Height'),
    maxHeightLabel,
    maxHeightTimeLabel,
    ui.label('Distance & Time'),
    flightDistLabel,
    flightTimeLabel,
    ui.label('Landing'),
    landingSpeedLabel,
    landingAngleLabel,
  ], 'golf-ball-app-stats-panel');

  function updateStatsPanel(sol: GolfBallSolution): void {
    maxHeightLabel.textContent = `Max height = ${sol.maxHeight.toFixed(2)} m`;
    maxHeightTimeLabel.textContent = `at t = ${sol.maxHeightTime.toFixed(2)} s`;
    flightDistLabel.textContent = `Distance = ${sol.flightDistance.toFixed(2)} m`;
    flightTimeLabel.textContent = `Flight time = ${sol.flightTime.toFixed(2)} s`;
    landingSpeedLabel.textContent = `Speed = ${sol.landingSpeed.toFixed(2)} m/s`;
    landingAngleLabel.textContent = `Angle = ${sol.landingAngle.toFixed(1)}\u00B0`;
  }
  updateStatsPanel(initSolution);

  // --- Controls ---

  // Launch parameters
  const ctrlV0 = ui.input.float('Initial speed v\u2080', {
    value: DEFAULTS.v0, nullable: false,
    min: RANGES.v0.min, max: RANGES.v0.max,
    tooltipText: 'Initial ball speed off the clubface. A typical driver produces 50\u201375 m/s.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlTheta = ui.input.float('Launch angle \u03B8', {
    value: DEFAULTS.theta, nullable: false,
    min: RANGES.theta.min, max: RANGES.theta.max,
    tooltipText: 'Launch angle relative to the ground. Optimal for distance is usually 10\u201315\u00B0 with drag, ~45\u00B0 without.',
    onValueChanged: () => debouncedRun(),
  });

  // Ball parameters
  const ctrlM = ui.input.float('Ball mass m', {
    value: DEFAULTS.m, nullable: false,
    min: RANGES.m.min, max: RANGES.m.max,
    tooltipText: 'Mass of the golf ball. Regulation balls weigh no more than 45.93 g (0.0459 kg).',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlD = ui.input.float('Ball diameter d', {
    value: DEFAULTS.d, nullable: false,
    min: RANGES.d.min, max: RANGES.d.max,
    tooltipText: 'Diameter of the ball. Regulation minimum is 42.67 mm (0.0427 m).',
    onValueChanged: () => debouncedRun(),
  });

  // Environment parameters
  const ctrlCd = ui.input.float('Drag coefficient Cd', {
    value: DEFAULTS.Cd, nullable: false,
    min: RANGES.Cd.min, max: RANGES.Cd.max,
    tooltipText: 'Aerodynamic drag coefficient. A smooth ball is ~0.5, a dimpled golf ball ~0.25.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlRho = ui.input.float('Air density \u03C1', {
    value: DEFAULTS.rho, nullable: false,
    min: RANGES.rho.min, max: RANGES.rho.max,
    tooltipText: 'Air density. Sea level \u2248 1.225 kg/m\u00B3. Higher altitude or warmer air \u2192 lower density \u2192 longer flight.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlG = ui.input.float('Gravity g', {
    value: DEFAULTS.g, nullable: false,
    min: RANGES.g.min, max: RANGES.g.max,
    tooltipText: 'Gravitational acceleration. Earth \u2248 9.81 m/s\u00B2. Moon \u2248 1.62, Mars \u2248 3.72, Jupiter \u2248 24.79.',
    onValueChanged: () => debouncedRun(),
  });

  // Set formats (block computations to avoid spurious runs)
  computationsBlocked = true;
  ctrlV0.format = '0.0';
  ctrlTheta.format = '0.0';
  ctrlM.format = '0.0000';
  ctrlD.format = '0.0000';
  ctrlCd.format = '0.00';
  ctrlRho.format = '0.000';
  ctrlG.format = '0.00';
  computationsBlocked = false;

  // --- Input map for validators ---
  const inputMap: Record<InputId, DG.InputBase> = {
    'ctrl_v0': ctrlV0,
    'ctrl_theta': ctrlTheta,
    'ctrl_m': ctrlM,
    'ctrl_d': ctrlD,
    'ctrl_Cd': ctrlCd,
    'ctrl_rho': ctrlRho,
    'ctrl_g': ctrlG,
  };

  // --- Gather current inputs ---
  function getInputs(): GolfBallParams {
    return {
      v0: ctrlV0.value ?? DEFAULTS.v0,
      theta: ctrlTheta.value ?? DEFAULTS.theta,
      m: ctrlM.value ?? DEFAULTS.m,
      d: ctrlD.value ?? DEFAULTS.d,
      Cd: ctrlCd.value ?? DEFAULTS.Cd,
      rho: ctrlRho.value ?? DEFAULTS.rho,
      g: ctrlG.value ?? DEFAULTS.g,
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

    ctrlV0.addValidator(validatorFor('ctrl_v0'));
    ctrlTheta.addValidator(validatorFor('ctrl_theta'));
    ctrlM.addValidator(validatorFor('ctrl_m'));
    ctrlD.addValidator(validatorFor('ctrl_d'));
    ctrlCd.addValidator(validatorFor('ctrl_Cd'));
    ctrlRho.addValidator(validatorFor('ctrl_rho'));
    ctrlG.addValidator(validatorFor('ctrl_g'));
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
  function updateDataFrame(result: GolfBallSolution): void {
    const newDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Time', result.t),
      DG.Column.fromFloat64Array('x', result.x),
      DG.Column.fromFloat64Array('y', result.y),
      DG.Column.fromFloat64Array('v', result.v),
    ]);
    newDf.name = 'Golf Ball Flight';
    view.dataFrame = newDf;
    trajectoryPlot.dataFrame = newDf;
    timeSeriesChart.dataFrame = newDf;
  }

  function clearResults(): void {
    const emptyDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Time', new Float64Array(0)),
      DG.Column.fromFloat64Array('x', new Float64Array(0)),
      DG.Column.fromFloat64Array('y', new Float64Array(0)),
      DG.Column.fromFloat64Array('v', new Float64Array(0)),
    ]);
    emptyDf.name = 'Golf Ball Flight';
    view.dataFrame = emptyDf;
    trajectoryPlot.dataFrame = emptyDf;
    timeSeriesChart.dataFrame = emptyDf;
  }

  // --- Maximize Height task ---
  const maximizeBtn = ui.bigButton('Maximize Height', () => showMaximizeDialog(),
    'Find the value of the selected parameter that produces the highest ball apex.');

  function setMaximizeEnabled(enabled: boolean): void {
    maximizeBtn.classList.toggle('golf-ball-app-btn--disabled', !enabled);
  }

  function showMaximizeDialog(): void {
    const choiceInput = ui.input.choice('Parameter to optimize', {
      items: OPTIMIZE_PARAMS as unknown as string[],
      value: OPTIMIZE_PARAMS[0],
      nullable: false,
      tooltipText: 'Choose which parameter to vary in the search for maximum ball height.',
    });

    const dlg = ui.dialog('Maximize Height')
      .add(choiceInput)
      .addButton('Run', () => { dlg.close(); runMaximizeHeight(choiceInput.value as OptimizeParam); })
      .show();

    ui.tooltip.bind(dlg.getButton('Run')!,
      'Run brute-force grid search to find the parameter value that maximizes peak altitude.');
  }

  async function runMaximizeHeight(paramName: OptimizeParam): Promise<void> {
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Invalid parameters. Fix inputs before optimizing.');
      return;
    }

    setMaximizeEnabled(false);

    const range = RANGES[paramName];
    const allValues = linspace(range.min, range.max, OPTIMIZE_STEPS);

    const workerCount = Math.max(1, (navigator.hardwareConcurrency ?? 4) - 2);
    const workerUrl = _package.webRoot + 'dist/optimize-worker.js';
    const nWorkers = Math.min(workerCount, OPTIMIZE_STEPS);

    // Distribute values round-robin
    const chunks: number[][] = Array.from({length: nWorkers}, () => []);
    for (let i = 0; i < allValues.length; i++)
      chunks[i % nWorkers].push(allValues[i]);

    activeWorkers = [];
    const pi = DG.TaskBarProgressIndicator.create('Maximizing height...', {cancelable: true});

    const resolvers = new Array<(value: OptimizeResult[]) => void>(nWorkers);
    const batchPromises = chunks.map((chunk, i) =>
      new Promise<OptimizeResult[]>((resolve, reject) => {
        resolvers[i] = resolve;
        let worker: Worker;
        try {
          worker = new Worker(workerUrl);
        } catch (_err) {
          reject(new Error('Failed to start parallel computations. Try again later.'));
          return;
        }
        activeWorkers.push(worker);

        worker.onmessage = (event: MessageEvent<OptimizeResult[]>) => {
          worker.terminate();
          activeWorkers = activeWorkers.filter((w) => w !== worker);
          resolve(event.data);
        };

        worker.onerror = (err) => {
          worker.terminate();
          activeWorkers = activeWorkers.filter((w) => w !== worker);
          reject(new Error(err.message ?? 'Worker error'));
        };

        const task: OptimizeTask = {paramName, values: chunk, baseParams: inputs};
        worker.postMessage(task);
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
    setMaximizeEnabled(true);

    if (pi.canceled)
      return;

    const results: OptimizeResult[] = [];
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

    if (errorCount > 0 && errorCount < OPTIMIZE_STEPS)
      grok.shell.warning(`${errorCount} of ${OPTIMIZE_STEPS} points failed. Result based on ${OPTIMIZE_STEPS - errorCount} points.`);

    if (results.length === 0) {
      grok.shell.error('Failed to compute any point. Check the parameters.');
      return;
    }

    // Find optimal
    let best = results[0];
    for (const r of results) {
      if (r.maxHeight > best.maxHeight)
        best = r;
    }

    // Batch update: block primary, write optimal value, unblock and run once
    const controlMap: Record<OptimizeParam, DG.InputBase> = {
      'theta': ctrlTheta,
      'v0': ctrlV0,
      'Cd': ctrlCd,
    };

    try {
      computationsBlocked = true;
      controlMap[paramName].value = best.value;
      computationsBlocked = false;
      runPrimary();
    } catch (_err) {
      computationsBlocked = false;
      grok.shell.warning('Optimal value found but could not update control automatically.');
    }
  }

  // --- Sensitivity Analysis task ---
  function setControlsEnabled(enabled: boolean): void {
    for (const input of Object.values(inputMap))
      input.enabled = enabled;
    maximizeBtn.classList.toggle('golf-ball-app-btn--disabled', !enabled);
  }

  async function runSensitivity(): Promise<void> {
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Invalid parameters. Fix inputs before running sensitivity analysis.');
      return;
    }

    setControlsEnabled(false);

    const workerUrl = _package.webRoot + 'dist/sensitivity-worker.js';
    const pi = DG.TaskBarProgressIndicator.create('Running sensitivity analysis...', {cancelable: true});

    activeWorkers = [];
    const resolvers = new Array<(value: SensitivityResult) => void>(SWEEP_PARAMS.length);

    const sweepPromises = SWEEP_PARAMS.map((paramName, i) =>
      new Promise<SensitivityResult>((resolve, reject) => {
        resolvers[i] = resolve;
        let worker: Worker;
        try {
          worker = new Worker(workerUrl);
        } catch (_err) {
          reject(new Error('Failed to start sensitivity worker.'));
          return;
        }
        activeWorkers.push(worker);

        worker.onmessage = (event: MessageEvent<SensitivityResult>) => {
          worker.terminate();
          activeWorkers = activeWorkers.filter((w) => w !== worker);
          resolve(event.data);
        };

        worker.onerror = (err) => {
          worker.terminate();
          activeWorkers = activeWorkers.filter((w) => w !== worker);
          reject(new Error(err.message ?? 'Worker error'));
        };

        const task: SensitivityTask = {
          paramName,
          paramMin: RANGES[paramName].min,
          paramMax: RANGES[paramName].max,
          steps: SENSITIVITY_STEPS,
          baseParams: inputs,
        };
        worker.postMessage(task);
      }),
    );

    const cancelSub = pi.onCanceled.subscribe(() => {
      terminateWorkers();
      for (const resolve of resolvers)
        resolve({paramName: '', paramValues: [], distances: [], error: 'canceled'});
    });

    const settled = await Promise.allSettled(sweepPromises);

    cancelSub.unsubscribe();
    pi.close();
    terminateWorkers();
    setControlsEnabled(true);

    if (pi.canceled)
      return;

    // Collect results
    const successResults: SensitivityResult[] = [];
    let errorCount = 0;
    for (const outcome of settled) {
      if (outcome.status === 'fulfilled' && !outcome.value.error)
        successResults.push(outcome.value);
      else
        errorCount++;
    }

    if (errorCount > 0 && successResults.length > 0)
      grok.shell.warning(`${errorCount} of ${SWEEP_PARAMS.length} parameter sweeps failed.`);

    if (successResults.length === 0) {
      grok.shell.error('Sensitivity analysis failed for all parameters.');
      return;
    }

    // Build DataFrame with interleaved columns: [param_values, param_distance, ...]
    const columns: DG.Column[] = [];
    for (const r of successResults) {
      columns.push(DG.Column.fromFloat64Array(r.paramName, new Float64Array(r.paramValues)));
      columns.push(DG.Column.fromFloat64Array(`Distance(${r.paramName})`, new Float64Array(r.distances)));
    }

    const sensDf = DG.DataFrame.fromColumns(columns);
    sensDf.name = 'Sensitivity Analysis';

    const sensView = grok.shell.addTableView(sensDf);
    sensView.name = 'Sensitivity Analysis';

    // Add line charts for each parameter
    for (const r of successResults) {
      sensView.addViewer('Line chart', {
        xColumnName: r.paramName,
        yColumnNames: [`Distance(${r.paramName})`],
        title: `Distance vs ${r.paramName}`,
        lineWidth: 2,
      });
    }
  }

  function terminateWorkers(): void {
    for (const w of activeWorkers)
      w.terminate();
    activeWorkers = [];
  }

  // --- Toolbar buttons ---
  const resetBtn = ui.iconFA('undo', () => {
    computationsBlocked = true;
    ctrlV0.value = DEFAULTS.v0;
    ctrlTheta.value = DEFAULTS.theta;
    ctrlM.value = DEFAULTS.m;
    ctrlD.value = DEFAULTS.d;
    ctrlCd.value = DEFAULTS.Cd;
    ctrlRho.value = DEFAULTS.rho;
    ctrlG.value = DEFAULTS.g;
    computationsBlocked = false;
    runPrimary();
  }, 'Reset all parameters to default values');

  const sensitivityBtn = ui.iconFA('chart-line', () => runSensitivity(),
    'Show how each parameter independently affects the total flight distance.');

  const helpBtn = ui.icons.help(
    () => window.open('https://en.wikipedia.org/wiki/Projectile_motion', '_blank'),
    'Learn more about projectile motion with drag. Opens in a separate page.',
  );
  helpBtn.classList.add('golf-ball-app-help-icon');

  view.setRibbonPanels([[resetBtn], [sensitivityBtn], [helpBtn]]);

  // --- Layout: left panel with form ---
  const form = ui.form([]);

  form.append(ui.h2('Launch'));
  form.append(ctrlV0.root);
  form.append(ctrlTheta.root);

  form.append(ui.h2('Ball'));
  form.append(ctrlM.root);
  form.append(ctrlD.root);

  form.append(ui.h2('Environment'));
  form.append(ctrlCd.root);
  form.append(ctrlRho.root);
  form.append(ctrlG.root);

  form.append(maximizeBtn);

  form.append(ui.h2('Flight Statistics'));
  form.append(statsPanel);

  const dockMng = view.dockManager;
  dockMng.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.6);

  // --- Trajectory chart (scatter plot: x vs y) ---
  trajectoryPlot = view.addViewer('Scatter plot', {
    xColumnName: 'x',
    yColumnName: 'y',
    title: 'Trajectory',
    markerDefaultSize: 2,
  });

  // --- Time series chart (line chart: Time vs x, y, v) ---
  timeSeriesChart = view.addViewer('Line chart', {
    xColumnName: 'Time',
    yColumnNames: ['x', 'y', 'v'],
    title: 'Time Series',
  });

  // Dock: trajectory top-right, time series bottom-right
  dockMng.dock(trajectoryPlot, DG.DOCK_TYPE.RIGHT, null, undefined, 0.7);
  const trajectoryNode = dockMng.findNode(trajectoryPlot.root);
  if (trajectoryNode != null)
    dockMng.dock(timeSeriesChart, DG.DOCK_TYPE.DOWN, trajectoryNode, undefined, 0.5);

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
