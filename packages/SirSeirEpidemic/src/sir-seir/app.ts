// SIR/SEIR Epidemic Simulation — Application (Coordinator + UI)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULTS, RANGES, validate, solve, POPULATION,
  EpidemicParams, EpidemicSolution, InputId,
} from './core';

import '../../css/sir-seir.css';

const DEBOUNCE_MS = 50;

export function sirSeirEpidemicApp(_package: DG.Package): void {
  // --- State ---
  let computationsBlocked = false;
  let debounceTimer: ReturnType<typeof setTimeout> | null = null;
  const subs: {unsubscribe(): void}[] = [];
  let lineChart!: DG.Viewer;
  let rEffChart!: DG.Viewer;
  let phaseChart!: DG.Viewer;

  // --- Initial DataFrame ---
  const initSolution = solve(DEFAULTS);
  const df = createDataFrame(initSolution, DEFAULTS.modelType);

  const view = grok.shell.addTableView(df);
  view.name = 'SIR / SEIR Epidemic Simulation';

  // --- Summary panel ---
  const peakLabel = ui.label('');
  const finalRecoveredLabel = ui.label('');
  const hitLabel = ui.label('');
  const r0Label = ui.label('');

  peakLabel.classList.add('sir-seir-app-summary-value');
  finalRecoveredLabel.classList.add('sir-seir-app-summary-value');
  hitLabel.classList.add('sir-seir-app-summary-value');
  r0Label.classList.add('sir-seir-app-summary-value');

  const summaryPanel = ui.divV([
    peakLabel,
    finalRecoveredLabel,
    hitLabel,
    r0Label,
  ], 'sir-seir-app-summary-panel');

  /** Updates summary panel labels from solution */
  function updateSummaryPanel(sol: EpidemicSolution, r0: number): void {
    peakLabel.textContent = `Peak infection: day ${sol.peakDay}, ${sol.peakCount} cases`;
    finalRecoveredLabel.textContent = `Final recovered: ${sol.finalRecoveredPct.toFixed(1)}% of population`;
    hitLabel.textContent = `Herd immunity threshold: ${sol.herdImmunityThreshold.toFixed(1)}% of population`;
    r0Label.textContent = `Basic reproduction number R\u2080 = ${r0.toFixed(1)}`;
  }
  updateSummaryPanel(initSolution, DEFAULTS.r0);

  // --- Controls ---

  // Model toggle (false = SIR, true = SEIR)
  const ctrlModelType = ui.input.toggle('Model', {
    value: false, nullable: false,
    tooltipText: 'Toggle between SIR (3 compartments) and SEIR (4 compartments, adds Exposed class with incubation period)',
    onValueChanged: () => {
      updateSigmaVisibility();
      debouncedRun();
    },
  });

  const ctrlR0 = ui.input.float('R\u2080 \u2014 basic reproduction number', {
    value: DEFAULTS.r0, nullable: false,
    min: RANGES.r0.min, max: RANGES.r0.max, step: 0.1,
    tooltipText: 'Basic reproduction number \u2014 average number of people one infected person infects in a fully susceptible population. Influenza \u2248 1.5, COVID-19 \u2248 2.5\u20133.5, Measles \u2248 12\u201318.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlGamma = ui.input.float('\u03B3 \u2014 recovery rate (1/\u03B3 = infectious period in days)', {
    value: DEFAULTS.gamma, nullable: false,
    min: RANGES.gamma.min, max: RANGES.gamma.max,
    tooltipText: '1/\u03B3 is the average duration of the infectious period in days. For example, \u03B3 = 0.1 means a person is infectious for ~10 days on average.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlSigma = ui.input.float('\u03C3 \u2014 incubation rate (1/\u03C3 = incubation period in days)', {
    value: DEFAULTS.sigma, nullable: false,
    min: RANGES.sigma.min, max: RANGES.sigma.max,
    tooltipText: '1/\u03C3 is the average incubation period in days. During this time the person is infected (E) but not yet infectious to others.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlVaccination = ui.input.float('Initial vaccination coverage (%)', {
    value: DEFAULTS.vaccination, nullable: false,
    min: RANGES.vaccination.min, max: RANGES.vaccination.max,
    tooltipText: 'Initial fraction of immune individuals. Herd immunity threshold = 1 \u2212 1/R\u2080. For R\u2080 = 3, at least 67% must be vaccinated to prevent an outbreak.',
    onValueChanged: () => debouncedRun(),
  });

  // Set formats (block computations to avoid spurious runs)
  computationsBlocked = true;
  ctrlR0.format = '0.0';
  ctrlGamma.format = '0.00';
  ctrlSigma.format = '0.00';
  ctrlVaccination.format = '0';
  computationsBlocked = false;

  // --- Sigma visibility ---
  function updateSigmaVisibility(): void {
    const isSEIR = ctrlModelType.value === true;
    ctrlSigma.root.classList.toggle('sir-seir-app-sigma-hidden', !isSEIR);
  }
  updateSigmaVisibility();

  // --- Input map for validators ---
  const inputMap: Record<InputId, DG.InputBase> = {
    'ctrl_r0': ctrlR0,
    'ctrl_gamma': ctrlGamma,
    'ctrl_sigma': ctrlSigma,
    'ctrl_vaccination': ctrlVaccination,
  };

  // --- Gather current inputs ---
  function getInputs(): EpidemicParams {
    return {
      modelType: ctrlModelType.value === true ? 'SEIR' : 'SIR',
      r0: ctrlR0.value ?? DEFAULTS.r0,
      gamma: ctrlGamma.value ?? DEFAULTS.gamma,
      sigma: ctrlSigma.value ?? DEFAULTS.sigma,
      vaccination: ctrlVaccination.value ?? DEFAULTS.vaccination,
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

    ctrlR0.addValidator(validatorFor('ctrl_r0'));
    ctrlGamma.addValidator(validatorFor('ctrl_gamma'));
    ctrlSigma.addValidator(validatorFor('ctrl_sigma'));
    ctrlVaccination.addValidator(validatorFor('ctrl_vaccination'));
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

      // Check for NaN/Inf
      if (hasNumericalIssue(result)) {
        clearResults();
        grok.shell.error('Numerical instability detected');
        return;
      }

      updateDataFrame(result, inputs.modelType);
      updateSummaryPanel(result, inputs.r0);
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

  // --- Numerical stability check ---
  function hasNumericalIssue(sol: EpidemicSolution): boolean {
    for (let i = 0; i < sol.t.length; i++) {
      if (!isFinite(sol.S[i]) || !isFinite(sol.I[i]) || !isFinite(sol.R[i]))
        return true;
      if (sol.E !== null && !isFinite(sol.E[i]))
        return true;
    }
    return false;
  }

  // --- DataFrame creation ---
  function createDataFrame(result: EpidemicSolution, modelType: 'SIR' | 'SEIR'): DG.DataFrame {
    const cols = [
      DG.Column.fromFloat64Array('Time', result.t),
      DG.Column.fromFloat64Array('S', result.S),
    ];

    if (modelType === 'SEIR' && result.E !== null)
      cols.push(DG.Column.fromFloat64Array('E', result.E));

    cols.push(DG.Column.fromFloat64Array('I', result.I));
    cols.push(DG.Column.fromFloat64Array('R', result.R));
    cols.push(DG.Column.fromFloat64Array('R_eff', result.rEff));

    const newDf = DG.DataFrame.fromColumns(cols);
    newDf.name = 'SIR-SEIR Epidemic';
    return newDf;
  }

  // --- Update DataFrame ---
  function updateDataFrame(result: EpidemicSolution, modelType: 'SIR' | 'SEIR'): void {
    const newDf = createDataFrame(result, modelType);
    view.dataFrame = newDf;
    lineChart.dataFrame = newDf;
    rEffChart.dataFrame = newDf;
    phaseChart.dataFrame = newDf;
  }

  function clearResults(): void {
    const cols = [
      DG.Column.fromFloat64Array('Time', new Float64Array(0)),
      DG.Column.fromFloat64Array('S', new Float64Array(0)),
      DG.Column.fromFloat64Array('I', new Float64Array(0)),
      DG.Column.fromFloat64Array('R', new Float64Array(0)),
      DG.Column.fromFloat64Array('R_eff', new Float64Array(0)),
    ];
    const emptyDf = DG.DataFrame.fromColumns(cols);
    emptyDf.name = 'SIR-SEIR Epidemic';
    view.dataFrame = emptyDf;
    lineChart.dataFrame = emptyDf;
    rEffChart.dataFrame = emptyDf;
    phaseChart.dataFrame = emptyDf;
  }

  // --- Preset buttons ---
  const btnInfluenza = ui.button('Influenza (R\u2080 \u2248 1.5)', () => applyPreset({
    modelType: ctrlModelType.value === true ? 'SEIR' : 'SIR',
    r0: 1.5, gamma: 0.143, sigma: 0.5,
  }), 'Influenza: R\u2080 \u2248 1.5, infectious period ~7 days, incubation ~2 days');

  const btnCovid = ui.button('COVID-19 (R\u2080 \u2248 3.0)', () => applyPreset({
    modelType: 'SEIR',
    r0: 3.0, gamma: 0.1, sigma: 0.2,
    switchToSEIR: true,
  }), 'COVID-19: R\u2080 \u2248 3.0, infectious period ~10 days, incubation ~5 days. Switches to SEIR model.');

  const btnMeasles = ui.button('Measles (R\u2080 \u2248 15)', () => applyPreset({
    modelType: 'SEIR',
    r0: 8.0, gamma: 0.125, sigma: 0.1,
    switchToSEIR: true,
  }), 'Measles: R\u2080 \u2248 15 (clamped to slider max 8.0), infectious period ~8 days, incubation ~10 days. Switches to SEIR model.');

  interface PresetConfig {
    modelType: 'SIR' | 'SEIR';
    r0: number;
    gamma: number;
    sigma: number;
    switchToSEIR?: boolean;
  }

  /** Applies a disease preset via batch update */
  function applyPreset(preset: PresetConfig): void {
    computationsBlocked = true;
    if (preset.switchToSEIR)
      ctrlModelType.value = true;
    ctrlR0.value = preset.r0;
    ctrlGamma.value = preset.gamma;
    ctrlSigma.value = preset.sigma;
    computationsBlocked = false;
    updateSigmaVisibility();
    runPrimary();
  }

  // --- Reset button ---
  const resetBtn = ui.iconFA('undo', () => {
    computationsBlocked = true;
    ctrlModelType.value = false;
    ctrlR0.value = DEFAULTS.r0;
    ctrlGamma.value = DEFAULTS.gamma;
    ctrlSigma.value = DEFAULTS.sigma;
    ctrlVaccination.value = DEFAULTS.vaccination;
    computationsBlocked = false;
    updateSigmaVisibility();
    runPrimary();
  }, 'Reset all parameters to default values');

  view.setRibbonPanels([[resetBtn]]);

  // --- Layout: left panel with form ---
  const form = ui.form([]);

  form.append(ui.h2('Model'));
  form.append(ctrlModelType.root);

  form.append(ui.h2('Parameters'));
  form.append(ctrlR0.root);
  form.append(ctrlGamma.root);
  form.append(ctrlSigma.root);
  form.append(ctrlVaccination.root);

  form.append(ui.h2('Disease Presets'));
  const presetGroup = ui.div([btnInfluenza, btnCovid, btnMeasles], 'sir-seir-app-preset-group');
  form.append(presetGroup);

  form.append(ui.h2('Epidemic Summary'));
  form.append(summaryPanel);

  const dockMng = view.dockManager;
  dockMng.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.6);

  // --- Epidemic dynamics line chart ---
  const yColumns = getInputs().modelType === 'SEIR' ? ['S', 'E', 'I', 'R'] : ['S', 'I', 'R'];
  lineChart = view.addViewer('Line chart', {
    xColumnName: 'Time',
    yColumnNames: yColumns,
    title: 'Epidemic Dynamics \u2014 Population Over Time',
  });

  // --- R_eff line chart ---
  rEffChart = view.addViewer('Line chart', {
    xColumnName: 'Time',
    yColumnNames: ['R_eff'],
    title: 'Effective Reproduction Number R_eff(t)',
  });

  // --- Phase portrait scatter plot ---
  phaseChart = view.addViewer('Scatter plot', {
    xColumnName: 'S',
    yColumnName: 'I',
    title: 'Phase Portrait \u2014 Susceptible vs Infectious',
    markerDefaultSize: 2,
  });

  // Dock viewers: line chart top-right, R_eff middle-right, phase portrait bottom-right
  dockMng.dock(lineChart, DG.DOCK_TYPE.RIGHT, null, undefined, 0.7);
  const lineChartNode = dockMng.findNode(lineChart.root);
  if (lineChartNode != null) {
    dockMng.dock(rEffChart, DG.DOCK_TYPE.DOWN, lineChartNode, undefined, 0.35);
    const rEffNode = dockMng.findNode(rEffChart.root);
    if (rEffNode != null)
      dockMng.dock(phaseChart, DG.DOCK_TYPE.DOWN, rEffNode, undefined, 0.5);
  }

  // --- Cleanup on close ---
  subs.push(grok.events.onViewRemoved.subscribe((v: any) => {
    if (v === view) {
      for (const sub of subs)
        sub.unsubscribe();
      if (debounceTimer !== null)
        clearTimeout(debounceTimer);
    }
  }));
}
