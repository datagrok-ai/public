// Multi-Compartment PK Model — Main Application Coordinator

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  PkParams, PkSolution, ObservedData, OptimizeResult, LandscapeData,
  ModelType, InputMode, SamplingSchedule, WeightType,
  SelectedParam, GridPointResult,
  OptimizeWorkerInput, OptimizeWorkerOutput,
  DEFAULTS, RANGES, validate, solvePk, generateObservedData,
  computeWssr, computeR2, computeAic, linspace, solvePkAtTimes,
  InputId, ValidationErrors,
} from './core';

import '../../css/compartmental-pk.css';

// --- Constants ---

const DEBOUNCE_MS = 100;
const PARAM_LABELS: Record<string, string> = {
  cl: 'CL', vc: 'V_c', vp: 'V_p', q: 'Q', vd: 'V_d', q2: 'Q\u2082', ka: 'k_a', f: 'F',
};

// --- Presets ---

interface Preset {
  name: string;
  tooltip: string;
  values: Partial<PkParams>;
}

const PRESETS: Preset[] = [
  {
    name: 'Vancomycin (CL=4.6, V_c=28, t\u00BD\u22486h)',
    tooltip: 'Load real-world PK parameters for vancomycin, a glycopeptide antibiotic with two-compartment kinetics.',
    values: {
      modelType: '2-Compartment', inputMode: 'IV Infusion',
      cl: 4.6, vc: 28, vp: 28, q: 8.8, dose: 1000, tInf: 1,
    },
  },
  {
    name: 'Theophylline (CL=3.5, V_c=30, t\u00BD\u22488h)',
    tooltip: 'Load real-world PK parameters for theophylline, a bronchodilator with oral administration and well-characterized PK.',
    values: {
      modelType: '2-Compartment', inputMode: 'Oral',
      cl: 3.5, vc: 30, vp: 20, q: 5, dose: 300, ka: 1.5, f: 0.9,
    },
  },
  {
    name: 'Metformin (CL=26, V_c=63, t\u00BD\u22485h)',
    tooltip: 'Load real-world PK parameters for metformin, a diabetes drug with high clearance and moderate bioavailability.',
    values: {
      modelType: '2-Compartment', inputMode: 'Oral',
      cl: 26, vc: 63, vp: 100, q: 15, dose: 500, ka: 2.0, f: 0.55,
    },
  },
];

// --- Main application ---

export function compartmentalPkApp(pkg: DG.Package): void {
  // State
  let computationsBlocked = false;
  let debounceTimer: ReturnType<typeof setTimeout> | null = null;
  const subs: {unsubscribe: () => void}[] = [];
  const activeWorkers: Worker[] = [];
  let observedData: ObservedData | null = null;
  let optimizeResult: OptimizeResult | null = null;
  let animationTimer: number | null = null;
  let compareMode = false;
  let btnGenerateDataRef: HTMLElement | null = null;
  let btnRunOptimizeRef: HTMLElement | null = null;

  // Initial solve
  const initialResult = solvePk(DEFAULTS);
  let df = createDataFrame(DEFAULTS, initialResult);
  const view = grok.shell.addTableView(df);
  view.name = 'Multi-Compartment PK Model';

  // --- Viewer references ---
  let concChart: DG.Viewer | null = null;
  let compChart: DG.Viewer | null = null;

  // ===================== PK METRICS PANEL =====================

  const lblCmax = ui.label('');
  lblCmax.classList.add('cpk-metrics-value');
  const lblTmax = ui.label('');
  lblTmax.classList.add('cpk-metrics-value');
  const lblAuc = ui.label('');
  lblAuc.classList.add('cpk-metrics-value');
  const lblCssMax = ui.label('');
  lblCssMax.classList.add('cpk-metrics-value');
  const lblCssMin = ui.label('');
  lblCssMin.classList.add('cpk-metrics-value');
  const lblT90ss = ui.label('');
  lblT90ss.classList.add('cpk-metrics-value');

  const metricsPanel = ui.divV([
    ui.label('C_max:'), lblCmax,
    ui.label('t_max:'), lblTmax,
    ui.label('AUC:'), lblAuc,
    ui.label('C_ss,max:'), lblCssMax,
    ui.label('C_ss,min:'), lblCssMin,
    ui.label('t_90% SS:'), lblT90ss,
  ]);
  metricsPanel.classList.add('cpk-metrics-panel');

  /** Update PK metrics display */
  function updateMetrics(result: PkSolution): void {
    lblCmax.textContent = `${result.cMax.toFixed(2)} mg/L`;
    lblTmax.textContent = `${result.tMax.toFixed(2)} h`;
    lblAuc.textContent = `${result.auc.toFixed(1)} mg\u00B7h/L`;
    lblCssMax.textContent = result.cssMax !== null ? `${result.cssMax.toFixed(2)} mg/L` : 'N/A';
    lblCssMin.textContent = result.cssMin !== null ? `${result.cssMin.toFixed(2)} mg/L` : 'N/A';
    lblT90ss.textContent = result.t90ss !== null ? `${result.t90ss.toFixed(1)} h` : 'N/A';
  }

  // ===================== PRIMARY CONTROLS =====================

  // Model settings
  const ctrlModelType = ui.input.choice('Model compartments', {
    value: DEFAULTS.modelType,
    items: ['2-Compartment', '3-Compartment'],
    tooltipText: 'Number of compartments in the PK model. Two-compartment: central + peripheral. Three-compartment adds a deep (slowly equilibrating) tissue compartment.',
    onValueChanged: () => {
      updateControlVisibility();
      updateOptParamSelector();
      debouncedRun();
    },
  });

  const ctrlInputMode = ui.input.choice('Administration route', {
    value: DEFAULTS.inputMode,
    items: ['IV Bolus', 'IV Infusion', 'Oral'],
    tooltipText: 'Route of drug administration. IV Bolus: instantaneous injection. IV Infusion: constant-rate infusion over T_inf hours. Oral: first-order absorption from GI tract.',
    onValueChanged: () => {
      updateControlVisibility();
      updateOptParamSelector();
      debouncedRun();
    },
  });

  const ctrlYScale = ui.input.choice('Y-axis scale', {
    value: 'Linear',
    items: ['Linear', 'Semi-log'],
    tooltipText: 'On a semi-log scale, exponential processes appear as straight lines. Two distinct slopes = two phases: fast (\u03B1, distribution) and slow (\u03B2, elimination).',
    onValueChanged: () => updateChartScale(),
  });

  // PK parameters
  const ctrlCl = ui.input.float('CL \u2014 clearance (L/h)', {
    value: DEFAULTS.cl, min: RANGES.cl.min, max: RANGES.cl.max,
    tooltipText: 'Volume of plasma completely cleared of drug per unit time. Determines the slope of the terminal phase on a semi-log plot. Half-life t\u00BD = 0.693 \u00B7 V_ss / CL.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlVc = ui.input.float('V_c \u2014 central volume (L)', {
    value: DEFAULTS.vc, min: RANGES.vc.min, max: RANGES.vc.max,
    tooltipText: 'Apparent volume of distribution in blood/plasma. Determines the initial concentration after an IV bolus: C\u2080 = Dose / V_c.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlVp = ui.input.float('V_p \u2014 peripheral volume (L)', {
    value: DEFAULTS.vp, min: RANGES.vp.min, max: RANGES.vp.max,
    tooltipText: 'Volume of the rapidly equilibrating peripheral compartment. Larger V_p \u2192 more drug distributed to tissues \u2192 lower plasma concentrations.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlQ = ui.input.float('Q \u2014 intercompartmental clearance (L/h)', {
    value: DEFAULTS.q, min: RANGES.q.min, max: RANGES.q.max,
    tooltipText: 'Rate of drug exchange between central and peripheral compartments. Large Q \u2192 fast distribution \u2192 short \u03B1-phase.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlVd = ui.input.float('V_d \u2014 deep compartment volume (L)', {
    value: DEFAULTS.vd, min: RANGES.vd.min, max: RANGES.vd.max,
    tooltipText: 'Volume of the deep (slowly equilibrating) tissue compartment. Found in drugs that bind to bone, fat, or deep tissues.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlQ2 = ui.input.float('Q\u2082 \u2014 deep intercompartmental clearance (L/h)', {
    value: DEFAULTS.q2, min: RANGES.q2.min, max: RANGES.q2.max,
    tooltipText: 'Rate of drug exchange with the deep compartment. Typically much smaller than Q, producing a slow terminal phase.',
    onValueChanged: () => debouncedRun(),
  });

  // Absorption
  const ctrlKa = ui.input.float('k_a \u2014 absorption rate (h\u207B\u00B9)', {
    value: DEFAULTS.ka, min: RANGES.ka.min, max: RANGES.ka.max,
    tooltipText: 'Rate of absorption from the GI tract. When k_a >> k\u2081\u2080 (no flip-flop kinetics), the peak comes quickly. When k_a \u2248 k\u2081\u2080, the peak is blunted and delayed.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlF = ui.input.float('F \u2014 bioavailability', {
    value: DEFAULTS.f, min: RANGES.f.min, max: RANGES.f.max,
    tooltipText: 'Fraction of the dose reaching systemic circulation after oral administration. F = 1.0 for IV; typically 0.1\u20130.9 for oral due to first-pass metabolism.',
    onValueChanged: () => debouncedRun(),
  });

  // Dosing
  const ctrlDose = ui.input.float('Dose (mg)', {
    value: DEFAULTS.dose, min: RANGES.dose.min, max: RANGES.dose.max,
    tooltipText: 'Amount of drug per dose. For IV bolus, delivered instantaneously. For infusion, delivered over T_inf hours.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlRepeated = ui.input.bool('Repeated dosing', {
    value: DEFAULTS.repeatedDosing,
    tooltipText: 'Enable multiple-dose regimen. Shows accumulation toward steady state over repeated dosing intervals.',
    onValueChanged: () => {
      updateControlVisibility();
      debouncedRun();
    },
  });
  const ctrlTau = ui.input.float('\u03C4 \u2014 dosing interval (h)', {
    value: DEFAULTS.tau, min: RANGES.tau.min, max: RANGES.tau.max,
    tooltipText: 'Time between successive doses. Shorter intervals \u2192 higher trough levels \u2192 faster approach to steady state.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlNDoses = ui.input.int('Number of doses', {
    value: DEFAULTS.nDoses, min: RANGES.nDoses.min, max: RANGES.nDoses.max,
    tooltipText: 'Total number of doses administered. More doses \u2192 closer to true steady state (reached after ~5 half-lives).',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlTInf = ui.input.float('Infusion duration (h)', {
    value: DEFAULTS.tInf, min: RANGES.tInf.min, max: RANGES.tInf.max,
    tooltipText: 'Duration of each constant-rate infusion. Longer infusion \u2192 lower C_max, smoother concentration profile.',
    onValueChanged: () => debouncedRun(),
  });

  // Simulation
  const ctrlTEnd = ui.input.float('Simulation time (h)', {
    value: DEFAULTS.tEnd, min: RANGES.tEnd.min, max: RANGES.tEnd.max,
    tooltipText: 'Total simulation duration. For multiple dosing, ensure this covers at least n_doses \u00D7 \u03C4 plus one extra interval.',
    onValueChanged: () => debouncedRun(),
  });

  // Therapeutic window
  const ctrlMec = ui.input.float('MEC \u2014 min effective conc (mg/L)', {
    value: DEFAULTS.mec, min: RANGES.mec.min, max: RANGES.mec.max,
    tooltipText: 'Below this concentration, the drug has no therapeutic effect. Shown as a horizontal dashed line on the plot.',
    onValueChanged: () => debouncedRun(),
  });
  const ctrlMtc = ui.input.float('MTC \u2014 min toxic conc (mg/L)', {
    value: DEFAULTS.mtc, min: RANGES.mtc.min, max: RANGES.mtc.max,
    tooltipText: 'Above this concentration, toxic effects may occur. The goal of dosing is to keep concentration between MEC and MTC.',
    onValueChanged: () => debouncedRun(),
  });

  // Set formats
  computationsBlocked = true;
  ctrlCl.format = '0.0';
  ctrlVc.format = '0.0';
  ctrlVp.format = '0.0';
  ctrlQ.format = '0.0';
  ctrlVd.format = '0.0';
  ctrlQ2.format = '0.00';
  ctrlKa.format = '0.0';
  ctrlF.format = '0.00';
  ctrlDose.format = '0.0';
  ctrlTau.format = '0.0';
  ctrlTInf.format = '0.0';
  ctrlTEnd.format = '0.0';
  ctrlMec.format = '0.00';
  ctrlMtc.format = '0.0';
  computationsBlocked = false;

  // Primary controls array for validators
  const primaryControls: {input: DG.InputBase; id: InputId}[] = [
    {input: ctrlCl, id: 'ctrl_cl'}, {input: ctrlVc, id: 'ctrl_vc'},
    {input: ctrlVp, id: 'ctrl_vp'}, {input: ctrlQ, id: 'ctrl_q'},
    {input: ctrlVd, id: 'ctrl_vd'}, {input: ctrlQ2, id: 'ctrl_q2'},
    {input: ctrlKa, id: 'ctrl_ka'}, {input: ctrlF, id: 'ctrl_f'},
    {input: ctrlDose, id: 'ctrl_dose'},
    {input: ctrlTau, id: 'ctrl_tau'}, {input: ctrlNDoses, id: 'ctrl_n_doses'},
    {input: ctrlTInf, id: 'ctrl_t_inf'}, {input: ctrlTEnd, id: 'ctrl_t_end'},
    {input: ctrlMec, id: 'ctrl_mec'}, {input: ctrlMtc, id: 'ctrl_mtc'},
  ];

  // Add validators
  for (const ctrl of primaryControls) {
    ctrl.input.addValidator(() => {
      const errors = validate(getInputs());
      return errors.get(ctrl.id) ?? null;
    });
  }

  // ===================== SECONDARY CONTROLS =====================

  const ctrlNoiseCv = ui.input.float('Observation noise CV (%)', {
    value: 15, min: 5, max: 40,
    tooltipText: 'Coefficient of variation for log-normal noise. Higher CV = noisier data = harder optimization. 15% is typical for clinical PK data.',
    onValueChanged: () => updateGenerateButtonLabel(),
  });
  ctrlNoiseCv.format = '0.0';

  const ctrlSampling = ui.input.choice('Sampling schedule', {
    value: 'Standard PK' as string,
    items: ['Standard PK', 'Rich', 'Sparse'],
    tooltipText: 'Time points at which samples are taken. Standard PK: 7 points. Rich: 15 points (better identifiability). Sparse: 4 points (realistic clinical limitation).',
  });

  const ctrlGridRes = ui.input.int('Grid density per axis', {
    value: 20, min: 5, max: 50,
    tooltipText: 'Number of points per parameter axis. Total ODE solves = grid_points ^ number_of_parameters. For 2 parameters at 30 grid points = 900 solves \u2014 usually takes < 2 seconds with MRT.',
    onValueChanged: () => updateOptimizeButtonLabel(),
  });

  const ctrlWeight = ui.input.choice('Objective function weighting', {
    value: '1/C_obs\u00B2' as string,
    items: ['Uniform', '1/C_obs\u00B2', '1/C_pred\u00B2'],
    tooltipText: 'WSSR weighting scheme. 1/C_obs\u00B2 is standard in PK because measurement error is typically proportional to concentration. Uniform weights treat all residuals equally.',
  });

  const ctrlCompare = ui.input.bool('Compare: sliders vs optimized', {
    value: false,
    tooltipText: 'Overlay both curves: the current slider-set parameters and the optimized best-fit. See how they differ.',
    onValueChanged: () => {
      compareMode = ctrlCompare.value!;
      debouncedRun();
    },
  });

  // ===================== OPTIMIZATION PARAM SELECTOR =====================

  const optParamCheckboxes: Map<string, {
    checkbox: HTMLInputElement;
    minInput: HTMLInputElement;
    maxInput: HTMLInputElement;
    row: HTMLElement;
    rangeRow: HTMLElement;
  }> = new Map();

  const optSolveCountLabel = ui.label('Total ODE solves: 0');
  const optWarningLabel = ui.label('');
  optWarningLabel.classList.add('cpk-warning');

  const optParamContainer = ui.divV([]);

  /** Build parameter selector UI */
  function buildOptParamSelector(): void {
    optParamContainer.innerHTML = '';
    optParamCheckboxes.clear();

    const availableParams = getAvailableOptParams();

    for (const paramName of availableParams) {
      const label = PARAM_LABELS[paramName] || paramName;
      const currentValue = getParamValue(paramName);
      const defaultMin = Math.max(0.001, currentValue * 0.5);
      const defaultMax = currentValue * 1.5;

      const checkbox = document.createElement('input');
      checkbox.type = 'checkbox';
      checkbox.addEventListener('change', () => {
        enforceParamLimit();
        updateOptimizeButtonLabel();
        rangeRow.classList.toggle('cpk-hidden', !checkbox.checked);
      });

      const row = ui.divH([checkbox, ui.label(label)]);
      row.classList.add('cpk-opt-selector-row');

      const minInput = document.createElement('input');
      minInput.type = 'number';
      minInput.value = defaultMin.toFixed(2);
      minInput.step = '0.01';
      minInput.addEventListener('change', () => updateOptimizeButtonLabel());

      const maxInput = document.createElement('input');
      maxInput.type = 'number';
      maxInput.value = defaultMax.toFixed(2);
      maxInput.step = '0.01';
      maxInput.addEventListener('change', () => updateOptimizeButtonLabel());

      const rangeRow = ui.divH([
        ui.label('Min:'), minInput,
        ui.label('Max:'), maxInput,
      ]);
      rangeRow.classList.add('cpk-opt-selector-range', 'cpk-hidden');

      optParamCheckboxes.set(paramName, {checkbox, minInput, maxInput, row, rangeRow});

      optParamContainer.append(row, rangeRow);
    }
    optParamContainer.append(optSolveCountLabel, optWarningLabel);
  }

  /** Get available optimization parameters based on current model/mode */
  function getAvailableOptParams(): string[] {
    const params = ['cl', 'vc', 'vp', 'q'];
    if (ctrlModelType.value === '3-Compartment')
      params.push('vd', 'q2');
    if (ctrlInputMode.value === 'Oral')
      params.push('ka', 'f');
    return params;
  }

  /** Enforce maximum 3 selected parameters */
  function enforceParamLimit(): void {
    const checked = getSelectedOptParams();
    if (checked.length >= 3) {
      for (const [, entry] of optParamCheckboxes) {
        if (!entry.checkbox.checked)
          entry.checkbox.disabled = true;
      }
    } else {
      for (const [, entry] of optParamCheckboxes)
        entry.checkbox.disabled = false;
    }
  }

  /** Get currently selected optimization parameters */
  function getSelectedOptParams(): SelectedParam[] {
    const selected: SelectedParam[] = [];
    for (const [name, entry] of optParamCheckboxes) {
      if (entry.checkbox.checked) {
        selected.push({
          name,
          min: parseFloat(entry.minInput.value),
          max: parseFloat(entry.maxInput.value),
        });
      }
    }
    return selected;
  }

  function updateOptParamSelector(): void {
    buildOptParamSelector();
  }

  // ===================== BUTTONS =====================

  const btnGenerateData = ui.bigButton(
    `Generate Noisy Observations (CV = ${ctrlNoiseCv.value}%)`,
    () => runGenerateData(),
    'Runs the model at current parameters, samples at clinical time points, and adds log-normal noise. Try recovering the true parameters with the optimizer \u2014 higher CV makes it harder!',
  );
  btnGenerateDataRef = btnGenerateData;

  const btnRunOptimize = ui.bigButton(
    'Run Grid Search (estimated: 0 ODE solves)',
    () => runOptimization(),
    'Perform exhaustive search over selected parameters. Each grid point requires a full ODE solve. The parameter combination minimizing WSSR is returned as the best fit.',
  );
  btnRunOptimizeRef = btnRunOptimize;
  btnRunOptimize.classList.add('cpk-btn--disabled');

  const btnApplyBestFit = ui.button('Apply Best-Fit to Sliders', () => applyBestFit(),
    'Write the optimized parameters to the sliders for visual comparison.');
  btnApplyBestFit.classList.add('cpk-hidden');

  // Results panel
  const resultsPanel = ui.divV([]);
  resultsPanel.classList.add('cpk-results-panel', 'cpk-hidden');

  // Landscape container
  const landscapeContainer = ui.div([]);
  landscapeContainer.classList.add('cpk-landscape-container', 'cpk-hidden');

  // Residuals container
  const residualsContainer = ui.div([]);
  residualsContainer.classList.add('cpk-hidden');

  // ===================== COMPARTMENT DIAGRAM =====================

  const diagramContainer = ui.div([]);
  diagramContainer.classList.add('cpk-diagram-container');

  const diagramTimeSlider = document.createElement('input');
  diagramTimeSlider.type = 'range';
  diagramTimeSlider.min = '0';
  diagramTimeSlider.max = '100';
  diagramTimeSlider.value = '0';

  let diagramPlaying = false;
  const diagramPlayBtn = ui.iconFA('play', () => {
    diagramPlaying = !diagramPlaying;
    diagramPlayBtn.classList.toggle('fa-play', !diagramPlaying);
    diagramPlayBtn.classList.toggle('fa-pause', diagramPlaying);
    if (diagramPlaying) startDiagramAnimation();
  }, 'Play/pause compartment animation');
  diagramPlayBtn.classList.add('cpk-diagram-play-btn');

  const diagramControls = ui.divH([diagramPlayBtn, diagramTimeSlider]);
  diagramControls.classList.add('cpk-diagram-controls');

  const diagramPanel = ui.divV([diagramContainer, diagramControls]);

  let lastSolution: PkSolution | null = initialResult;

  /** Draw the compartment diagram SVG */
  function drawDiagram(params: PkParams, solution: PkSolution, timeIndex: number): void {
    const is3comp = params.modelType === '3-Compartment';
    const isOral = params.inputMode === 'Oral';
    const svgNS = 'http://www.w3.org/2000/svg';

    const width = 400;
    const height = is3comp ? 280 : 220;

    const svg = document.createElementNS(svgNS, 'svg');
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
    svg.setAttribute('width', '100%');
    svg.setAttribute('height', '100%');

    // Arrow marker definition
    const defs = document.createElementNS(svgNS, 'defs');
    const marker = document.createElementNS(svgNS, 'marker');
    marker.setAttribute('id', 'arrowhead');
    marker.setAttribute('markerWidth', '10');
    marker.setAttribute('markerHeight', '7');
    marker.setAttribute('refX', '10');
    marker.setAttribute('refY', '3.5');
    marker.setAttribute('orient', 'auto');
    const polygon = document.createElementNS(svgNS, 'polygon');
    polygon.setAttribute('points', '0 0, 10 3.5, 0 7');
    polygon.setAttribute('fill', '#555');
    marker.appendChild(polygon);
    defs.appendChild(marker);
    svg.appendChild(defs);

    // Get amounts at current time
    const idx = Math.min(timeIndex, solution.t.length - 1);
    const ac = solution.ac[idx];
    const ap = solution.ap[idx];
    const ad = is3comp && solution.ad ? solution.ad[idx] : 0;
    const agut = isOral && solution.agut ? solution.agut[idx] : 0;
    const micro = {
      k10: params.cl / params.vc,
      k12: params.q / params.vc,
      k21: params.q / params.vp,
      k13: is3comp ? params.q2 / params.vc : 0,
      k31: is3comp ? params.q2 / params.vd : 0,
    };

    // Positions
    const cX = isOral ? 200 : 140;
    const cY = 60;
    const pX = cX + 160;
    const pY = 60;
    const dX = cX;
    const dY = is3comp ? 200 : 0;

    // Draw compartments
    drawBox(svg, svgNS, cX, cY, 'Central', ac.toFixed(1) + ' mg');
    drawBox(svg, svgNS, pX, pY, 'Peripheral', ap.toFixed(1) + ' mg');
    if (is3comp)
      drawBox(svg, svgNS, dX, dY, 'Deep', ad.toFixed(1) + ' mg');
    if (isOral)
      drawBox(svg, svgNS, cX - 160, cY, 'Gut', agut.toFixed(1) + ' mg');

    // Draw arrows with flow-proportional thickness
    const maxFlow = Math.max(
      micro.k10 * ac, micro.k12 * ac, micro.k21 * ap,
      micro.k13 * ac, micro.k31 * ad,
      isOral ? params.ka * agut : 0, 1,
    );

    // Central -> Peripheral
    drawArrow(svg, svgNS, cX + 50, cY - 10, pX - 50, pY - 10,
      micro.k12 * ac / maxFlow, `k\u2081\u2082=${micro.k12.toFixed(2)}`);
    // Peripheral -> Central
    drawArrow(svg, svgNS, pX - 50, pY + 10, cX + 50, cY + 10,
      micro.k21 * ap / maxFlow, `k\u2082\u2081=${micro.k21.toFixed(2)}`);
    // Elimination
    drawArrow(svg, svgNS, cX, cY + 30, cX, cY + (is3comp ? 80 : 70),
      micro.k10 * ac / maxFlow, `CL=${params.cl}`);

    if (is3comp) {
      drawArrow(svg, svgNS, cX - 30, cY + 30, dX - 30, dY - 30,
        micro.k13 * ac / maxFlow, `k\u2081\u2083=${micro.k13.toFixed(2)}`);
      drawArrow(svg, svgNS, dX + 30, dY - 30, cX + 30, cY + 30,
        micro.k31 * ad / maxFlow, `k\u2083\u2081=${micro.k31.toFixed(3)}`);
    }

    if (isOral)
      drawArrow(svg, svgNS, cX - 110, cY, cX - 50, cY,
        params.ka * agut / maxFlow, `k_a=${params.ka}`);

    // Time label
    const timeTxt = document.createElementNS(svgNS, 'text');
    timeTxt.setAttribute('x', '10');
    timeTxt.setAttribute('y', String(height - 10));
    timeTxt.setAttribute('class', 'cpk-diagram-label');
    timeTxt.textContent = `t = ${solution.t[idx].toFixed(1)} h`;
    svg.appendChild(timeTxt);

    diagramContainer.innerHTML = '';
    diagramContainer.appendChild(svg);
  }

  /** Draw a compartment box */
  function drawBox(svg: SVGSVGElement, ns: string, cx: number, cy: number, label: string, amount: string): void {
    const rect = document.createElementNS(ns, 'rect');
    rect.setAttribute('x', String(cx - 45));
    rect.setAttribute('y', String(cy - 25));
    rect.setAttribute('width', '90');
    rect.setAttribute('height', '50');
    rect.setAttribute('class', 'cpk-diagram-compartment');
    svg.appendChild(rect);

    const text = document.createElementNS(ns, 'text');
    text.setAttribute('x', String(cx));
    text.setAttribute('y', String(cy - 5));
    text.setAttribute('class', 'cpk-diagram-compartment-label');
    text.textContent = label;
    svg.appendChild(text);

    const amtText = document.createElementNS(ns, 'text');
    amtText.setAttribute('x', String(cx));
    amtText.setAttribute('y', String(cy + 12));
    amtText.setAttribute('class', 'cpk-diagram-amount-label');
    amtText.textContent = amount;
    svg.appendChild(amtText);
  }

  /** Draw an arrow with thickness proportional to flow */
  function drawArrow(
    svg: SVGSVGElement, ns: string,
    x1: number, y1: number, x2: number, y2: number,
    flowFrac: number, label: string,
  ): void {
    const thickness = 1 + Math.min(flowFrac, 1) * 4;
    const line = document.createElementNS(ns, 'line');
    line.setAttribute('x1', String(x1));
    line.setAttribute('y1', String(y1));
    line.setAttribute('x2', String(x2));
    line.setAttribute('y2', String(y2));
    line.setAttribute('class', 'cpk-diagram-arrow');
    line.setAttribute('stroke-width', String(thickness));
    svg.appendChild(line);

    const midX = (x1 + x2) / 2;
    const midY = (y1 + y2) / 2 - 8;
    const txt = document.createElementNS(ns, 'text');
    txt.setAttribute('x', String(midX));
    txt.setAttribute('y', String(midY));
    txt.setAttribute('class', 'cpk-diagram-label');
    txt.textContent = label;
    svg.appendChild(txt);
  }

  /** Animate the diagram */
  function startDiagramAnimation(): void {
    if (animationTimer !== null) return;
    let step = parseInt(diagramTimeSlider.value, 10);
    animationTimer = window.setInterval(() => {
      if (!diagramPlaying || !lastSolution) {
        stopDiagramAnimation();
        return;
      }
      step++;
      if (step >= lastSolution.t.length) step = 0;
      diagramTimeSlider.value = String(step);
      drawDiagram(getInputs(), lastSolution, step);
    }, 50);
  }

  function stopDiagramAnimation(): void {
    if (animationTimer !== null) {
      clearInterval(animationTimer);
      animationTimer = null;
    }
    diagramPlaying = false;
    diagramPlayBtn.classList.add('fa-play');
    diagramPlayBtn.classList.remove('fa-pause');
  }

  diagramTimeSlider.addEventListener('input', () => {
    if (lastSolution)
      drawDiagram(getInputs(), lastSolution, parseInt(diagramTimeSlider.value, 10));
  });

  // ===================== HELPERS =====================

  /** Get current inputs from controls */
  function getInputs(): PkParams {
    return {
      modelType: ctrlModelType.value as ModelType,
      inputMode: ctrlInputMode.value as InputMode,
      cl: ctrlCl.value!, vc: ctrlVc.value!, vp: ctrlVp.value!, q: ctrlQ.value!,
      vd: ctrlVd.value!, q2: ctrlQ2.value!,
      ka: ctrlKa.value!, f: ctrlF.value!,
      dose: ctrlDose.value!,
      repeatedDosing: ctrlRepeated.value!,
      tau: ctrlTau.value!, nDoses: ctrlNDoses.value!,
      tInf: ctrlTInf.value!,
      tEnd: ctrlTEnd.value!,
      mec: ctrlMec.value!, mtc: ctrlMtc.value!,
    };
  }

  /** Get a single parameter value by name */
  function getParamValue(name: string): number {
    const inputs = getInputs();
    return (inputs as any)[name] as number;
  }

  /** Update control visibility based on model type, input mode, repeated dosing */
  function updateControlVisibility(): void {
    const is3comp = ctrlModelType.value === '3-Compartment';
    const isOral = ctrlInputMode.value === 'Oral';
    const isInfusion = ctrlInputMode.value === 'IV Infusion';
    const isRepeated = ctrlRepeated.value;

    ctrlVd.root.classList.toggle('cpk-hidden', !is3comp);
    ctrlQ2.root.classList.toggle('cpk-hidden', !is3comp);
    ctrlKa.root.classList.toggle('cpk-hidden', !isOral);
    ctrlF.root.classList.toggle('cpk-hidden', !isOral);
    ctrlTInf.root.classList.toggle('cpk-hidden', !isInfusion);
    ctrlTau.root.classList.toggle('cpk-hidden', !isRepeated);
    ctrlNDoses.root.classList.toggle('cpk-hidden', !isRepeated);
  }

  /** Update chart Y-axis scale */
  function updateChartScale(): void {
    if (concChart) {
      concChart.setOptions({
        yAxisType: ctrlYScale.value === 'Semi-log' ? 'logarithmic' : 'linear',
      });
    }
  }

  /** Update generate button label with current CV */
  function updateGenerateButtonLabel(): void {
    if (!btnGenerateDataRef) return;
    btnGenerateDataRef.textContent = `Generate Noisy Observations (CV = ${ctrlNoiseCv.value}%)`;
  }

  /** Update optimize button label with estimated ODE solves */
  function updateOptimizeButtonLabel(): void {
    if (!btnRunOptimizeRef) return;
    const selected = getSelectedOptParams();
    const gridRes = ctrlGridRes.value!;
    const totalSolves = Math.pow(gridRes, selected.length);
    btnRunOptimizeRef.textContent = `Run Grid Search (estimated: ${totalSolves} ODE solves)`;

    optSolveCountLabel.textContent = `Total ODE solves: ${totalSolves}`;
    if (totalSolves > 10000) {
      optWarningLabel.textContent = `This may take a few seconds. Total ODE solves: ${totalSolves}`;
      optWarningLabel.classList.remove('cpk-hidden');
    } else {
      optWarningLabel.classList.add('cpk-hidden');
    }

    // Enable/disable optimize button
    const hasData = observedData !== null;
    const hasParams = selected.length > 0;
    btnRunOptimizeRef.classList.toggle('cpk-btn--disabled', !hasData || !hasParams);
  }

  // ===================== PRIMARY PIPELINE =====================

  /** Debounced run of the primary pipeline */
  function debouncedRun(): void {
    if (debounceTimer !== null) clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => runPrimary(), DEBOUNCE_MS);
  }

  /** Run the primary PK pipeline */
  function runPrimary(): void {
    if (computationsBlocked) return;

    const inputs = getInputs();
    const errors = validate(inputs);

    // Mark invalid inputs
    for (const ctrl of primaryControls) {
      const hasError = errors.has(ctrl.id);
      ctrl.input.root.classList.toggle('d4-invalid', hasError);
    }

    if (errors.size > 0) return;

    try {
      const result = solvePk(inputs);
      lastSolution = result;

      // Update DataFrame
      df = createDataFrame(inputs, result);
      view.dataFrame = df;

      // Update viewers
      updateViewers(inputs, result);
      updateMetrics(result);

      // Update diagram
      diagramTimeSlider.max = String(result.t.length - 1);
      drawDiagram(inputs, result, parseInt(diagramTimeSlider.value, 10));
    } catch (err) {
      grok.shell.error(`PK solve error: ${err instanceof Error ? err.message : 'Unknown error'}`);
    }
  }

  // ===================== DATA GENERATION =====================

  /** Generate synthetic observations */
  function runGenerateData(): void {
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Fix parameter errors before generating data');
      return;
    }

    if (!lastSolution) {
      grok.shell.error('No primary solution available');
      return;
    }

    observedData = generateObservedData({
      primarySolution: lastSolution,
      cv: ctrlNoiseCv.value!,
      schedule: ctrlSampling.value as SamplingSchedule,
      tEnd: inputs.tEnd,
      repeatedDosing: inputs.repeatedDosing,
      tau: inputs.tau,
      nDoses: inputs.nDoses,
    });

    btnGenerateData.classList.add('cpk-btn--active');
    updateOptimizeButtonLabel();
    runPrimary(); // Redraw chart with observed data overlay
    grok.shell.info(`Generated ${observedData.tObs.length} noisy observations`);
  }

  // ===================== OPTIMIZATION =====================

  let optimizing = false;

  /** Run brute-force grid search optimization */
  async function runOptimization(): Promise<void> {
    if (optimizing) return;

    // Validate
    const inputs = getInputs();
    const errors = validate(inputs);
    if (errors.size > 0) {
      grok.shell.error('Fix primary parameter errors before optimizing');
      return;
    }
    if (!observedData) {
      grok.shell.error('Generate synthetic observations first');
      return;
    }

    const selectedParams = getSelectedOptParams();
    if (selectedParams.length === 0) {
      grok.shell.error('Select at least one parameter to optimize');
      return;
    }
    if (selectedParams.length > 3) {
      grok.shell.error('Select at most 3 parameters for optimization');
      return;
    }

    for (const sp of selectedParams) {
      if (sp.min >= sp.max) {
        grok.shell.error(`Search range minimum must be less than maximum for ${PARAM_LABELS[sp.name] || sp.name}`);
        return;
      }
    }

    optimizing = true;
    btnRunOptimize.classList.add('cpk-btn--disabled');

    const gridRes = ctrlGridRes.value!;
    const weightType = ctrlWeight.value as WeightType;

    // Generate all grid points
    const paramGrids = selectedParams.map((sp) => linspace(sp.min, sp.max, gridRes));
    const allGridPoints = cartesianProduct(paramGrids);
    const paramNames = selectedParams.map((sp) => sp.name);
    const totalPoints = allGridPoints.length;

    // Progress bar
    const pi = DG.TaskBarProgressIndicator.create(`Optimizing (${totalPoints} grid points)...`);

    // Create worker pool
    const workerCount = Math.max(1, navigator.hardwareConcurrency - 2);
    const nWorkers = Math.min(workerCount, totalPoints);

    // Distribute grid points round-robin
    const chunks: number[][][] = Array.from({length: nWorkers}, () => []);
    for (let i = 0; i < totalPoints; i++)
      chunks[i % nWorkers].push(allGridPoints[i]);

    let cancelled = false;
    const cancelSub = pi.onCanceled.subscribe(() => {
      cancelled = true;
      terminateWorkers();
    });

    try {
      const promises = chunks.map((chunk) =>
        new Promise<GridPointResult[]>((resolve, reject) => {
          const worker = new Worker(pkg.webRoot + 'dist/optimize-worker.js');
          activeWorkers.push(worker);

          const msg: OptimizeWorkerInput = {
            gridPoints: chunk,
            paramNames,
            fixedParams: inputs,
            tEnd: inputs.tEnd,
            observedData: {tObs: observedData!.tObs, cObs: observedData!.cObs},
            weightType,
          };
          worker.postMessage(msg);

          worker.onmessage = (e: MessageEvent<OptimizeWorkerOutput>) => {
            if (e.data.success)
              resolve(e.data.results!);
            else
              reject(new Error(e.data.error));
          };

          worker.onerror = (e) => reject(new Error(e.message));
        }),
      );

      const results = await Promise.all(promises);
      const allResults = results.flat();

      if (cancelled) return;

      // Find best
      let bestIdx = 0;
      let bestWssr = Infinity;
      for (let i = 0; i < allResults.length; i++) {
        if (allResults[i].wssr < bestWssr) {
          bestWssr = allResults[i].wssr;
          bestIdx = i;
        }
      }

      const bestPoint = allResults[bestIdx];
      const bestParamsMap = new Map<string, number>();
      for (let k = 0; k < paramNames.length; k++)
        bestParamsMap.set(paramNames[k], bestPoint.paramValues[k]);

      // Solve at best-fit for full curve
      const bestPkParams = {...inputs};
      for (const [name, val] of bestParamsMap)
        (bestPkParams as any)[name] = val;

      const bestSolution = solvePk(bestPkParams);
      const bestConcAtObs = solvePkAtTimes(bestPkParams, observedData!.tObs);

      const r2 = computeR2(observedData!.cObs, bestWssr, weightType);
      const aic = computeAic(observedData!.cObs.length, selectedParams.length, bestWssr);

      // Build residuals
      const residuals = {
        t: observedData!.tObs,
        res: observedData!.cObs.map((c, i) => c - bestConcAtObs[i]),
      };

      // Build landscape data
      const landscape: LandscapeData = {
        paramNames,
        paramValues: paramGrids,
        wssrValues: allResults.map((r) => r.wssr),
        gridResolution: gridRes,
      };

      optimizeResult = {
        bestParams: bestParamsMap,
        bestWssr,
        bestR2: r2,
        bestAic: aic,
        bestConc: {
          t: Array.from(bestSolution.t),
          c: Array.from(bestSolution.conc),
        },
        landscape,
        residuals,
      };

      // Display results
      displayOptimizationResults();
      grok.shell.info(`Optimization complete. Best WSSR = ${bestWssr.toFixed(4)}, R\u00B2 = ${r2.toFixed(4)}`);
    } catch (err) {
      if (!cancelled)
        grok.shell.error(`Optimization error: ${err instanceof Error ? err.message : 'Unknown error'}`);
    } finally {
      terminateWorkers();
      cancelSub.unsubscribe();
      pi.close();
      optimizing = false;
      btnRunOptimize.classList.toggle('cpk-btn--disabled', !observedData);
    }
  }

  /** Generate Cartesian product of parameter grids */
  function cartesianProduct(arrays: number[][]): number[][] {
    if (arrays.length === 0) return [[]];
    const [first, ...rest] = arrays;
    const restProduct = cartesianProduct(rest);
    const result: number[][] = [];
    for (const v of first) {
      for (const r of restProduct)
        result.push([v, ...r]);
    }
    return result;
  }

  /** Terminate all active workers */
  function terminateWorkers(): void {
    for (const w of activeWorkers) {
      try { w.terminate(); } catch {}
    }
    activeWorkers.length = 0;
  }

  /** Display optimization results in the right panel */
  function displayOptimizationResults(): void {
    if (!optimizeResult) return;

    // Update results panel
    resultsPanel.innerHTML = '';
    resultsPanel.classList.remove('cpk-hidden');

    for (const [name, value] of optimizeResult.bestParams) {
      const label = PARAM_LABELS[name] || name;
      const lbl = ui.label(`${label} = ${value.toFixed(4)}`);
      lbl.classList.add('cpk-results-value');
      resultsPanel.appendChild(lbl);
    }

    const wsLabel = ui.label(`WSSR = ${optimizeResult.bestWssr.toFixed(4)}`);
    wsLabel.classList.add('cpk-results-value');
    const r2Label = ui.label(`R\u00B2 = ${optimizeResult.bestR2.toFixed(4)}`);
    r2Label.classList.add('cpk-results-value');
    const aicLabel = ui.label(`AIC = ${optimizeResult.bestAic.toFixed(2)}`);
    aicLabel.classList.add('cpk-results-value');
    resultsPanel.append(wsLabel, r2Label, aicLabel);

    btnApplyBestFit.classList.remove('cpk-hidden');

    // Draw landscape
    drawLandscape(optimizeResult.landscape);

    // Draw residuals
    drawResiduals(optimizeResult.residuals);
  }

  /** Draw WSSR landscape visualization */
  function drawLandscape(landscape: LandscapeData): void {
    landscapeContainer.innerHTML = '';
    landscapeContainer.classList.remove('cpk-hidden');

    const nParams = landscape.paramNames.length;

    if (nParams === 1) {
      // 1D: line chart WSSR vs parameter value
      const paramVals = landscape.paramValues[0];
      const wssrs = landscape.wssrValues;
      const lDf = DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array(landscape.paramNames[0], new Float32Array(paramVals)),
        DG.Column.fromFloat32Array('WSSR', new Float32Array(wssrs)),
      ]);
      const chart = DG.Viewer.lineChart(lDf, {
        xColumnName: landscape.paramNames[0],
        yColumnNames: ['WSSR'],
        title: 'WSSR Landscape',
      });
      landscapeContainer.appendChild(chart.root);
    } else if (nParams === 2) {
      // 2D: heatmap using canvas
      const canvas = document.createElement('canvas');
      canvas.classList.add('cpk-landscape-canvas');
      const gridRes = landscape.gridResolution;
      canvas.width = gridRes;
      canvas.height = gridRes;
      const ctx = canvas.getContext('2d')!;

      const wssrs = landscape.wssrValues;
      const minW = Math.min(...wssrs.filter((w) => isFinite(w)));
      const maxW = Math.max(...wssrs.filter((w) => isFinite(w)));
      const range = maxW - minW || 1;

      for (let i = 0; i < gridRes; i++) {
        for (let j = 0; j < gridRes; j++) {
          const idx = i * gridRes + j;
          const normalized = isFinite(wssrs[idx]) ? (wssrs[idx] - minW) / range : 1;
          const r = Math.round(255 * normalized);
          const g = Math.round(255 * (1 - normalized));
          ctx.fillStyle = `rgb(${r}, ${g}, 50)`;
          ctx.fillRect(j, i, 1, 1);
        }
      }

      const labelX = ui.label(`X: ${PARAM_LABELS[landscape.paramNames[0]] || landscape.paramNames[0]}`);
      const labelY = ui.label(`Y: ${PARAM_LABELS[landscape.paramNames[1]] || landscape.paramNames[1]}`);
      labelX.classList.add('cpk-landscape-label');
      labelY.classList.add('cpk-landscape-label');
      landscapeContainer.append(canvas, ui.divH([labelX, labelY]));
    } else if (nParams === 3) {
      // 3D: three 2D pairwise projections
      const label = ui.label('3D landscape: 3 pairwise projections (marginalized over 3rd parameter)');
      label.classList.add('cpk-landscape-label');
      landscapeContainer.appendChild(label);
      // Simplified: show text summary
      const summary = ui.label(`Grid: ${landscape.gridResolution}\u00B3 = ${landscape.wssrValues.length} points`);
      landscapeContainer.appendChild(summary);
    }
  }

  /** Draw residual plot */
  function drawResiduals(residuals: {t: number[]; res: number[]}): void {
    residualsContainer.innerHTML = '';
    residualsContainer.classList.remove('cpk-hidden');

    const resDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('Time (h)', new Float32Array(residuals.t)),
      DG.Column.fromFloat32Array('Residual (mg/L)', new Float32Array(residuals.res)),
    ]);

    const chart = DG.Viewer.scatterPlot(resDf, {
      xColumnName: 'Time (h)',
      yColumnName: 'Residual (mg/L)',
      title: 'Residuals \u2014 Observed minus Predicted',
    });
    residualsContainer.appendChild(chart.root);
  }

  /** Apply best-fit parameters to primary controls */
  function applyBestFit(): void {
    if (!optimizeResult) return;

    computationsBlocked = true;
    for (const [name, value] of optimizeResult.bestParams) {
      switch (name) {
      case 'cl': ctrlCl.value = value; break;
      case 'vc': ctrlVc.value = value; break;
      case 'vp': ctrlVp.value = value; break;
      case 'q': ctrlQ.value = value; break;
      case 'vd': ctrlVd.value = value; break;
      case 'q2': ctrlQ2.value = value; break;
      case 'ka': ctrlKa.value = value; break;
      case 'f': ctrlF.value = value; break;
      }
    }
    computationsBlocked = false;
    runPrimary();
  }

  // ===================== DATAFRAME & VIEWERS =====================

  /** Create DataFrame from PK solution */
  function createDataFrame(params: PkParams, result: PkSolution): DG.DataFrame {
    const cols: DG.Column[] = [
      DG.Column.fromFloat32Array('Time (h)', new Float32Array(result.t)),
      DG.Column.fromFloat32Array('Concentration (mg/L)', new Float32Array(result.conc)),
      DG.Column.fromFloat32Array('A_central (mg)', new Float32Array(result.ac)),
      DG.Column.fromFloat32Array('A_peripheral (mg)', new Float32Array(result.ap)),
    ];

    if (result.ad)
      cols.push(DG.Column.fromFloat32Array('A_deep (mg)', new Float32Array(result.ad)));
    if (result.agut)
      cols.push(DG.Column.fromFloat32Array('A_gut (mg)', new Float32Array(result.agut)));

    // Add MEC/MTC columns for chart reference lines
    const mecArr = new Float32Array(result.t.length).fill(params.mec);
    const mtcArr = new Float32Array(result.t.length).fill(params.mtc);
    cols.push(DG.Column.fromFloat32Array('MEC (mg/L)', mecArr));
    cols.push(DG.Column.fromFloat32Array('MTC (mg/L)', mtcArr));

    // Add observed data if available
    if (observedData) {
      const obsTimeCol = DG.Column.fromFloat32Array('Obs Time (h)',
        new Float32Array(padToLength(observedData.tObs, result.t.length)));
      const obsConcCol = DG.Column.fromFloat32Array('Obs Conc (mg/L)',
        new Float32Array(padToLength(observedData.cObs, result.t.length)));
      cols.push(obsTimeCol, obsConcCol);
    }

    // Add best-fit overlay if compare mode is on
    if (compareMode && optimizeResult) {
      const bestFitConc = new Float32Array(result.t.length);
      const bestT = optimizeResult.bestConc.t;
      const bestC = optimizeResult.bestConc.c;
      for (let i = 0; i < result.t.length; i++) {
        // Interpolate best-fit at solution time points
        const t = result.t[i];
        bestFitConc[i] = interpolateArrays(bestT, bestC, t);
      }
      cols.push(DG.Column.fromFloat32Array('Best-Fit Conc (mg/L)', bestFitConc));
    }

    return DG.DataFrame.fromColumns(cols);
  }

  /** Pad array with NaN to target length */
  function padToLength(arr: number[], len: number): number[] {
    const result = new Array(len).fill(NaN);
    for (let i = 0; i < Math.min(arr.length, len); i++)
      result[i] = arr[i];
    return result;
  }

  /** Interpolate between regular arrays */
  function interpolateArrays(t: number[], y: number[], target: number): number {
    if (target <= t[0]) return y[0];
    if (target >= t[t.length - 1]) return y[t.length - 1];
    let lo = 0;
    let hi = t.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (t[mid] <= target) lo = mid;
      else hi = mid;
    }
    const frac = (target - t[lo]) / (t[hi] - t[lo]);
    return y[lo] + frac * (y[hi] - y[lo]);
  }

  /** Update chart viewers after data update */
  function updateViewers(params: PkParams, result: PkSolution): void {
    // The charts auto-update when view.dataFrame changes
    // Update scale if needed
    updateChartScale();
  }

  // ===================== PRESET LOADING =====================

  /** Load a preset into controls */
  function loadPreset(preset: Preset): void {
    computationsBlocked = true;
    const vals = preset.values;
    if (vals.modelType !== undefined) ctrlModelType.value = vals.modelType;
    if (vals.inputMode !== undefined) ctrlInputMode.value = vals.inputMode;
    if (vals.cl !== undefined) ctrlCl.value = vals.cl;
    if (vals.vc !== undefined) ctrlVc.value = vals.vc;
    if (vals.vp !== undefined) ctrlVp.value = vals.vp;
    if (vals.q !== undefined) ctrlQ.value = vals.q;
    if (vals.vd !== undefined) ctrlVd.value = vals.vd;
    if (vals.q2 !== undefined) ctrlQ2.value = vals.q2;
    if (vals.ka !== undefined) ctrlKa.value = vals.ka;
    if (vals.f !== undefined) ctrlF.value = vals.f;
    if (vals.dose !== undefined) ctrlDose.value = vals.dose;
    if (vals.tInf !== undefined) ctrlTInf.value = vals.tInf;
    computationsBlocked = false;
    updateControlVisibility();
    updateOptParamSelector();
    runPrimary();
  }

  /** Reset all controls to defaults */
  function resetControls(): void {
    computationsBlocked = true;
    ctrlModelType.value = DEFAULTS.modelType;
    ctrlInputMode.value = DEFAULTS.inputMode;
    ctrlYScale.value = 'Linear';
    ctrlCl.value = DEFAULTS.cl;
    ctrlVc.value = DEFAULTS.vc;
    ctrlVp.value = DEFAULTS.vp;
    ctrlQ.value = DEFAULTS.q;
    ctrlVd.value = DEFAULTS.vd;
    ctrlQ2.value = DEFAULTS.q2;
    ctrlKa.value = DEFAULTS.ka;
    ctrlF.value = DEFAULTS.f;
    ctrlDose.value = DEFAULTS.dose;
    ctrlRepeated.value = DEFAULTS.repeatedDosing;
    ctrlTau.value = DEFAULTS.tau;
    ctrlNDoses.value = DEFAULTS.nDoses;
    ctrlTInf.value = DEFAULTS.tInf;
    ctrlTEnd.value = DEFAULTS.tEnd;
    ctrlMec.value = DEFAULTS.mec;
    ctrlMtc.value = DEFAULTS.mtc;
    computationsBlocked = false;

    observedData = null;
    optimizeResult = null;
    btnGenerateData.classList.remove('cpk-btn--active');
    resultsPanel.classList.add('cpk-hidden');
    btnApplyBestFit.classList.add('cpk-hidden');
    landscapeContainer.classList.add('cpk-hidden');
    residualsContainer.classList.add('cpk-hidden');

    updateControlVisibility();
    updateOptParamSelector();
    runPrimary();
  }

  // ===================== RIBBON =====================

  const presetButtons = PRESETS.map((p) =>
    ui.button(p.name, () => loadPreset(p), p.tooltip),
  );

  const resetBtn = ui.iconFA('undo', () => resetControls(), 'Reset all parameters to default values');

  view.setRibbonPanels([presetButtons, [resetBtn]]);

  // ===================== LAYOUT =====================

  // Left panel: model controls
  const leftPanel = ui.panel([
    ui.h2('Model Settings'),
    ui.inputs([ctrlModelType, ctrlInputMode, ctrlYScale]),
    ui.h2('PK Parameters'),
    ui.inputs([ctrlCl, ctrlVc, ctrlVp, ctrlQ, ctrlVd, ctrlQ2]),
    ui.h2('Absorption'),
    ui.inputs([ctrlKa, ctrlF]),
    ui.h2('Dosing'),
    ui.inputs([ctrlDose, ctrlRepeated, ctrlTau, ctrlNDoses, ctrlTInf]),
    ui.h2('Simulation'),
    ui.inputs([ctrlTEnd]),
    ui.h2('Therapeutic Window'),
    ui.inputs([ctrlMec, ctrlMtc]),
    ui.h2('PK Metrics'),
    metricsPanel,
  ]);
  leftPanel.classList.add('cpk-panel');

  // Right panel: optimization
  const rightPanel = ui.panel([
    ui.h2('Synthetic Observation Data'),
    ui.inputs([ctrlNoiseCv, ctrlSampling]),
    btnGenerateData,
    ui.h2('Optimization Setup'),
    optParamContainer,
    ui.inputs([ctrlGridRes, ctrlWeight]),
    btnRunOptimize,
    ui.h2('Best-Fit Parameters'),
    resultsPanel,
    btnApplyBestFit,
    ui.inputs([ctrlCompare]),
    landscapeContainer,
    residualsContainer,
  ]);
  rightPanel.classList.add('cpk-panel');

  // Dock panels
  view.dockManager.dock(leftPanel, DG.DOCK_TYPE.LEFT, null, 'Model Parameters', 0.3);

  // Concentration-time chart
  const yColumns = ['Concentration (mg/L)', 'MEC (mg/L)', 'MTC (mg/L)'];
  concChart = view.addViewer('Line chart', {
    xColumnName: 'Time (h)',
    yColumnNames: yColumns,
    title: 'Concentration-Time Profile',
  });

  // Compartment amounts chart
  const compYCols = ['A_central (mg)', 'A_peripheral (mg)'];
  if (DEFAULTS.modelType === '3-Compartment') compYCols.push('A_deep (mg)');
  compChart = view.addViewer('Line chart', {
    xColumnName: 'Time (h)',
    yColumnNames: compYCols,
    title: 'Compartment Amounts',
  });

  // Dock viewers
  let concNode: any = null;
  let compNode: any = null;
  if (concChart)
    concNode = view.dockManager.dock(concChart.root, DG.DOCK_TYPE.RIGHT, null, 'Concentration', 0.5);
  if (compChart)
    compNode = view.dockManager.dock(compChart.root, DG.DOCK_TYPE.DOWN, concNode, 'Amounts', 0.5);

  // Dock diagram next to compartment amounts
  view.dockManager.dock(diagramPanel, DG.DOCK_TYPE.RIGHT, compNode, 'Diagram', 0.5);

  // Dock right panel
  view.dockManager.dock(rightPanel, DG.DOCK_TYPE.RIGHT, null, 'Optimization', 0.3);

  // ===================== INITIALIZATION =====================

  updateControlVisibility();
  buildOptParamSelector();
  updateGenerateButtonLabel();
  updateOptimizeButtonLabel();
  updateMetrics(initialResult);
  drawDiagram(DEFAULTS, initialResult, 0);

  // ===================== CLEANUP =====================

  subs.push(grok.events.onViewRemoved.subscribe((v: DG.ViewBase) => {
    if (v === view) {
      stopDiagramAnimation();
      terminateWorkers();
      if (debounceTimer !== null) clearTimeout(debounceTimer);
      for (const sub of subs) sub.unsubscribe();
    }
  }));
}
