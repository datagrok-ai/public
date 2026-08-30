import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {mrt, ODEs} from 'diff-grok';

// ── Types ─────────────────────────────────────────────────────────────────────

interface LVParams {
  alpha: number;  // prey birth rate
  beta: number;   // predation rate
  delta: number;  // predator growth efficiency
  gamma: number;  // predator death rate
  x0: number;     // initial prey
  y0: number;     // initial predators
  T: number;      // simulation time
}

// ── Solver ────────────────────────────────────────────────────────────────────

function solveLV(p: LVParams): {t: Float64Array; x: Float64Array; y: Float64Array} {
  const step = Math.max(0.05, p.T / 2000);
  const odes: ODEs = {
    name: 'LotkaVolterra',
    arg: {name: 't', start: 0, finish: p.T, step},
    initial: [p.x0, p.y0],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      out[0] = p.alpha * y[0] - p.beta * y[0] * y[1];
      out[1] = p.delta * y[0] * y[1] - p.gamma * y[1];
    },
    tolerance: 1e-6,
    solutionColNames: ['prey', 'predators'],
  };
  const sol = mrt(odes);
  return {t: sol[0], x: sol[1], y: sol[2]};
}

// ── DataFrame factories ────────────────────────────────────────────────────────

function makeTimeSeriesDf(t: Float64Array, x: Float64Array, y: Float64Array): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromFloat64Array('t', t.slice()),
    DG.Column.fromFloat64Array('prey', x.slice()),
    DG.Column.fromFloat64Array('predators', y.slice()),
  ]);
}

function makePhaseDf(x: Float64Array, y: Float64Array): DG.DataFrame {
  const n = x.length;
  const labels = new Array<string>(n).fill('trajectory');
  labels[0] = 'start';
  return DG.DataFrame.fromColumns([
    DG.Column.fromFloat64Array('prey', x.slice()),
    DG.Column.fromFloat64Array('predators', y.slice()),
    DG.Column.fromStrings('point', labels),
  ]);
}

function makeTableDf(t: Float64Array, x: Float64Array, y: Float64Array): DG.DataFrame {
  const n = t.length;
  const stride = Math.max(1, Math.floor(n / 200));
  const rows = Math.ceil(n / stride);
  const st = new Float64Array(rows);
  const sx = new Float64Array(rows);
  const sy = new Float64Array(rows);
  for (let r = 0; r < rows; r++) {
    const i = Math.min(r * stride, n - 1);
    st[r] = t[i]; sx[r] = x[i]; sy[r] = y[i];
  }
  return DG.DataFrame.fromColumns([
    DG.Column.fromFloat64Array('t', st),
    DG.Column.fromFloat64Array('x (prey)', sx),
    DG.Column.fromFloat64Array('y (pred)', sy),
  ]);
}

// ── Slider definitions ─────────────────────────────────────────────────────────

interface SliderDef {
  key: keyof LVParams;
  label: string;
  min: number; max: number; step: number; decimals: number;
  tooltip: string;
}

const SLIDER_DEFS: SliderDef[] = [
  {
    key: 'alpha', label: 'α — prey birth rate',
    min: 0.1, max: 3.0, step: 0.01, decimals: 2,
    tooltip: 'Prey birth rate (α): how fast prey reproduce without predators.\n' +
      'Higher α → larger oscillation amplitude, higher y* = α/β equilibrium.',
  },
  {
    key: 'beta', label: 'β — predation rate',
    min: 0.01, max: 0.5, step: 0.001, decimals: 3,
    tooltip: 'Predation rate (β): likelihood of a predator catching prey per encounter.\n' +
      'Higher β → prey equilibrium x* = γ/δ unchanged but oscillations shift.',
  },
  {
    key: 'delta', label: 'δ — predator efficiency',
    min: 0.01, max: 0.5, step: 0.001, decimals: 3,
    tooltip: 'Predator growth efficiency (δ): fraction of eaten prey converted to new predators.\n' +
      'Higher δ → lower prey equilibrium x* = γ/δ, more predators.',
  },
  {
    key: 'gamma', label: 'γ — predator death rate',
    min: 0.1, max: 3.0, step: 0.01, decimals: 2,
    tooltip: 'Predator death rate (γ): natural mortality of predators.\n' +
      'Higher γ → fewer predators, higher prey equilibrium x* = γ/δ.',
  },
  {
    key: 'x0', label: 'x₀ — initial prey',
    min: 1, max: 100, step: 1, decimals: 0,
    tooltip: 'Initial prey population (x₀): starting point on the phase portrait.\n' +
      'Near x* = γ/δ → small oscillations; far away → large swings.',
  },
  {
    key: 'y0', label: 'y₀ — initial predators',
    min: 1, max: 50, step: 1, decimals: 0,
    tooltip: 'Initial predator population (y₀): starting point on the phase portrait.\n' +
      'Near y* = α/β → small oscillations; far away → large swings.',
  },
  {
    key: 'T', label: 'T — simulation time',
    min: 10, max: 500, step: 5, decimals: 0,
    tooltip: 'Total simulation time T. Increase to observe more oscillation cycles.\n' +
      'Lotka-Volterra is conservative: cycle period ≈ 2π / √(α·γ).',
  },
];

// ── Slider widget factory ──────────────────────────────────────────────────────

function makeSlider(
  def: SliderDef,
  value: number,
  onChange: (v: number) => void,
): {root: HTMLElement; setValue: (v: number) => void} {
  const valSpan = document.createElement('span');
  valSpan.className = 'lv-val';
  valSpan.textContent = value.toFixed(def.decimals);

  const labelRow = document.createElement('div');
  labelRow.className = 'lv-label-row';
  labelRow.appendChild(document.createTextNode(def.label));
  labelRow.appendChild(valSpan);

  const range = document.createElement('input');
  range.type = 'range';
  range.className = 'lv-range';
  range.min = String(def.min);
  range.max = String(def.max);
  range.step = String(def.step);
  range.value = String(value);
  range.addEventListener('input', () => {
    const v = parseFloat(range.value);
    valSpan.textContent = v.toFixed(def.decimals);
    onChange(v);
  });

  const root = document.createElement('div');
  root.className = 'lv-slider-row';
  root.title = def.tooltip;
  root.appendChild(labelRow);
  root.appendChild(range);

  return {
    root,
    setValue(v: number) {
      range.value = String(v);
      valSpan.textContent = v.toFixed(def.decimals);
    },
  };
}

// ── CSS ────────────────────────────────────────────────────────────────────────

const CSS = `
.lv-root {
  display: flex; height: 100%; width: 100%; overflow: hidden; font-size: 12px; box-sizing: border-box;
}
.lv-left {
  width: 274px; min-width: 274px; display: flex; flex-direction: column; gap: 6px;
  padding: 10px 10px 10px 12px; overflow-y: auto; border-right: 1px solid var(--grey-2,#e0e0e0);
}
.lv-panel-title {
  font-size: 13px; font-weight: 600; color: var(--grey-8,#222); margin-bottom: 2px;
}
.lv-group-label {
  font-size: 11px; font-weight: 600; color: var(--blue-2,#1565c0); text-transform: uppercase;
  letter-spacing: 0.05em; margin-top: 4px; margin-bottom: 2px;
  border-bottom: 1px solid var(--grey-2,#e0e0e0); padding-bottom: 2px;
}
.lv-slider-row { padding: 3px 0; cursor: default; }
.lv-label-row { display: flex; justify-content: space-between; margin-bottom: 2px; color: var(--grey-6,#444); }
.lv-val { font-weight: 700; color: var(--blue-1,#1976d2); min-width: 38px; text-align: right; }
.lv-range { width: 100%; cursor: pointer; accent-color: var(--blue-1,#1976d2); }
.lv-stats-block { border-top: 1px solid var(--grey-2,#e0e0e0); padding-top: 6px; }
.lv-stats-title { font-weight: 600; color: var(--grey-7,#333); margin-bottom: 3px; }
.lv-stats-line { color: var(--grey-6,#555); line-height: 1.6; font-family: monospace; font-size: 11.5px; }
.lv-progress-wrap {
  display: none; height: 18px; background: var(--grey-1,#f0f0f0); border-radius: 9px;
  overflow: hidden; position: relative; margin-top: 2px;
}
.lv-progress-bar { height: 100%; background: var(--blue-1,#1976d2); transition: width 0.08s linear; width: 0; }
.lv-progress-pct {
  position: absolute; right: 8px; top: 0; line-height: 18px; font-size: 11px; font-weight: 600;
  color: var(--grey-7,#333);
}
.lv-center {
  flex: 1; min-width: 0; display: flex; flex-direction: column;
}
.lv-chart-section {
  flex: 1; min-height: 0; display: flex; flex-direction: column;
}
.lv-section-title {
  font-size: 11px; font-weight: 600; color: var(--grey-6,#555); text-transform: uppercase;
  letter-spacing: 0.06em; padding: 4px 10px 0; flex-shrink: 0;
}
.lv-chart-wrap {
  flex: 1; min-height: 0; position: relative; overflow: hidden;
}
.lv-chart-wrap > * { position: absolute !important; inset: 0 !important; width: 100% !important; height: 100% !important; }
.lv-divider { height: 1px; background: var(--grey-2,#e0e0e0); flex-shrink: 0; }
.lv-right {
  width: 230px; min-width: 230px; border-left: 1px solid var(--grey-2,#e0e0e0);
  display: flex; flex-direction: column; overflow: hidden;
}
.lv-right-title {
  font-size: 11px; font-weight: 600; color: var(--grey-6,#555); text-transform: uppercase;
  letter-spacing: 0.06em; padding: 4px 10px 0; flex-shrink: 0;
}
.lv-grid-wrap { flex: 1; min-height: 0; overflow: hidden; }
.lv-grid-wrap > * { width: 100% !important; height: 100% !important; }
`;

// ── App entry point ────────────────────────────────────────────────────────────

export function runLotkaVolterra(): void {
  // ── State ──────────────────────────────────────────────────────────────────
  const params: LVParams = {alpha: 1.0, beta: 0.1, delta: 0.075, gamma: 1.5, x0: 10, y0: 5, T: 100};

  // ── Initial solve ──────────────────────────────────────────────────────────
  let sol = solveLV(params);
  let tsDf = makeTimeSeriesDf(sol.t, sol.x, sol.y);
  let phaseDf = makePhaseDf(sol.x, sol.y);
  let tableDf = makeTableDf(sol.t, sol.x, sol.y);

  // ── Viewers ────────────────────────────────────────────────────────────────
  const lineChart = DG.Viewer.fromType('Line chart', tsDf, {xColumnName: 't'});
  const phaseScatter = DG.Viewer.fromType('Scatter plot', phaseDf, {
    xColumnName: 'prey',
    yColumnName: 'predators',
    colorColumnName: 'point',
  });
  const gridViewer = DG.Viewer.fromType('Grid', tableDf);

  // ── Stats ──────────────────────────────────────────────────────────────────
  const equilLine = document.createElement('div');
  equilLine.className = 'lv-stats-line';
  const summaryLine = document.createElement('div');
  summaryLine.className = 'lv-stats-line';

  function updateStats(): void {
    const xStar = params.gamma / params.delta;
    const yStar = params.alpha / params.beta;
    equilLine.textContent = `x* = ${xStar.toFixed(2)},  y* = ${yStar.toFixed(2)}`;
    let maxX = -Infinity;
    let maxY = -Infinity;
    for (let i = 0; i < sol.x.length; i++) { if (sol.x[i] > maxX) maxX = sol.x[i]; }
    for (let i = 0; i < sol.y.length; i++) { if (sol.y[i] > maxY) maxY = sol.y[i]; }
    summaryLine.textContent = `max prey = ${maxX.toFixed(1)}\nmax pred = ${maxY.toFixed(1)}\nsteps    = ${sol.t.length}`;
  }
  updateStats();

  // ── Debounced solver update ────────────────────────────────────────────────
  let debounceTimer = 0;
  function scheduleUpdate(): void {
    clearTimeout(debounceTimer);
    debounceTimer = window.setTimeout(() => {
      try {
        sol = solveLV(params);
        tsDf = makeTimeSeriesDf(sol.t, sol.x, sol.y);
        phaseDf = makePhaseDf(sol.x, sol.y);
        tableDf = makeTableDf(sol.t, sol.x, sol.y);
        lineChart.dataFrame = tsDf;
        phaseScatter.dataFrame = phaseDf;
        gridViewer.dataFrame = tableDf;
        updateStats();
      } catch (err) {
        console.error('[LotkaVolterra] solver error:', err);
      }
    }, 60);
  }

  // ── Sliders ────────────────────────────────────────────────────────────────
  const sliderRefs = new Map<keyof LVParams, {setValue: (v: number) => void}>();
  const modelGroup = document.createElement('div');
  modelGroup.className = 'lv-group-label';
  modelGroup.textContent = 'Model Coefficients';
  const initGroup = document.createElement('div');
  initGroup.className = 'lv-group-label';
  initGroup.textContent = 'Initial Conditions';

  const slidersDiv = document.createElement('div');
  slidersDiv.appendChild(modelGroup);
  for (const def of SLIDER_DEFS.slice(0, 4)) {
    const {root, setValue} = makeSlider(def, params[def.key], (v) => {
      (params as unknown as Record<string, number>)[def.key] = v;
      scheduleUpdate();
    });
    sliderRefs.set(def.key, {setValue});
    slidersDiv.appendChild(root);
  }
  slidersDiv.appendChild(initGroup);
  for (const def of SLIDER_DEFS.slice(4)) {
    const {root, setValue} = makeSlider(def, params[def.key], (v) => {
      (params as unknown as Record<string, number>)[def.key] = v;
      scheduleUpdate();
    });
    sliderRefs.set(def.key, {setValue});
    slidersDiv.appendChild(root);
  }

  // ── Equilibrium & stats blocks ─────────────────────────────────────────────
  const equilBlock = document.createElement('div');
  equilBlock.className = 'lv-stats-block';
  const equilTitle = document.createElement('div');
  equilTitle.className = 'lv-stats-title';
  equilTitle.textContent = 'Equilibrium';
  equilBlock.appendChild(equilTitle);
  equilBlock.appendChild(equilLine);

  const summaryBlock = document.createElement('div');
  summaryBlock.className = 'lv-stats-block';
  const summaryTitle = document.createElement('div');
  summaryTitle.className = 'lv-stats-title';
  summaryTitle.textContent = 'Summary';
  summaryBlock.appendChild(summaryTitle);
  summaryBlock.appendChild(summaryLine);

  // ── Optimizer button & progress bar ───────────────────────────────────────
  const progressWrap = document.createElement('div');
  progressWrap.className = 'lv-progress-wrap';
  const progressBar = document.createElement('div');
  progressBar.className = 'lv-progress-bar';
  const progressPct = document.createElement('span');
  progressPct.className = 'lv-progress-pct';
  progressPct.textContent = '0%';
  progressWrap.appendChild(progressBar);
  progressWrap.appendChild(progressPct);

  let worker: Worker | null = null;
  const optimizeBtn = ui.bigButton('Optimize Max Prey', () => {
    if (worker) {
      worker.terminate();
      worker = null;
      optimizeBtn.textContent = 'Optimize Max Prey';
      progressWrap.style.display = 'none';
      return;
    }

    progressWrap.style.display = 'block';
    progressBar.style.width = '0%';
    progressPct.textContent = '0%';
    optimizeBtn.textContent = 'Cancel';

    worker = new Worker(new URL('./optimizer.worker.ts', import.meta.url));

    worker.postMessage({
      alphaRange: [0.1, 3.0] as [number, number],
      betaRange:  [0.01, 0.5] as [number, number],
      deltaRange: [0.01, 0.5] as [number, number],
      gammaRange: [0.1, 3.0] as [number, number],
      x0: params.x0,
      y0: params.y0,
      T: Math.min(params.T, 100),
    });

    worker.onmessage = (ev: MessageEvent) => {
      if (ev.data.type === 'progress') {
        const pct = Math.round(ev.data.progress * 100);
        progressBar.style.width = `${pct}%`;
        progressPct.textContent = `${pct}%`;
      } else if (ev.data.type === 'result') {
        const {alpha, beta, delta, gamma} = ev.data;
        params.alpha = alpha; params.beta = beta;
        params.delta = delta; params.gamma = gamma;
        sliderRefs.get('alpha')!.setValue(alpha);
        sliderRefs.get('beta')!.setValue(beta);
        sliderRefs.get('delta')!.setValue(delta);
        sliderRefs.get('gamma')!.setValue(gamma);
        scheduleUpdate();
        worker = null;
        optimizeBtn.textContent = 'Optimize Max Prey';
        progressWrap.style.display = 'none';
      }
    };

    worker.onerror = (ev) => {
      console.error('[LotkaVolterra] worker error:', ev);
      worker = null;
      optimizeBtn.textContent = 'Optimize Max Prey';
      progressWrap.style.display = 'none';
    };
  });

  // ── Left panel ─────────────────────────────────────────────────────────────
  const leftPanel = document.createElement('div');
  leftPanel.className = 'lv-left';

  const panelTitle = document.createElement('div');
  panelTitle.className = 'lv-panel-title';
  panelTitle.textContent = 'Parameters';

  leftPanel.appendChild(panelTitle);
  leftPanel.appendChild(slidersDiv);
  leftPanel.appendChild(equilBlock);
  leftPanel.appendChild(summaryBlock);
  leftPanel.appendChild(progressWrap);
  leftPanel.appendChild(optimizeBtn);

  // ── Center panel ───────────────────────────────────────────────────────────
  const tsTitle = document.createElement('div');
  tsTitle.className = 'lv-section-title';
  tsTitle.textContent = 'Population Dynamics';
  const tsWrap = document.createElement('div');
  tsWrap.className = 'lv-chart-wrap';
  tsWrap.appendChild(lineChart.root);
  const tsSection = document.createElement('div');
  tsSection.className = 'lv-chart-section';
  tsSection.appendChild(tsTitle);
  tsSection.appendChild(tsWrap);

  const divider = document.createElement('div');
  divider.className = 'lv-divider';

  const phaseTitle = document.createElement('div');
  phaseTitle.className = 'lv-section-title';
  phaseTitle.textContent = 'Phase Portrait  (● = start)';
  const phaseWrap = document.createElement('div');
  phaseWrap.className = 'lv-chart-wrap';
  phaseWrap.appendChild(phaseScatter.root);
  const phaseSection = document.createElement('div');
  phaseSection.className = 'lv-chart-section';
  phaseSection.appendChild(phaseTitle);
  phaseSection.appendChild(phaseWrap);

  const centerPanel = document.createElement('div');
  centerPanel.className = 'lv-center';
  centerPanel.appendChild(tsSection);
  centerPanel.appendChild(divider);
  centerPanel.appendChild(phaseSection);

  // ── Right panel ────────────────────────────────────────────────────────────
  const rightTitle = document.createElement('div');
  rightTitle.className = 'lv-right-title';
  rightTitle.textContent = 'Data (sampled)';
  const gridWrap = document.createElement('div');
  gridWrap.className = 'lv-grid-wrap';
  gridWrap.appendChild(gridViewer.root);

  const rightPanel = document.createElement('div');
  rightPanel.className = 'lv-right';
  rightPanel.appendChild(rightTitle);
  rightPanel.appendChild(gridWrap);

  // ── Root ───────────────────────────────────────────────────────────────────
  const appRoot = document.createElement('div');
  appRoot.className = 'lv-root';
  appRoot.appendChild(leftPanel);
  appRoot.appendChild(centerPanel);
  appRoot.appendChild(rightPanel);

  // Inject CSS once
  if (!document.getElementById('lv-styles')) {
    const style = document.createElement('style');
    style.id = 'lv-styles';
    style.textContent = CSS;
    document.head.appendChild(style);
  }

  // ── Open view ──────────────────────────────────────────────────────────────
  const view = grok.shell.newView('Lotka–Volterra', [appRoot]);
  view.root.style.cssText = 'padding:0;overflow:hidden;display:flex;height:100%;';

  // Clean up worker if view is closed
  (view as any).onClosed?.subscribe(() => { if (worker) { worker.terminate(); worker = null; } });
}
