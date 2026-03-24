// CR3BP — Application (Coordinator + UI)

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULTS, RANGES, validate, solve,
  CR3BPParams, CR3BPSolution, InputId,
} from './core';

import '../../css/cr3bp.css';

const DEBOUNCE_MS = 100;
const ROUTH_CRITERION = 0.0385;

export function threeBodyProblemApp(_package: DG.Package): void {
  // --- State ---
  let computationsBlocked = false;
  let debounceTimer: ReturnType<typeof setTimeout> | null = null;
  const subs: {unsubscribe(): void}[] = [];
  let orbitViewer!: DG.ScatterPlotViewer;
  let jacobiViewer!: DG.Viewer;
  let zvcViewer!: DG.ScatterPlotViewer;
  let currentSolution: CR3BPSolution | null = null;

  // --- Initial computation ---
  const initSolution = solve(DEFAULTS);
  currentSolution = initSolution;

  const df = createTrajectoryDataFrame(initSolution);
  let zvcDf = createZVCDataFrame(initSolution);

  const view = grok.shell.addTableView(df);
  view.name = 'Circular Restricted Three-Body Problem';

  // --- Lagrange info panel ---
  const lagrangeInfoDiv = ui.divV([], 'cr3bp-app-lagrange-panel');

  function updateLagrangePanel(sol: CR3BPSolution, mu?: number): void {
    lagrangeInfoDiv.innerHTML = '';
    const muVal = mu ?? (ctrlMu ? (ctrlMu.value ?? DEFAULTS.mu) : DEFAULTS.mu);
    const isStable = muVal < ROUTH_CRITERION;

    for (const lp of sol.lagrangePoints) {
      const entry = ui.divH([
        ui.label(`${lp.name}: (${lp.x.toFixed(4)}, ${lp.y.toFixed(4)})`),
      ], 'cr3bp-app-lagrange-entry');
      lagrangeInfoDiv.append(entry);
    }

    const stabilityNote = ui.label(
      isStable
        ? `L4, L5 stable (\u03BC < ${ROUTH_CRITERION})`
        : `L4, L5 unstable (\u03BC \u2265 ${ROUTH_CRITERION})`,
    );
    stabilityNote.classList.add(isStable ? 'cr3bp-app-lagrange-stable' : 'cr3bp-app-lagrange-unstable');
    lagrangeInfoDiv.append(stabilityNote);

    const collinearNote = ui.label('L1, L2, L3 always unstable');
    collinearNote.classList.add('cr3bp-app-lagrange-unstable');
    lagrangeInfoDiv.append(collinearNote);
  }
  updateLagrangePanel(initSolution, DEFAULTS.mu);

  // --- Controls ---

  // System Parameters
  const ctrlMu = ui.input.float('\u03BC \u2014 mass ratio m\u2082/(m\u2081+m\u2082)', {
    value: DEFAULTS.mu, nullable: false,
    min: RANGES.mu.min, max: RANGES.mu.max,
    tooltipText: 'Ratio of the smaller body\'s mass to the total. For Earth\u2013Moon, \u03BC \u2248 0.012. Increasing \u03BC enlarges the gravitational sphere of influence of the second body \u2014 the topology of forbidden regions changes.',
    onValueChanged: () => debouncedRun(),
  });

  // Initial Conditions — Position
  const ctrlX0 = ui.input.float('Initial position x\u2080', {
    value: DEFAULTS.x0, nullable: false,
    min: RANGES.x0.min, max: RANGES.x0.max,
    tooltipText: 'Starting x coordinate of the test body in the rotating frame. Can also be set by clicking on the orbit plot.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlY0 = ui.input.float('Initial position y\u2080', {
    value: DEFAULTS.y0, nullable: false,
    min: RANGES.y0.min, max: RANGES.y0.max,
    tooltipText: 'Starting y coordinate of the test body in the rotating frame. Can also be set by clicking on the orbit plot.',
    onValueChanged: () => debouncedRun(),
  });

  // Initial Conditions — Velocity
  const ctrlVx0 = ui.input.float('Initial velocity vx\u2080 (rotating frame)', {
    value: DEFAULTS.vx0, nullable: false,
    min: RANGES.vx0.min, max: RANGES.vx0.max,
    tooltipText: 'Initial velocity in the x direction. Can be set by dragging from the click point on the orbit plot.',
    onValueChanged: () => debouncedRun(),
  });

  const ctrlVy0 = ui.input.float('Initial velocity vy\u2080 (rotating frame)', {
    value: DEFAULTS.vy0, nullable: false,
    min: RANGES.vy0.min, max: RANGES.vy0.max,
    tooltipText: 'Initial velocity in the y direction. Can be set by dragging from the click point on the orbit plot.',
    onValueChanged: () => debouncedRun(),
  });

  // Integration
  const ctrlT = ui.input.float('Integration time (orbital periods)', {
    value: DEFAULTS.T, nullable: false,
    min: RANGES.T.min, max: RANGES.T.max,
    tooltipText: 'Total integration time in dimensionless units. One lunar orbital period \u2248 2\u03C0 \u2248 6.28 time units.',
    onValueChanged: () => debouncedRun(),
  });

  // Set formats (block computations to avoid spurious runs from format-triggered events)
  computationsBlocked = true;
  ctrlMu.format = '0.00000';
  ctrlX0.format = '0.000';
  ctrlY0.format = '0.000';
  ctrlVx0.format = '0.000';
  ctrlVy0.format = '0.000';
  ctrlT.format = '0.0';
  computationsBlocked = false;

  // --- Input map for validators ---
  const inputMap: Record<InputId, DG.InputBase> = {
    'ctrl_mu': ctrlMu,
    'ctrl_x0': ctrlX0,
    'ctrl_y0': ctrlY0,
    'ctrl_vx0': ctrlVx0,
    'ctrl_vy0': ctrlVy0,
    'ctrl_T': ctrlT,
  };

  // --- Gather current inputs ---
  function getInputs(): CR3BPParams {
    return {
      mu: ctrlMu.value ?? DEFAULTS.mu,
      x0: ctrlX0.value ?? DEFAULTS.x0,
      y0: ctrlY0.value ?? DEFAULTS.y0,
      vx0: ctrlVx0.value ?? DEFAULTS.vx0,
      vy0: ctrlVy0.value ?? DEFAULTS.vy0,
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

    ctrlMu.addValidator(validatorFor('ctrl_mu'));
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
      currentSolution = result;
      updateDataFrame(result);
      updateZVCDataFrame(result);
      updateLagrangePanel(result);
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

  // --- DataFrame creation ---
  function createTrajectoryDataFrame(result: CR3BPSolution): DG.DataFrame {
    const newDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Time', result.t),
      DG.Column.fromFloat64Array('x', result.x),
      DG.Column.fromFloat64Array('y', result.y),
      DG.Column.fromFloat64Array('vx', result.vx),
      DG.Column.fromFloat64Array('vy', result.vy),
      DG.Column.fromFloat64Array('Velocity', result.vmag),
      DG.Column.fromFloat64Array('Jacobi Constant', result.cj),
    ]);
    newDf.name = 'CR3BP Trajectory';
    return newDf;
  }

  function createZVCDataFrame(result: CR3BPSolution): DG.DataFrame {
    const grid = result.zvcGrid;
    const forbiddenCol = new Float64Array(grid.xArr.length);
    forbiddenCol.fill(1.0);
    const newDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('x_grid', grid.xArr),
      DG.Column.fromFloat64Array('y_grid', grid.yArr),
      DG.Column.fromFloat64Array('Forbidden', forbiddenCol),
    ]);
    newDf.name = 'CR3BP Forbidden Regions';
    return newDf;
  }

  function updateDataFrame(result: CR3BPSolution): void {
    const newDf = createTrajectoryDataFrame(result);
    view.dataFrame = newDf;
    orbitViewer.dataFrame = newDf;
    jacobiViewer.dataFrame = newDf;
  }

  function updateZVCDataFrame(result: CR3BPSolution): void {
    zvcDf = createZVCDataFrame(result);
    zvcViewer.dataFrame = zvcDf;
  }

  function clearResults(): void {
    const emptyDf = DG.DataFrame.fromColumns([
      DG.Column.fromFloat64Array('Time', new Float64Array(0)),
      DG.Column.fromFloat64Array('x', new Float64Array(0)),
      DG.Column.fromFloat64Array('y', new Float64Array(0)),
      DG.Column.fromFloat64Array('vx', new Float64Array(0)),
      DG.Column.fromFloat64Array('vy', new Float64Array(0)),
      DG.Column.fromFloat64Array('Velocity', new Float64Array(0)),
      DG.Column.fromFloat64Array('Jacobi Constant', new Float64Array(0)),
    ]);
    emptyDf.name = 'CR3BP Trajectory';
    view.dataFrame = emptyDf;
    orbitViewer.dataFrame = emptyDf;
    jacobiViewer.dataFrame = emptyDf;
  }

  // --- Preset handling ---
  function applyPreset(values: Partial<CR3BPParams>): void {
    computationsBlocked = true;
    if (values.mu !== undefined) ctrlMu.value = values.mu;
    if (values.x0 !== undefined) ctrlX0.value = values.x0;
    if (values.y0 !== undefined) ctrlY0.value = values.y0;
    if (values.vx0 !== undefined) ctrlVx0.value = values.vx0;
    if (values.vy0 !== undefined) ctrlVy0.value = values.vy0;
    if (values.T !== undefined) ctrlT.value = values.T;
    computationsBlocked = false;
    runPrimary();
  }

  // --- Ribbon buttons ---
  const presetEarthMoon = ui.button('Earth\u2013Moon', () => applyPreset({mu: 0.01215}),
    'Set mass parameter to the Earth\u2013Moon system value');
  const presetSunJupiter = ui.button('Sun\u2013Jupiter', () => applyPreset({mu: 0.000953}),
    'Set mass parameter to the Sun\u2013Jupiter system value');
  const presetPlutoCharon = ui.button('Pluto\u2013Charon', () => applyPreset({mu: 0.1}),
    'Set mass parameter to the Pluto\u2013Charon system value');

  const presetHalo = ui.button('L1 Halo Orbit', () => applyPreset({
    mu: 0.01215, x0: 0.8369, y0: 0, vx0: 0, vy0: -0.0559, T: 6.19,
  }), 'A periodic orbit around L1 \u2014 used by real space missions (e.g., JWST). Requires precise initial conditions.');

  const presetFreeReturn = ui.button('Free-Return', () => applyPreset({
    mu: 0.01215, x0: 0.994, y0: 0, vx0: 0, vy0: -2.0016, T: 11.12,
  }), 'A figure-eight trajectory that loops around the Moon and returns to the vicinity of Earth without propulsion. Used by Apollo missions as an abort trajectory.');

  const resetBtn = ui.iconFA('undo', () => applyPreset(DEFAULTS),
    'Reset all parameters to default values');

  view.setRibbonPanels([
    [presetEarthMoon, presetSunJupiter, presetPlutoCharon],
    [presetHalo, presetFreeReturn],
    [resetBtn],
  ]);

  // --- Layout: left panel with form ---
  const form = ui.form([]);

  form.append(ui.h2('Initial Conditions & System Parameters'));

  form.append(ui.h3('System'));
  form.append(ctrlMu.root);

  form.append(ui.h3('Position'));
  form.append(ctrlX0.root);
  form.append(ctrlY0.root);

  form.append(ui.h3('Velocity'));
  form.append(ctrlVx0.root);
  form.append(ctrlVy0.root);

  form.append(ui.h3('Integration'));
  form.append(ctrlT.root);

  form.append(ui.h2('Lagrange Points'));
  form.append(lagrangeInfoDiv);

  form.classList.add('cr3bp-app-input-panel');

  const dockMng = view.dockManager;
  dockMng.dock(form, DG.DOCK_TYPE.LEFT, null, undefined, 0.6);

  // --- Orbit scatter plot ---
  orbitViewer = view.addViewer('Scatter plot', {
    xColumnName: 'x',
    yColumnName: 'y',
    colorColumnName: 'Velocity',
    title: 'Orbit in Rotating Frame (x, y)',
    markerDefaultSize: 2,
  }) as DG.ScatterPlotViewer;

  // --- ZVC scatter plot (separate DataFrame) ---
  zvcViewer = DG.Viewer.scatterPlot(zvcDf, {
    xColumnName: 'x_grid',
    yColumnName: 'y_grid',
    title: 'Forbidden Regions (Hill surfaces at current C_J)',
    markerDefaultSize: 1,
    colorColumnName: 'Forbidden',
  }) as DG.ScatterPlotViewer;

  // --- Jacobi constant line chart ---
  jacobiViewer = view.addViewer('Line chart', {
    xColumnName: 'Time',
    yColumnNames: ['Jacobi Constant'],
    title: 'Jacobi Constant Drift (integration accuracy)',
  });

  // --- Docking: orbit top-right, ZVC bottom-left, Jacobi bottom-right ---
  dockMng.dock(orbitViewer, DG.DOCK_TYPE.RIGHT, null, undefined, 0.5);
  const orbitNode = dockMng.findNode(orbitViewer.root);
  if (orbitNode != null) {
    const zvcNode = dockMng.dock(zvcViewer, DG.DOCK_TYPE.DOWN, orbitNode, undefined, 0.5);
    if (zvcNode != null)
      dockMng.dock(jacobiViewer, DG.DOCK_TYPE.RIGHT, zvcNode, undefined, 0.5);
  }

  // --- Orbit overlay (bodies + Lagrange points) ---
  subs.push(orbitViewer.onAfterDrawScene.subscribe(() => {
    if (!currentSolution) return;
    const canvas = orbitViewer.canvas;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const mu = getInputs().mu;
    const lps = currentSolution.lagrangePoints;

    // Draw Body 1 (Earth) at (-μ, 0)
    drawMarker(ctx, orbitViewer, -mu, 0, 8, '#4169E1', 'Earth');

    // Draw Body 2 (Moon) at (1-μ, 0)
    drawMarker(ctx, orbitViewer, 1 - mu, 0, 6, '#808080', 'Moon');

    // Draw Lagrange points
    for (const lp of lps)
      drawMarker(ctx, orbitViewer, lp.x, lp.y, 4, '#FF6600', lp.name);
  }));

  // --- Click-to-launch on orbit plot ---
  let isDragging = false;
  let dragStartWorld: {x: number; y: number} | null = null;

  const setupClickToLaunch = () => {
    const canvas = orbitViewer.canvas;
    if (!canvas) return;

    const handleMouseDown = (e: MouseEvent) => {
      if (e.button !== 0) return;
      const rect = canvas.getBoundingClientRect();
      const scaleX = canvas.width / rect.width;
      const scaleY = canvas.height / rect.height;
      const screenX = (e.clientX - rect.left) * scaleX;
      const screenY = (e.clientY - rect.top) * scaleY;
      const world = orbitViewer.screenToWorld(screenX, screenY);
      if (!world) return;
      dragStartWorld = {x: world.x, y: world.y};
      isDragging = true;
    };

    const handleMouseUp = (e: MouseEvent) => {
      if (!isDragging || !dragStartWorld) return;
      isDragging = false;

      const rect = canvas.getBoundingClientRect();
      const scaleX = canvas.width / rect.width;
      const scaleY = canvas.height / rect.height;
      const screenX = (e.clientX - rect.left) * scaleX;
      const screenY = (e.clientY - rect.top) * scaleY;
      const worldEnd = orbitViewer.screenToWorld(screenX, screenY);

      computationsBlocked = true;
      ctrlX0.value = dragStartWorld.x;
      ctrlY0.value = dragStartWorld.y;
      if (worldEnd) {
        const dx = worldEnd.x - dragStartWorld.x;
        const dy = worldEnd.y - dragStartWorld.y;
        // Only set velocity if drag distance is meaningful
        if (Math.sqrt(dx * dx + dy * dy) > 0.01) {
          ctrlVx0.value = dx;
          ctrlVy0.value = dy;
        } else {
          ctrlVx0.value = 0;
          ctrlVy0.value = 0;
        }
      } else {
        ctrlVx0.value = 0;
        ctrlVy0.value = 0;
      }
      computationsBlocked = false;
      dragStartWorld = null;
      runPrimary();
    };

    canvas.addEventListener('dblclick', (e: MouseEvent) => {
      const rect = canvas.getBoundingClientRect();
      const scaleX = canvas.width / rect.width;
      const scaleY = canvas.height / rect.height;
      const screenX = (e.clientX - rect.left) * scaleX;
      const screenY = (e.clientY - rect.top) * scaleY;
      const world = orbitViewer.screenToWorld(screenX, screenY);
      if (!world) return;

      computationsBlocked = true;
      ctrlX0.value = world.x;
      ctrlY0.value = world.y;
      ctrlVx0.value = 0;
      ctrlVy0.value = 0;
      computationsBlocked = false;
      runPrimary();
    });

    canvas.addEventListener('mousedown', handleMouseDown);
    canvas.addEventListener('mouseup', handleMouseUp);
  };

  // Defer click-to-launch setup until the viewer is rendered
  setTimeout(() => setupClickToLaunch(), 500);

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

/** Draws a labeled marker on the scatter plot canvas */
function drawMarker(
  ctx: CanvasRenderingContext2D,
  viewer: DG.ScatterPlotViewer,
  worldX: number, worldY: number,
  size: number, color: string, label: string,
): void {
  const pt = viewer.worldToScreen(worldX, worldY);
  if (!pt) return;

  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(pt.x, pt.y, size, 0, 2 * Math.PI);
  ctx.fill();

  ctx.fillStyle = color;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  ctx.fillText(label, pt.x + size + 2, pt.y - 2);
}
