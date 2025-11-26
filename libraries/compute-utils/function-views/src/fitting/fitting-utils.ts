/* eslint-disable valid-jsdoc */
// Fitting utilities

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Extremum, InconsistentTables, ValueBoundsData} from './optimizer-misc';

import '../../css/fitting-view.css';
import '../../css/sens-analysis.css';
import {GRID_SIZE, HELP_LINK, INDICES, TIMEOUT, TARGET_DATAFRAME_INFO, TITLE, NAME,
  MIN_RADAR_COLS_COUNT, RAND_MOD, RAND_MULT, ReproSettings,
  REPRO_DEFAULT,
  SEED_DEFAULT,
  EarlyStoppingSettings,
  EARLY_STOP_DEFAULT,
  COST_FUNC_THRESH,
  STOP_AFTER_DEFAULT,
  STOP_AFTER_MIN,
  USE_ABOVE_THRESHOLD_DEFAULT} from './constants';
import {deepCopy} from '../../../shared-utils/utils';

/** Returns indices corresponding to the closest items */
export function getIndices(expArg: DG.Column, simArg: DG.Column): Uint32Array {
  const expArgRaw = expArg.getRawData();
  const simArgRaw = simArg.getRawData();
  const simCount = simArg.length;
  const elemsCount = Math.min(expArg.length, simCount);
  const indeces = new Uint32Array(elemsCount);
  let idxSim = 0;
  let difPrev = 0;
  let difCur = 0;

  for (let idxExp = 0; idxExp < elemsCount; ++idxExp) {
    while (true) {
      difPrev = Math.abs(simArgRaw[idxSim] - expArgRaw[idxExp]);
      ++idxSim;

      if (idxSim < simCount) {
        difCur = Math.abs(simArgRaw[idxSim] - expArgRaw[idxExp]);

        if (difCur > difPrev) {
          --idxSim;
          break;
        } else
          difPrev = difCur;
      } else {
        --idxSim;
        break;
      }
    }

    indeces[idxExp] = idxSim;
  }

  return indeces;
};

/** Return errors of approximation */
export function getErrors(expArg: DG.Column | null, expFuncs: DG.Column[],
  simDf: DG.DataFrame, toScale: boolean): Float32Array {
  if (expArg === null)
    throw new InconsistentTables('no argument column in the target output dataframe');

  const arg = expArg.name;
  const simArg = simDf.col(arg);

  if (simArg === null)
    throw new InconsistentTables(`no "${arg}" column in the output dataframe "${simDf.name}"`);

  const indices = getIndices(expArg, simArg);
  const expColsCount = expFuncs.length;
  const errors = new Float32Array(expColsCount * indices.length);
  let errIdx = 0;

  for (let idx = 0; idx < expColsCount; ++idx) {
    const expCol = expFuncs[idx];
    const simCol = simDf.col(expCol.name);

    if (simCol === null)
      throw new InconsistentTables(`no "${expCol.name}" column in the output dataframe "${simDf.name}"`);

    const simRaw = simCol.getRawData();
    const expRaw = expCol.getRawData();

    if (toScale) {
      const expScale = Math.max(Math.abs(expCol.stats.max), Math.abs(expCol.stats.min));
      const coef = (expScale > 0) ? expScale : 1;

      indices.forEach((simIdx, expIdx) => {
        errors[errIdx] = (simRaw[simIdx] - expRaw[expIdx]) / coef;
        ++errIdx;
      });
    } else {
      indices.forEach((simIdx, expIdx) => {
        errors[errIdx] = simRaw[simIdx] - expRaw[expIdx];
        ++errIdx;
      });
    }
  }

  return errors;
} // getErrors


/** Get call funcCall with the specified inputs */
export function makeGetCalledFuncCall(func: DG.Func, inputs: Record<string, any>, variedInputNames: string[], useClone: boolean) {
  const funcCall = func.prepare(inputs);
  const resetSharedCall = () => {
    for (const param of funcCall.inputParams.values())
      funcCall.inputs[param.name] = inputs[param.name];
    for (const param of funcCall.outputParams.values())
      funcCall.outputs[param.name] = undefined;
  };

  return async function getCalledFuncCall(x: Float64Array): Promise<DG.FuncCall> {
    resetSharedCall();
    x.forEach((val, idx) => funcCall.inputs[variedInputNames[idx]] = val);
    await funcCall.call(undefined, undefined, {processed: true, report: false});
    return useClone ? deepCopy(funcCall) : funcCall;
  };
} // makeGetCalledFuncCall


/** Convert bounds data to inputs */
export function getInputsData(inputBounds: Record<string, ValueBoundsData>) {
  const variedInputNames: string[] = [];
  const fixedInputs: Record<string, any> = {};
  for (const [name, bound] of Object.entries(inputBounds)) {
    if (bound.type === 'const')
      fixedInputs[name] = bound.value;
    else
      variedInputNames.push(name);
  }
  return {variedInputNames, fixedInputs};
} // getInputsData


/** Return widget for show/hide group of inputs */
export function getCategoryWidget(category: string, roots: HTMLElement[],
  expandHandler?: (r: HTMLElement, isExpanded: boolean, category: string) => void) {
  const updateWgts = (isExpanded: boolean) => {
    chevronToOpen.hidden = isExpanded;
    chevronToClose.hidden = !isExpanded;
    roots.forEach((r) => expandHandler ? expandHandler(r, isExpanded, category) : (r.hidden = !isExpanded));
  };

  const chevronToOpen = ui.iconFA('chevron-right', () => updateWgts(true), 'Open category');
  chevronToOpen.classList.add('fit-view-chevron-right');
  chevronToOpen.hidden = true;

  const chevronToClose = ui.iconFA('chevron-down', () => updateWgts(false), 'Close category');
  chevronToClose.classList.add('fit-view-chevron-down');

  return ui.divH(
    [chevronToOpen, chevronToClose, ui.label(category)],
    'fit-view-inputs-category',
  );
}

function isFullyVisible(element: HTMLElement) {
  const rect = element.getBoundingClientRect();
  return (
    rect.top >= 0 &&
    rect.left >= 0 &&
    rect.bottom <= (window.innerHeight || document.documentElement.clientHeight) &&
    rect.right <= (window.innerWidth || document.documentElement.clientWidth)
  );
}

/** Return the Show Info widget */
export function getShowInfoWidget(root: HTMLElement, dfName: string) {
  let isInfoShown = false;
  let popup: HTMLDivElement;
  let closeIcn: HTMLElement;
  const info = `# **${dfName}**\n\n${TARGET_DATAFRAME_INFO}`;

  const infoIcon = ui.icons.info(() => {
    if (isInfoShown) {
      isInfoShown = false;
      popup!.remove();
    } else {
      isInfoShown = true;
      popup = ui.hints.addHint(root, ui.markdown(info), ui.hints.POSITION.RIGHT);

      if (!isFullyVisible(popup)) {
        popup.remove();
        popup = ui.hints.addHint(root, ui.markdown(info), ui.hints.POSITION.TOP);
      }

      closeIcn = popup.querySelector('i') as HTMLElement;
      closeIcn.onclick = () => isInfoShown = false;
      grok.shell.v.root.appendChild(popup);
    }
  }, 'Click to see details');

  infoIcon.classList.add('sa-switch-input');

  return infoIcon;
} // getShowInfoWidget

/** Return the open help widget */
export function getHelpIcon(): HTMLElement {
  const icon = ui.icons.help(() => window.open(HELP_LINK, '_blank'), 'Open help in a new tab');
  icon.classList.add('fit-view-help-icon');

  return icon;
}

/** Return fitting widget */
export function getFittingWgt(): HTMLElement {
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('width', '24');
  svg.setAttribute('height', '24');
  svg.setAttribute('viewBox', '0 0 24 24');

  const addAttrs = (element: SVGPathElement | SVGCircleElement, attrs: Map<string, string>) => {
    attrs.forEach((value, attr) => {
      element.setAttribute(attr, value);
    });
  };

  const path = document.createElementNS(svgNS, 'path');
  addAttrs(path, new Map([
    ['d', 'M3 4 Q12 28, 20 4'],
    ['stroke', '#40607F'],
    ['stroke-width', '1'],
    ['fill', 'none'],
  ]));

  const circle1 = document.createElementNS(svgNS, 'circle');
  addAttrs(circle1, new Map([
    ['cx', '19'],
    ['cy', '6'],
    ['r', '1.5'],
    ['fill', '#40607F'],
  ]));

  const circle2 = document.createElementNS(svgNS, 'circle');
  addAttrs(circle2, new Map([
    ['cx', '3'],
    ['cy', '10'],
    ['r', '1.5'],
    ['fill', '#40607F'],
  ]));

  const circle3 = document.createElementNS(svgNS, 'circle');
  addAttrs(circle3, new Map([
    ['cx', '10'],
    ['cy', '11'],
    ['r', '1.5'],
    ['fill', '#40607F'],
  ]));

  const circle4 = document.createElementNS(svgNS, 'circle');
  addAttrs(circle4, new Map([
    ['cx', '17'],
    ['cy', '16'],
    ['r', '1.5'],
    ['fill', '#40607F'],
  ]));

  [path, circle1, circle2, circle3, circle4].forEach((element) => svg.appendChild(element));

  const span = ui.span(['Fit']);
  span.classList.add('fit-view-ribbon-text');

  svg.classList.add('fit-view-svg-icon');
  const div = ui.div([svg, span]);

  ui.tooltip.bind(div, 'Fit parameters. Opens a separate view');

  return div;
}

/** Return sensitivity analysis widget */
export function getSensAnWgt(): HTMLElement {
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');
  svg.setAttribute('width', '24');
  svg.setAttribute('height', '24');
  svg.setAttribute('viewBox', '0 0 24 24');

  const path = document.createElementNS(svgNS, 'path');
  // eslint-disable-next-line max-len
  path.setAttribute('d', 'M5 5.5L5 16.5C5 17.3438 5.65625 18 6.5 18H19.5C19.75 18 20 18.25 20 18.5C20 18.7813 19.75 19 19.5 19H6.5C5.09375 19 4 17.9063 4 16.5L4 5.5C4 5.25 4.21875 5 4.5 5C4.75 5 5 5.25 5 5.5ZM15.5 7C15.2188 7 15 6.78125 15 6.5C15 6.25 15.2188 6 15.5 6L18.5 6C18.75 6 19 6.25 19 6.5V9.5C19 9.78125 18.75 10 18.5 10C18.2188 10 18 9.78125 18 9.5V7.71875L13.3438 12.375C13.1562 12.5625 12.8125 12.5625 12.625 12.375L10.5 10.2188L7.84375 12.875C7.65625 13.0625 7.3125 13.0625 7.125 12.875C6.9375 12.6875 6.9375 12.3438 7.125 12.1563L10.125 9.15625C10.2188 9.0625 10.3438 9 10.5 9C10.625 9 10.75 9.0625 10.8438 9.15625L13 11.3125L17.2812 7L15.5 7ZM15.5 16C15.2188 16 15 15.7813 15 15.5C15 15.25 15.2188 15 15.5 15H17.2812L15.125 12.875L15.8438 12.1563L18 14.3125V12.5C18 12.25 18.2188 12 18.5 12C18.75 12 19 12.25 19 12.5V15.5C19 15.7813 18.75 16 18.5 16H15.5Z');
  path.setAttribute('fill', '#40607F');
  svg.appendChild(path);
  const span = ui.span(['Sensitivity']);
  span.classList.add('fit-view-ribbon-text');
  svg.classList.add('sensitivity-analysis-svg-icon');
  const div = ui.div([svg, span]);
  ui.tooltip.bind(div, 'Run sensitivity analysis. Opens a separate view');

  return div;
}

/** Return optimization widget */
export function getOptimizationAnWgt(): HTMLElement {
  const span = ui.span(['Optimize']);
  span.classList.add('fit-view-ribbon-text');
  const icn = ui.iconFA('arrow-up');
  const div = ui.div([icn, span]);
  ui.tooltip.bind(div, 'Run optimization. Opens a separate view');

  return div;
}

/** Return dataframe with loss function vals */
export function getLossFuncDf(extr: Extremum): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.INT, TITLE.ITER, [...Array(extr.iterCount).keys()].map((i) => i + 1)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, TITLE.LOSS, extr.iterCosts.slice(0, extr.iterCount)),
  ]);
}

/** Transform RGB-color to hex number */
export function rgbToHex(rgbString: string) {
  // Extract the RGB values using a regular expression
  const match = rgbString.match(/rgb\((\d+),\s*(\d+),\s*(\d+)\)/);
  if (!match)
    throw new Error('Invalid RGB format. Expected format: \'rgb(r, g, b)\'');


  // Parse the RGB values
  const r = parseInt(match[1], 10);
  const g = parseInt(match[2], 10);
  const b = parseInt(match[3], 10);

  // Convert each component to a 2-digit hexadecimal string
  const toHex = (value: number) => value.toString(16).padStart(2, '0');

  const hexR = toHex(r);
  const hexG = toHex(g);
  const hexB = toHex(b);

  // Return the concatenated hex color string
  return `#${hexR}${hexG}${hexB}`;
}

/** Return RGB-color lightened with respect to the given percentage */
export function lightenRGB(rgbString: string, percentage: number) {
  // Extract the RGB values using a regular expression
  const match = rgbString.match(/rgb\((\d+),\s*(\d+),\s*(\d+)\)/);
  if (!match)
    throw new Error('Invalid RGB format. Expected format: \'rgb(r, g, b)\'');


  // Parse the RGB values
  const r = parseInt(match[1], 10);
  const g = parseInt(match[2], 10);
  const b = parseInt(match[3], 10);

  // Calculate the lighter values
  const lighten = (value: number) => Math.min(255, Math.round(value + (255 - value) * (percentage / 100)));

  const newR = lighten(r);
  const newG = lighten(g);
  const newB = lighten(b);

  // Return the new RGB string
  return `rgb(${newR}, ${newG}, ${newB})`;
}

/** Return options for single scalar goodness of fit bar chart */
function getSingleScalarBarChartOpts(valColName: string, catColName: string): Partial<DG.IBarChartSettings> {
  return {
    showValueAxis: false,
    showStackSelector: false,
    showCategorySelector: false,
    valueColumnName: valColName,
    splitColumnName: catColName,
    colorColumnName: valColName,
    colorAggrType: DG.AGG.AVG,
    linearColorScheme: [DG.Color.categoricalPalette[INDICES.LIGHT_OLIVE], DG.Color.categoricalPalette[INDICES.OLIVE]],
    valueAggrType: DG.AGG.AVG,
    showValueSelector: false,
    barSortType: 'by category',
  };
}

/** Return options for the first chart of double scalar goodness of fit bar chart */
function getFirstDoubleScalarBarChartOpts(valColName: string, catColName: string): Partial<DG.IBarChartSettings> {
  return {
    showValueAxis: false,
    showStackSelector: false,
    showCategorySelector: false,
    valueColumnName: valColName,
    splitColumnName: catColName,
    colorColumnName: valColName,
    colorAggrType: DG.AGG.AVG,
    linearColorScheme: [DG.Color.categoricalPalette[INDICES.LIGHT_OLIVE], DG.Color.categoricalPalette[INDICES.OLIVE]],
    valueAggrType: DG.AGG.AVG,
    showValueSelector: false,
    description: valColName,
    descriptionVisibilityMode: 'Always',
    maxBarHeight: GRID_SIZE.BAR_HEIGHT,
    barSortType: 'by category',
  };
}

/** Return options for the second chart of double scalar goodness of fit bar chart */
function getSecondDoubleScalarBarChartOpts(valColName: string, catColName: string): Partial<DG.IBarChartSettings> {
  return {
    showValueAxis: false,
    showStackSelector: false,
    showCategorySelector: false,
    valueColumnName: valColName,
    splitColumnName: catColName,
    colorColumnName: valColName,
    colorAggrType: DG.AGG.AVG,
    linearColorScheme: [DG.Color.categoricalPalette[INDICES.GREY], DG.Color.categoricalPalette[INDICES.LIGHT_GREY]],
    valueAggrType: DG.AGG.AVG,
    showValueSelector: false,
    description: valColName,
    descriptionVisibilityMode: 'Always',
    maxBarHeight: GRID_SIZE.BAR_HEIGHT,
    barSortType: 'by category',
  };
}

/** Return goodness of fit viewer for scalar targets */
export function getScalarsGoodnessOfFitViewer(table: DG.DataFrame): HTMLElement {
  const cols = table.columns;
  const colsCount = cols.length;

  if (colsCount < 2)
    throw new Error(`Incorrect scalars fit table. Columns: ${colsCount}, expected: > 1`);

  switch (colsCount) {
  case 2:
    return DG.Viewer.barChart(table, getSingleScalarBarChartOpts(
      cols.byIndex(INDICES.SINGLE_VAL).name,
      cols.byIndex(INDICES.SINGLE_CAT).name,
    )).root;

  case 3:
    const catColName = cols.byIndex(INDICES.DOUBLE_CAT).name;

    const first = DG.Viewer.barChart(table, getFirstDoubleScalarBarChartOpts(
      cols.byIndex(INDICES.DOUBLE_VAL1).name,
      catColName,
    )).root;

    const second = DG.Viewer.barChart(table, getSecondDoubleScalarBarChartOpts(
      cols.byIndex(INDICES.DOUBLE_VAL2).name,
      catColName,
    )).root;

    const container = ui.splitV([first, second], {style: {height: `${GRID_SIZE.ROW_HEIGHT}px`}});
    setTimeout(() => container.style.height = '100%', TIMEOUT.STYLE_TIMEOUT);

    return container;

  default:
    return DG.Viewer.grid(table).root;
  }
} // getScalarsGoodnessOfFitViewer

/** Return radar tooltip */
export function getRadarTooltip(): HTMLElement {
  const colors = DG.Color.categoricalPalette;
  const simLabel = ui.label(NAME.SIMULATION);
  simLabel.style.color = rgbToHex(DG.Color.toRgb(colors[INDICES.BLUE]));
  const targetLabel = ui.label(NAME.TARGET);
  targetLabel.style.color = rgbToHex(DG.Color.toRgb(colors[INDICES.ORANGE]));

  return ui.divV([simLabel, targetLabel]);
}

/** Return true iff the radar viewer can be used */
export function toUseRadar(table: DG.DataFrame): boolean {
  // check that table contains results on scalar outputs <=> it has the Category column
  if (table.col(NAME.CATEGORY) == null)
    return false;

  return (table.columns.length >= MIN_RADAR_COLS_COUNT);
}

/** Return seeded random generator */
export function seededRandom(seed: number): () => number {
  let state = seed % RAND_MOD;
  if (state <= 0) state += RAND_MOD - 1;

  return function(): number {
    state = (state * RAND_MULT) % RAND_MOD;
    return (state - 1) / (RAND_MOD - 1);
  };
}

export const defaultRandomSeedSettings: ReproSettings = {
  reproducible: REPRO_DEFAULT,
  seed: SEED_DEFAULT,
};

/** Returns random seed settings */
export function getRandomSeedSettings(overrides: Record<string, any> = {}) {
  const settings = {...defaultRandomSeedSettings, ...overrides};

  const reprInput = ui.input.bool('reproducible', {
    value: settings.reproducible,
    onValueChanged: (val) => {
      seedInput.root.hidden = !val;
      settings.reproducible = val;
    },
    tooltipText: 'Enable to get reproducible results',
  });

  const seedInput = ui.input.int('random seed', {
    value: settings.seed,
    nullable: false,
    tooltipText: 'Numeric value used to initialize the initial points of the optimization method',
    onValueChanged: (val) => {
      settings.seed = val ?? SEED_DEFAULT;
    },
  });

  seedInput.root.hidden = !settings.reproducible;

  return {
    reproducibility: reprInput,
    seed: seedInput,
    settings: settings,
  };
} // getRandomSeedSettings

export const defaultEarlyStoppingSettings: EarlyStoppingSettings = {
  useEarlyStopping: EARLY_STOP_DEFAULT,
  costFuncThreshold: COST_FUNC_THRESH,
  stopAfter: STOP_AFTER_DEFAULT,
  useAboveThresholdPoints: USE_ABOVE_THRESHOLD_DEFAULT,
};

/** Returns early stopping intpus & settings */
export function getEarlyStoppingInputs(overrides: Record<string, any> = {}) {
  const settings = {...defaultEarlyStoppingSettings, ...overrides};

  const stopAtFirstInput = ui.input.int('stop after', {
    value: settings.stopAfter,
    onValueChanged: (val) => {
      settings.stopAfter = val;
    },
    tooltipText: 'Stop the optimization once this number of valid points is found',
    showPlusMinus: true,
    min: STOP_AFTER_MIN,
    nullable: false,
    property: DG.Property.fromOptions({units: 'point(s)'}),
  });
  stopAtFirstInput.root.hidden = !settings.useEarlyStopping;

  const earlyStopInput = ui.input.bool('early stopping', {
    value: settings.useEarlyStopping,
    onValueChanged: (val) => {
      thresholdInput.root.hidden = !val;
      stopAtFirstInput.root.hidden = !val;
      settings.useEarlyStopping = val;
    },
    tooltipText: 'Enable to stop fitting once the loss function reaches the specified threshold',
  });

  const thresholdInput = ui.input.float('threshold', {
    value: settings.costFuncThreshold,
    nullable: false,
    tooltipText: 'Loss function value at which to stop fitting',
    onValueChanged: (val) => {
      settings.costFuncThreshold = val ?? COST_FUNC_THRESH;
    },
  });

  thresholdInput.root.hidden = !settings.useEarlyStopping;

  return {
    stopAtFirst: stopAtFirstInput,
    earlyStopping: earlyStopInput,
    threshold: thresholdInput,
    settings: settings,
  };
} // getEarlyStoppingInputs
