/* eslint-disable valid-jsdoc */
// Fitting utilities

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {InconsistentTables} from './optimizer-misc';

import '../../css/fitting-view.css';
import '../../css/sens-analysis.css';
import {TARGET_DATAFRAME_INFO} from './constants';

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

/** Return widget for show/hide group of inputs */
export function getCategoryWidget(category: string, roots: HTMLElement[]) {
  const updateWgts = (isExpanded: boolean) => {
    chevronToOpen.hidden = isExpanded;
    chevronToClose.hidden = !isExpanded;

    roots.forEach((r) => r.hidden = !isExpanded);
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

  const div = ui.div([svg]);
  div.classList.add('fit-view-svg-icon');
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

  const div = ui.div([svg]);
  div.classList.add('sensitivity-analysis-svg-icon');
  ui.tooltip.bind(div, 'Run sensitivity analysis. Opens a separate view');

  return div;
}
