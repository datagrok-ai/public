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

function moveElementIntoView(element: HTMLElement) {
  if (!element) return;

  const rect = element.getBoundingClientRect();
  const offsetX = Math.max(0, -rect.left) + Math.max(0, rect.right - window.innerWidth);
  const offsetY = Math.max(0, -rect.top) + Math.max(0, rect.bottom - window.innerHeight);

  // Apply translation via CSS transform
  element.style.position = 'absolute';
  element.style.transform = `translate(${offsetX}px, ${offsetY}px)`;
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
