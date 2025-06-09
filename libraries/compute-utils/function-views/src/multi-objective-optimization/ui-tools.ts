/* eslint-disable valid-jsdoc */
// Optimization specific features
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../../css/optimization-view.css';
import {OPT_TYPE} from './defs';

export const HELP_URL = '/help/compute/function-analysis#optimization';
export const HELP_FULL_URL = `https://datagrok.ai${HELP_URL}`;

export const STARTING_HELP = `# Optimization

Use optimization to find ...`;

export function getOptTypeInput(): DG.ChoiceInput<OPT_TYPE | null> {
  const input = ui.input.choice<OPT_TYPE>('', {
    value: OPT_TYPE.MAX,
    nullable: false,
    items: [OPT_TYPE.MIN, OPT_TYPE.MAX],
    tooltipText: 'Type of optimization',
  });

  input.input.style.border = 'white';
  input.input.classList.add('optimization-type-input');

  return input;
}

export enum TITLE {
  CONSTR = 'Inputs constraints',
  METHOD = 'Method',
  USING = 'Using',
  SET = 'with settings',
  PARETO_FRONT = 'Pareto front',
};

export enum METHOD {
  NELDER_MEAD = 'Nelder-Mead',
  MOEAD = 'MOEA/D',
};

export const METHOD_HINT = new Map([
  [METHOD.MOEAD, 'Multi-objective evolutionary algorithm based on decomposition'],
  [METHOD.NELDER_MEAD, 'The Nelder-Mead method'],
]);

export const INPUT_COLOR_PALETTE = [DG.Color.lightLightGray, DG.Color.lightGray, DG.Color.gray];

export function getOutputPalette(type: OPT_TYPE): number[] {
  const palette = [DG.Color.darkGreen, DG.Color.yellow, DG.Color.darkRed];

  if (type === OPT_TYPE.MIN)
    return palette;

  return palette.reverse();
}

/** Return the open help widget */
export function getHelpIcon(): HTMLElement {
  const icon = ui.icons.help(() => window.open(HELP_FULL_URL, '_blank'), 'Open help in a new tab');
  icon.classList.add('optimization-view-help-icon');

  return icon;
}

/** Return color scale element */
export function getColorScaleDiv(type: OPT_TYPE): HTMLElement {
  const minLbl = ui.label('min');
  const midLbl = ui.label('. . .');
  const maxLbl = ui.label('max');
  const palette = getOutputPalette(type);

  const elems = [minLbl, midLbl, maxLbl].map((el, idx) => {
    if (idx !== 1) {
      el.style.fontWeight = 'bold';
      el.style.color = DG.Color.toRgb(palette[idx]);
    }

    el.style.marginRight = '5px';

    return el;
  });

  return ui.divH(elems);
}

export enum DOCK_RATIO {
  SINGLE_SINGLE_SCATTER = 0.3,
  MULT_SINGLE_PC = 0.5,
};
