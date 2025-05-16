/* eslint-disable valid-jsdoc */
// Optimization specific features
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../../css/optimization-view.css';

export const HELP_URL = '/help/compute/function-analysis#optimization';
export const HELP_FULL_URL = `https://datagrok.ai${HELP_URL}`;

export const STARTING_HELP = `# Optimization

Use optimization to find ...`;

export enum OPT_TYPE {
  MIN = 'Minimize',
  MAX = 'Maximize',
};

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
  SUBJ = 'Subject to',
  METHOD = 'Method',
};

export enum METHOD {
  NELDER_MEAD = 'Nelder-Mead',
  MOEAD = 'MOEA/D',
}

/** Return the open help widget */
export function getHelpIcon(): HTMLElement {
  const icon = ui.icons.help(() => window.open(HELP_FULL_URL, '_blank'), 'Open help in a new tab');
  icon.classList.add('fit-view-help-icon');

  return icon;
}
