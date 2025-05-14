/* eslint-disable valid-jsdoc */
// utils.ts

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {HELP_LINK} from './constants';

import '../../css/sens-analysis.css';

// Error messeges
enum ERROR_MSG {
  INCORRECT_SIZE_TYPE = 'Size must be an integer.',
  SIZE_MUST_BE_POSITIVE = 'Size must be positive.',
};

enum FUNCCALL_CONSTS {
  BATCH_SIZE = 100,
  MS_TO_SLEEP = 10,
};

const sleep = (ms: number) => {
  return new Promise((resolve, reject) => setTimeout(resolve, ms));
};

// Check variable that is supposed to be a size
export function checkSize(size: any): void {
  // check if it's a number
  if (!Number.isInteger(size))
    throw new Error(ERROR_MSG.INCORRECT_SIZE_TYPE);

  // check positivity
  if (size < 1)
    throw new Error(ERROR_MSG.SIZE_MUST_BE_POSITIVE);
}

export async function getCalledFuncCalls(funccalls: DG.FuncCall[]): Promise<DG.FuncCall[]> {
  const calledFuncCalls = [] as DG.FuncCall[];
  const pi = DG.TaskBarProgressIndicator.create(`Analyzing...`);

  // the feature cancelable should be added
  (pi as any).cancelable = true;

  let idx = 0;
  const funcCallsCount = funccalls.length;
  let percentage = 0;

  for (let i = 0; i < funcCallsCount; ++i) {
    calledFuncCalls.push(await funccalls[i].call());

    if (i === 1) {
      percentage = Math.floor(100 * i / funcCallsCount);
      pi.update(percentage, `Analyzing... (${percentage}%)`);
      await sleep(FUNCCALL_CONSTS.MS_TO_SLEEP * 10);
    }

    ++idx;

    // computation cancel check
    if ((pi as any).canceled)
      break;

    if (idx >= FUNCCALL_CONSTS.BATCH_SIZE) {
      percentage = Math.floor(100 * i / funcCallsCount);
      pi.update(percentage, `Analyzing... (${percentage}%)`);
      await sleep(FUNCCALL_CONSTS.MS_TO_SLEEP);
      idx = 0;
    }
  }

  pi.close();

  return calledFuncCalls;
}

/** Return the open help widget */
export function getHelpIcon(): HTMLElement {
  const icon = ui.icons.help(() => window.open(HELP_LINK, '_blank'));
  icon.classList.add('sensitivity-analysis-help-icon');

  return icon;
}
