/* eslint-disable valid-jsdoc */

import {IVP, IVP2WebWorker} from '@datagrok/diff-studio-tools';
import {LOSS} from '../constants';

/** Return fitted params of Diff Studio model using the Nelder-Mead method */
export async function getFittedParams(
  loss: LOSS,
  ivp: IVP,
  ivp2ww: IVP2WebWorker,
  settingsNames: string[],
  settingVals: number[],
  variedInputNames: string[],
  minVals: Float32Array,
  maxVals: Float32Array,
  fixedInputNames: string[],
  fixedVals: Float32Array,
  outputNames: string[],
  outputVals: Float32Array[]): Promise<Float32Array[]> {
  const res: Float32Array[] = [];

  return res;
}
