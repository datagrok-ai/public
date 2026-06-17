/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {getIVP, getIvp2WebWorker, getPipelineCreator} from 'diff-grok';
import {DiffGrok} from '../../fitting-view';

/** Function option key holding the path to the model's `.ivp` text file. */
export const DIFF_STUDIO_MODEL_KEY = 'diffStudioModel';

/** Build the web-worker fitting/SA payload for a Diff Studio model from the `.ivp`
 *  file linked on the function via `meta.diffStudioModel`. Returns `undefined` for
 *  non-IVP functions so callers fall back to the generic (serial) path. */
export async function buildDiffGrokFromFunc(func: DG.Func): Promise<DiffGrok | undefined> {
  const path = func.options[DIFF_STUDIO_MODEL_KEY];
  if (!path)
    return undefined;

  const ivpText = await grok.dapi.files.readAsText(path);
  const ivp = getIVP(ivpText);
  return {ivp, ivpWW: getIvp2WebWorker(ivp), pipelineCreator: getPipelineCreator(ivp)};
} // buildDiffGrokFromFunc
