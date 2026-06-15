import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {DF_NAME} from './constants';
import {getIVP, getScriptLines} from './scripting-tools';

/** Solve a Diff Studio `.ivp` model (referenced by its App Data path) with the given input
 *  values and return the result dataframe. Used by the Model-Hub model wrappers that the
 *  func-gen build plugin generates from `#meta.role: model` `.ivp` files. */
export async function runDiffStudioModel(
  modelPath: string, inputs: Record<string, number>): Promise<DG.DataFrame> {
  const text = await grok.dapi.files.readAsText(modelPath);
  const ivp = getIVP(text);
  const script = DG.Script.create(getScriptLines(ivp).join('\n'));
  const call = script.prepare(inputs);
  await call.call();
  return call.outputs[DF_NAME];
}
