import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DIRECTION} from './consts';

// DEALING WITH BUG: https://reddata.atlassian.net/browse/GROK-13015
export const getDataFrame = (call: DG.FuncCall, name: string, direction: DIRECTION): DG.DataFrame => {
  if (direction === DIRECTION.INPUT)
    return call.inputs[name] instanceof DG.DataFrame ? call.inputs[name] : call.inputs[name].dataFrame;
  else
    return call.outputs[name] instanceof DG.DataFrame ? call.outputs[name] : call.outputs[name].dataFrame;
};

export const boundImportFunction = (func: DG.Func): string | undefined => {
  return func.options['getRealData'];
};
