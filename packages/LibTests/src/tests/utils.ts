import * as DG from 'datagrok-api/dg';
import cloneDeepWith from 'lodash.clonedeepwith';

export function getFuncCallIO(fc: DG.FuncCall) {
  return {
    inputs: [...fc.inputs.entries()],
    outputs: [...fc.outputs.entries()],
  };
}
