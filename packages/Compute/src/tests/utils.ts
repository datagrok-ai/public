import * as DG from 'datagrok-api/dg';

export function getFuncCallIO(fc: DG.FuncCall) {
  return {
    inputs: [...fc.inputs.entries()],
    outputs: [...fc.outputs.entries()],
  };
}
