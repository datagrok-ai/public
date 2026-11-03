import * as DG from 'datagrok-api/dg';
import {waitUntil} from '@datagrok-libraries/test/src/test';

export function getFuncCallIO(fc: DG.FuncCall) {
  return {
    inputs: [...fc.inputs.entries()],
    outputs: [...fc.outputs.entries()],
  };
}

/** Resolves once the view has rendered an input for each input of its function. */
export async function inputsRendered(view: {funcCall: DG.FuncCall}, inputsMap: Record<string, unknown>): Promise<void> {
  await waitUntil(() => view.funcCall.func.inputs.every((p) => inputsMap[p.name] != null),
    () => `inputs not rendered: ${view.funcCall.func.inputs.map((p) => p.name).filter((n) => inputsMap[n] == null)}`);
}

/** Resolves once the view's call holds [values]: an input passes its value on after a validation debounce. */
export async function inputsApplied(view: {funcCall: DG.FuncCall}, values: Record<string, any>): Promise<void> {
  // a dataframe comes back in a new wrapper on every read
  const holds = (name: string) => values[name] instanceof DG.DataFrame ?
    view.funcCall.inputs[name]?.dart === values[name].dart : view.funcCall.inputs[name] === values[name];
  await waitUntil(() => Object.keys(values).every(holds),
    () => `the call did not take the inputs: ${JSON.stringify(Object.keys(values).filter((name) => !holds(name)))}`);
}

/** Retries [assertion] until it passes; after [timeout] ms throws what it threw last. */
export async function eventually(assertion: () => void, timeout: number = 10000): Promise<void> {
  let last: unknown = null;
  const passes = () => {
    try {
      assertion();
      return true;
    } catch (e) {
      last = e;
      return false;
    }
  };
  await waitUntil(passes, 'unreached', timeout).catch(() => {
    throw last;
  });
}
