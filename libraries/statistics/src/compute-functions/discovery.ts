import * as DG from 'datagrok-api/dg';
import {ComputeFunctions} from './types';
import {funcTypeNames} from './consts';

export function discoverComputeFunctions(tag: string): ComputeFunctions {
  const funcs = Array.from(new Set(
    DG.Func.find({meta: {role: tag}}).concat(DG.Func.find({tags: [tag]})),
  ));
  const functions = funcs.filter((f) => f.type === funcTypeNames.function);
  const scripts = funcs.filter((f) => f.type === funcTypeNames.script) as DG.Script[];
  const queries = funcs.filter((f) => f.type === funcTypeNames.query) as DG.DataQuery[];
  return {functions, scripts, queries};
}
