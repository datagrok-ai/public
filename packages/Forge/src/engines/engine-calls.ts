import * as DG from 'datagrok-api/dg';
import {Engine, EngineRole, Hyperparameters} from './engine';
import {ForgeError} from '../forge-error';

export async function isApplicable(engine: Engine, features: DG.DataFrame, target: DG.Column): Promise<boolean> {
  const func = roleFunc(engine, 'isApplicable', 'does not say which data it supports');
  return await callCheck(engine, func, features, target);
}

export async function isInteractive(engine: Engine, features: DG.DataFrame, target: DG.Column): Promise<boolean> {
  const func = engine.functions.isInteractive;
  return func === undefined ? false : await callCheck(engine, func, features, target);
}

export async function train(engine: Engine, features: DG.DataFrame, target: DG.Column,
  hyperparameters: Hyperparameters): Promise<Uint8Array> {
  const func = roleFunc(engine, 'train', 'cannot train models');
  const [df, predictColumn] = dataArgs(engine, features, target);
  const result: unknown = await func.apply({...hyperparameters, df, predictColumn});
  if (result instanceof Uint8Array)
    return result;
  if (result instanceof DG.FileInfo)
    return result.data;
  throw new ForgeError(`The method '${engine.name}' returned no model.`);
}

export async function apply(engine: Engine, features: DG.DataFrame, blob: Uint8Array): Promise<DG.Column> {
  const func = roleFunc(engine, 'apply', 'cannot apply models');
  const result: unknown = await func.apply({df: features, model: blob});
  if (result instanceof DG.DataFrame && result.columns.length > 0)
    return result.columns.byIndex(0);
  throw new ForgeError(`The method '${engine.name}' returned no predictions.`);
}

function roleFunc(engine: Engine, role: EngineRole, problem: string): DG.Func {
  const func = engine.functions[role];
  if (func === undefined)
    throw new ForgeError(`The method '${engine.name}' ${problem}. Choose another method.`);
  return func;
}

async function callCheck(engine: Engine, func: DG.Func, features: DG.DataFrame, target: DG.Column): Promise<boolean> {
  const result: unknown = await func.apply(dataArgs(engine, features, target));
  return result === true;
}

/** Function engines take the feature table and the target column; script engines one table and the target name. */
function dataArgs(engine: Engine, features: DG.DataFrame, target: DG.Column): [DG.DataFrame, DG.Column | string] {
  if (engine.kind === 'function')
    return [features, target];
  const table = features.clone();
  table.columns.add(target.clone());
  return [table, target.name];
}
