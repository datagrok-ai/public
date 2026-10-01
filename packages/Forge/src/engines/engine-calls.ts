import * as DG from 'datagrok-api/dg';
import {Engine} from './engine';
import {ForgeError} from '../forge-error';

export async function isApplicable(engine: Engine, features: DG.DataFrame, target: DG.Column): Promise<boolean> {
  const func = engine.functions.isApplicable;
  if (func === undefined)
    throw new ForgeError(`The engine '${engine.name}' does not say which data it supports. Choose another engine.`);
  return await callCheck(engine, func, features, target);
}

export async function isInteractive(engine: Engine, features: DG.DataFrame, target: DG.Column): Promise<boolean> {
  const func = engine.functions.isInteractive;
  return func === undefined ? false : await callCheck(engine, func, features, target);
}

async function callCheck(engine: Engine, func: DG.Func, features: DG.DataFrame, target: DG.Column): Promise<boolean> {
  if (engine.kind === 'function') {
    const result: unknown = await func.apply([features, target]);
    return result === true;
  }
  const table = features.clone();
  table.columns.add(target.clone());
  const result: unknown = await func.apply([table, target.name]);
  return result === true;
}
