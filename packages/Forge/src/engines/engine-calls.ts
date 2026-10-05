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

/** The first column of the engine's result, taken out of the result frame, so a table it is added to is its parent. */
export async function apply(engine: Engine, features: DG.DataFrame, blob: Uint8Array): Promise<DG.Column> {
  const func = roleFunc(engine, 'apply', 'cannot apply models');
  const result: unknown = await func.apply({df: features, model: blob});
  if (!(result instanceof DG.DataFrame) || result.columns.length === 0)
    throw new ForgeError(`The method '${engine.name}' returned no predictions.`);
  const prediction = result.columns.byIndex(0);
  result.columns.remove(prediction, false);
  return prediction;
}

/** The part of a progress indicator a loop of engine calls uses: the cancel flag and the update. */
export type LoopProgress = Pick<DG.ProgressIndicator, 'canceled' | 'update'>;

/** Lets the browser handle pending events, such as a click on a progress's cancel, between engine calls: awaited
 * engine calls that finish without I/O resume as microtasks, so a loop of them never reaches the event loop.
 * A message, unlike `setTimeout`, is not throttled to about a second in a background tab. */
export function yieldToEventLoop(): Promise<void> {
  return new Promise((resolve) => {
    const channel = new MessageChannel();
    channel.port1.onmessage = () => {
      channel.port1.close();
      resolve();
    };
    channel.port2.postMessage(null);
  });
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
