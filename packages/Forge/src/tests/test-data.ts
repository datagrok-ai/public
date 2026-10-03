import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {defaultHyperparameters, Engine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {ModelInsert} from '../generated/db';
import {TrainingRequest} from '../training/train-model';

export const IRIS = 'System:DemoFiles/iris.csv';
export const MEASUREMENTS = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];
export const XGBOOST_FIELDS: Pick<ModelInsert, 'engine_name' | 'engine_namespace' | 'engine_kind'> =
  {engine_name: 'XGBoost', engine_namespace: 'Eda', engine_kind: 'function'};

export function engineByName(engines: Engine[], name: string): Engine {
  const engine = engines.find((e) => e.name === name);
  if (engine === undefined)
    throw new Error(`Engine '${name}' is not discovered`);
  return engine;
}

let xgboostEngine: Engine | undefined;

export function xgboost(): Engine {
  xgboostEngine ??= engineByName(EngineRegistry.discover(), 'XGBoost');
  return xgboostEngine;
}

export async function openIris(): Promise<DG.DataFrame> {
  const iris = await grok.data.files.openTable(IRIS);
  iris.name = `forge-test-iris-${Date.now()}`;
  return iris;
}

export function requestOf(features: DG.Column[] | DG.DataFrame, target: DG.Column): TrainingRequest {
  const engine = xgboost();
  const frame = features instanceof DG.DataFrame ? features : DG.DataFrame.fromColumns(features.map((c) => c.clone()));
  return {engine, features: frame, target, hyperparameters: defaultHyperparameters(engine), seed: 42, folds: 5};
}
