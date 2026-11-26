import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Sync MLflow models
//description: Writes the MLflow models you may use into the Forge catalog; returns what was not refreshed
//output: string result
//meta.mlname: MLflow
//meta.mlrole: sync
export async function syncMlflow() : Promise<string> {
  return await PackageFunctions.syncMlflow();
}

//name: Apply MLflow
//description: Predicts with a published MLflow model in its container, which the first application starts
//input: dataframe df { description: The model's input columns, in its order. }
//input: dynamic model { description: The model reference as UTF-8 json bytes. }
//output: dataframe result
//meta.mlname: MLflow
//meta.mlrole: apply
export async function applyMlflow(df: DG.DataFrame, model: Uint8Array) : Promise<any> {
  return await PackageFunctions.applyMlflow(df, model);
}
