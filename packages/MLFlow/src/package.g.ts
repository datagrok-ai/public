import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: MLflow models
//description: Lists the MLflow models synced from the MLFlow connections: id, name, target, features, created
//output: dataframe result
//meta.mlname: MLflow
//meta.mlrole: models
export async function mlflowModels() : Promise<any> {
  return await PackageFunctions.mlflowModels();
}

//name: Apply MLflow
//description: Predicts with an MLflow model in its container, which the first application starts
//input: dataframe df { description: The model's input columns, in its order. }
//input: dynamic model { description: The model id from MLflow models, as UTF-8 bytes. }
//output: dataframe result
//meta.mlname: MLflow
//meta.mlrole: apply
export async function applyMlflow(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyMlflow(df, model);
}
