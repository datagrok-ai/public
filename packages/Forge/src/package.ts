/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import {ForgeApp} from './ui/forge-app';
export * from './package.g';

export const _package = new DG.Package();

//name: Forge
//description: Predictive modeling: engines and the model catalog
//tags: app
//output: view result
export async function forgeApp(): Promise<DG.ViewBase> {
  return ForgeApp.create();
}

//name: forgeModels
//description: Opens the Forge app
//top-menu: ML | Forge | Models
export async function forgeModels(): Promise<void> {
  await ForgeApp.open();
}
